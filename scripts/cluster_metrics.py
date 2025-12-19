#!/usr/bin/env python3

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score, calinski_harabasz_score, davies_bouldin_score


def load_features(path: str) -> tuple[pd.DataFrame, pd.Index]:
    df = pd.read_csv(path, sep="\t")
    if "sample_id" not in df.columns:
        raise ValueError("features TSV must have a 'sample_id' column")
    sample_ids = df["sample_id"].astype(str)
    X = df.drop(columns=["sample_id"]).apply(pd.to_numeric, errors="coerce")
    X = X.fillna(X.mean(axis=0))
    return X, sample_ids


def load_clusters(path: str) -> pd.Series:
    df = pd.read_csv(path, sep="\t")
    if "sample_id" not in df.columns or "cluster" not in df.columns:
        raise ValueError("clusters TSV must have columns: sample_id, cluster")
    return df.set_index(df["sample_id"].astype(str))["cluster"].astype(int)


def safe_cluster_metrics(X: np.ndarray, labels: np.ndarray) -> dict:
    # metrics are undefined for 0/1 cluster, or if every point is its own cluster
    uniq = np.unique(labels)
    n_clusters = len(uniq) - (1 if -1 in uniq else 0)
    if n_clusters < 2:
        return {
            "n_clusters": int(n_clusters),
            "silhouette": None,
            "calinski_harabasz": None,
            "davies_bouldin": None,
        }

    # For DBSCAN, ignore noise (-1) for CH/DB; silhouette can be computed on non-noise too
    mask = labels != -1
    X_use = X[mask]
    y_use = labels[mask]
    if len(np.unique(y_use)) < 2:
        return {
            "n_clusters": int(n_clusters),
            "silhouette": None,
            "calinski_harabasz": None,
            "davies_bouldin": None,
        }

    return {
        "n_clusters": int(n_clusters),
        "silhouette": float(silhouette_score(X_use, y_use)),
        "calinski_harabasz": float(calinski_harabasz_score(X_use, y_use)),
        "davies_bouldin": float(davies_bouldin_score(X_use, y_use)),
    }


def main() -> None:
    ap = argparse.ArgumentParser(description="Compute clustering metrics and a KMeans k-sweep.")
    ap.add_argument("--features", required=True)
    ap.add_argument("--clusters", required=True)
    ap.add_argument("--k-min", type=int, default=1)
    ap.add_argument("--k-max", type=int, default=12)
    ap.add_argument("--out-k-sweep", required=True)
    ap.add_argument("--out-selected", required=True)
    args = ap.parse_args()

    X_df, sample_ids = load_features(args.features)
    clusters = load_clusters(args.clusters)

    # Align
    common = sample_ids[sample_ids.isin(clusters.index)]
    if len(common) == 0:
        raise ValueError("No overlapping sample_id between features and clusters")

    X = X_df.loc[common.index].values
    labels = clusters.loc[common.values].values

    # Selected metrics (based on provided labels)
    selected = safe_cluster_metrics(X, labels)
    selected["input_clusters"] = Path(args.clusters).name
    selected["input_features"] = Path(args.features).name

    # KMeans sweep
    kmin = int(args.k_min)
    kmax = int(args.k_max)
    if kmin < 1 or kmax < kmin:
        raise ValueError("Invalid k range")

    rows = []
    for k in range(kmin, kmax + 1):
        model = KMeans(n_clusters=k, n_init="auto", random_state=42)
        y = model.fit_predict(X)
        inertia = float(model.inertia_)
        sil = None
        ch = None
        db = None
        if k >= 2:
            sil = float(silhouette_score(X, y))
            ch = float(calinski_harabasz_score(X, y))
            db = float(davies_bouldin_score(X, y))
        rows.append({
            "k": k,
            "inertia": inertia,
            "silhouette": sil,
            "calinski_harabasz": ch,
            "davies_bouldin": db,
        })

    
    df = pd.DataFrame(rows)
    df.to_csv(args.out_k_sweep, sep=",", index=False)

    # Plots (3): elbow (inertia), silhouette, Davies-Bouldin
    try:
        import matplotlib.pyplot as plt

        def plot_curve(metric: str, title: str, ylabel: str, out_png: str) -> None:
            plt.figure(figsize=(7, 4.5))
            plt.plot(df["k"], df[metric], marker="o")
            plt.xticks(df["k"].tolist())
            plt.title(title)
            plt.xlabel("k")
            plt.ylabel(ylabel)
            plt.tight_layout()
            plt.savefig(out_png, dpi=200)
            plt.close()

        plot_curve("inertia", "Elbow method (KMeans inertia)", "inertia", "elbow.png")
        plot_curve("silhouette", "Silhouette score (higher is better)", "silhouette", "silhouette.png")
        plot_curve("davies_bouldin", "Davies-Bouldin index (lower is better)", "davies_bouldin", "davies_bouldin.png")
        plot_curve(
    "calinski_harabasz",
    "Calinski-Harabasz index (higher is better)",
    "calinski_harabasz",
    "calinski.png",
)

    except Exception as e:
        # If matplotlib isn't available or plotting fails, keep pipeline usable.
        Path("plot_warning.txt").write_text(f"Plotting failed: {e}\\n")

    Path(args.out_selected).write_text(json.dumps(selected, indent=2))


if __name__ == "__main__":
    main()
