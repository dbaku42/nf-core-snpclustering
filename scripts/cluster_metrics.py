#!/usr/bin/env python3
import argparse
import os

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score, calinski_harabasz_score, pairwise_distances


def load_features_and_clusters(features_path, clusters_path):
    if not os.path.exists(features_path):
        raise SystemExit(f"ERROR: features file not found: {features_path}")
    if not os.path.exists(clusters_path):
        raise SystemExit(f"ERROR: clusters file not found: {clusters_path}")

    df = pd.read_csv(features_path, sep="\t")
    if df.shape[1] < 2:
        raise SystemExit("ERROR: matrix must have at least one feature column")

    X_raw = df.iloc[:, 1:]
    X = X_raw.apply(pd.to_numeric, errors="coerce").fillna(-1.0).to_numpy()

    cl = pd.read_csv(clusters_path)
    if "cluster" not in cl.columns:
        raise SystemExit("ERROR: clusters CSV must have a 'cluster' column")

    labels = cl["cluster"].to_numpy()

    n = min(X.shape[0], labels.shape[0])
    if X.shape[0] != labels.shape[0]:
        print(
            f"WARNING: X has {X.shape[0]} rows, labels has {labels.shape[0]} rows; truncating to {n}."
        )
    X = X[:n, :]
    labels = labels[:n]

    return X, labels


def dunn_index_from_distance(D, labels):
    """Dunn index using full distance matrix D (n x n)."""
    unique = np.unique(labels)
    unique = [u for u in unique if u != -1]  # ignore DBSCAN noise if present

    if len(unique) < 2:
        return np.nan

    max_intra = 0.0
    for u in unique:
        idx = np.where(labels == u)[0]
        if len(idx) < 2:
            continue
        sub = D[np.ix_(idx, idx)]
        diam = np.max(sub)
        if diam > max_intra:
            max_intra = diam

    if max_intra == 0:
        return np.nan

    min_inter = np.inf
    for i, u in enumerate(unique):
        idx_i = np.where(labels == u)[0]
        for v in unique[i + 1:]:
            idx_j = np.where(labels == v)[0]
            sub = D[np.ix_(idx_i, idx_j)]
            dmin = np.min(sub)
            if dmin < min_inter:
                min_inter = dmin

    if not np.isfinite(min_inter):
        return np.nan

    return float(min_inter / max_intra)


def plot_metrics_bar(metrics_dict, out_png):
    names = list(metrics_dict.keys())
    values = [metrics_dict[k] for k in names]

    cmap = plt.get_cmap("prism")
    colors = [cmap(i / max(len(names) - 1, 1)) for i in range(len(names))]

    plt.figure(figsize=(5, 4))
    plt.bar(names, values, color=colors)
    plt.ylabel("value")
    plt.title("Clustering metrics (chosen k)")
    plt.xticks(rotation=30, ha="right")
    plt.tight_layout()
    plt.savefig(out_png, dpi=150)
    plt.close()


def plot_metric_vs_k(k_vals, metric_vals, metric_name, out_png):
    cmap = plt.get_cmap("prism")
    color = cmap(0.6)

    plt.figure(figsize=(5, 4))
    plt.plot(k_vals, metric_vals, marker="o", color=color)
    plt.xlabel("k")
    plt.ylabel(metric_name)
    plt.title(f"{metric_name} vs k")
    plt.grid(True, linestyle="--", alpha=0.3)
    plt.tight_layout()
    plt.savefig(out_png, dpi=150)
    plt.close()


def main():
    ap = argparse.ArgumentParser(description="Cluster metrics (elbow, silhouette, CH, Dunn)")
    ap.add_argument("--features", required=True, help="TSV matrix (samples x features)")
    ap.add_argument("--clusters", required=True, help="CSV with 'cluster' column")
    ap.add_argument("--k_min", type=int, default=2, help="min k for scan")
    ap.add_argument("--k_max", type=int, default=8, help="max k for scan")
    ap.add_argument("--out_prefix", required=True, help="Output prefix")
    args = ap.parse_args()

    X, labels = load_features_and_clusters(args.features, args.clusters)
    k_chosen = len(np.unique(labels[labels != -1]))  # ignore noise if any

    # Precompute distance matrix for Dunn
    print("Computing pairwise distances for Dunn index ...")
    D = pairwise_distances(X, metric="euclidean")

    # Metrics for chosen clustering
    metrics = {}

    if len(np.unique(labels)) > 1:
        try:
            sil = silhouette_score(X, labels)
        except Exception:
            sil = np.nan
    else:
        sil = np.nan

    try:
        ch = calinski_harabasz_score(X, labels)
    except Exception:
        ch = np.nan

    try:
        dunn = dunn_index_from_distance(D, labels)
    except Exception:
        dunn = np.nan

    metrics["silhouette"] = float(sil) if np.isfinite(sil) else np.nan
    metrics["calinski_harabasz"] = float(ch) if np.isfinite(ch) else np.nan
    metrics["dunn"] = float(dunn) if np.isfinite(dunn) else np.nan

    # Save metrics.tsv
    metrics_tsv = f"{args.out_prefix}_metrics.tsv"
    with open(metrics_tsv, "w") as fh:
        fh.write("metric\tvalue\n")
        for k, v in metrics.items():
            fh.write(f"{k}\t{v}\n")
        fh.write(f"n_clusters\t{k_chosen}\n")

    # Barplot of metrics (not normalized)
    plot_metrics_bar(metrics, f"{args.out_prefix}_metrics.png")

    # Scan k = k_min..k_max
    k_vals = []
    inertia_vals = []
    sil_vals = []
    ch_vals = []
    dunn_vals = []

    for k in range(args.k_min, args.k_max + 1):
        if k < 2:
            continue
        print(f"Scanning k={k} ...")
        try:
            model = KMeans(n_clusters=k, random_state=0, n_init="auto")
        except TypeError:
            model = KMeans(n_clusters=k, random_state=0)
        labels_k = model.fit_predict(X)
        inertia = model.inertia_

        if len(np.unique(labels_k)) > 1:
            try:
                sil_k = silhouette_score(X, labels_k)
            except Exception:
                sil_k = np.nan
            try:
                ch_k = calinski_harabasz_score(X, labels_k)
            except Exception:
                ch_k = np.nan
            try:
                dunn_k = dunn_index_from_distance(D, labels_k)
            except Exception:
                dunn_k = np.nan
        else:
            sil_k = np.nan
            ch_k = np.nan
            dunn_k = np.nan

        k_vals.append(k)
        inertia_vals.append(float(inertia))
        sil_vals.append(float(sil_k) if np.isfinite(sil_k) else np.nan)
        ch_vals.append(float(ch_k) if np.isfinite(ch_k) else np.nan)
        dunn_vals.append(float(dunn_k) if np.isfinite(dunn_k) else np.nan)

    # Save k_scan.tsv
    k_scan_tsv = f"{args.out_prefix}_k_scan.tsv"
    with open(k_scan_tsv, "w") as fh:
        fh.write("k\tinertia\tsilhouette\tcalinski_harabasz\tdunn\n")
        for i, k in enumerate(k_vals):
            fh.write(
                f"{k}\t{inertia_vals[i]}\t{sil_vals[i]}\t{ch_vals[i]}\t{dunn_vals[i]}\n"
            )

    # Plot elbow + metrics vs k (each separate, no normalization)
    plot_metric_vs_k(k_vals, inertia_vals, "inertia", f"{args.out_prefix}_elbow.png")
    plot_metric_vs_k(
        k_vals, sil_vals, "silhouette", f"{args.out_prefix}_silhouette_vs_k.png"
    )
    plot_metric_vs_k(
        k_vals, ch_vals, "calinski_harabasz", f"{args.out_prefix}_calinski_vs_k.png"
    )
    plot_metric_vs_k(
        k_vals, dunn_vals, "dunn", f"{args.out_prefix}_dunn_vs_k.png"
    )


if __name__ == "__main__":
    main()
