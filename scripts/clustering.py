#!/usr/bin/env python3
import argparse
import os
import numpy as np
import pandas as pd
from sklearn.cluster import KMeans, DBSCAN
from sklearn.metrics import silhouette_score
from sklearn.preprocessing import StandardScaler


def load_matrix(path):
    if not os.path.exists(path):
        raise SystemExit(f"ERROR: input file not found: {path}")
    df = pd.read_csv(path, sep="\t")
    if df.shape[1] < 2:
        raise SystemExit("ERROR: matrix must have at least one feature column")

    samples = df.iloc[:, 0].astype(str).tolist()
    X_raw = df.iloc[:, 1:]

    # Try numeric conversion first
    X_num = X_raw.apply(pd.to_numeric, errors="coerce")

    # If there are NaNs, interpret as genotype codes (0|1, 1/1, ./.)
    if X_num.isna().any().any():
        def gt_to_dosage(gt):
            if pd.isna(gt):
                return np.nan
            gt = str(gt)
            if gt in ("./.", "."):
                return np.nan
            gt = gt.replace("|", "/")
            alleles = gt.split("/")
            try:
                return sum(1 for a in alleles if a == "1")
            except Exception:
                return np.nan

        X_num = X_raw.applymap(gt_to_dosage)

    X_num = X_num.fillna(-1.0)
    X = X_num.to_numpy(dtype=float)
    return samples, X


def main():
    ap = argparse.ArgumentParser(description="Clustering on genotype matrix")
    ap.add_argument("--input", required=True, help="TSV matrix with first column 'sample'")
    ap.add_argument("--n_clusters", type=int, default=4, help="Number of clusters for KMeans")
    ap.add_argument("--method", choices=["kmeans", "dbscan"], default="kmeans", help="Clustering method")
    ap.add_argument("--db_eps", type=float, default=0.5, help="DBSCAN eps")
    ap.add_argument("--db_min_samples", type=int, default=5, help="DBSCAN min_samples")
    ap.add_argument("--standardize", action="store_true", help="Standardize features before clustering")
    ap.add_argument("--out_prefix", required=True, help="Output prefix")
    args = ap.parse_args()

    samples, X = load_matrix(args.input)

    if args.standardize:
        X = StandardScaler().fit_transform(X)

    if args.method == "kmeans":
        if args.n_clusters < 1:
            raise SystemExit("ERROR: n_clusters must be >= 1 for KMeans")
        try:
            model = KMeans(n_clusters=args.n_clusters, random_state=0, n_init="auto")
        except TypeError:
            model = KMeans(n_clusters=args.n_clusters, random_state=0)
        labels = model.fit_predict(X)
    else:
        model = DBSCAN(eps=args.db_eps, min_samples=args.db_min_samples)
        labels = model.fit_predict(X)

    uniq = set(labels)
    if len(uniq) > 1:
        try:
            sil = silhouette_score(X, labels)
            print(f"Silhouette score: {sil:.4f}")
        except Exception as e:
            print(f"Could not compute silhouette: {e}")

    df_out = pd.DataFrame({"sample": samples, "cluster": labels})
    out_csv = f"{args.out_prefix}_clusters.csv"
    df_out.to_csv(out_csv, index=False)
    print(f"Saved clusters to {out_csv}")


if __name__ == "__main__":
    main()
