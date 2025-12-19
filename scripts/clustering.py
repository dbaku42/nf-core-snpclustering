#!/usr/bin/env python3
import argparse
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.cluster import KMeans, DBSCAN

PC_COL_RE = re.compile(r"^PC\d+$", re.IGNORECASE)


def read_table_robust(path: str) -> pd.DataFrame:
    # Prova TSV, se fallisce o viene 1 sola colonna -> whitespace
    try:
        df = pd.read_csv(path, sep="\t", dtype=str)
        if df.shape[1] == 1:
            raise ValueError("single col -> likely whitespace")
        return df
    except Exception:
        return pd.read_csv(path, sep=r"\s+", dtype=str)


def build_sample_id(df: pd.DataFrame) -> tuple[pd.Series, pd.DataFrame]:
    cols = list(df.columns)

    if "sample_id" in cols:
        sid = df["sample_id"].astype(str)
        df2 = df.drop(columns=["sample_id"])
        return sid, df2

    # flashpca2 tipico: FID IID PC1 PC2 ...
    fid_candidates = [c for c in cols if str(c).upper() == "FID"]
    iid_candidates = [c for c in cols if str(c).upper() == "IID"]
    if fid_candidates and iid_candidates:
        fid = fid_candidates[0]
        iid = iid_candidates[0]
        sid = df[fid].astype(str) + ":" + df[iid].astype(str)
        df2 = df.drop(columns=[fid, iid])
        return sid, df2

    # Caso comune: prima colonna = id (string), resto numerico
    if len(cols) >= 2:
        # prova: se la seconda colonna sembra numerica -> prima è id
        second = pd.to_numeric(df[cols[1]], errors="coerce")
        if second.notna().mean() > 0.95:
            sid = df[cols[0]].astype(str)
            df2 = df.drop(columns=[cols[0]])
            return sid, df2

    # Fallback: usa index (meno ideale, ma non rompe la pipeline)
    sid = df.index.astype(str)
    return sid, df


def select_pc_matrix(df_feats: pd.DataFrame, require_pca: bool = True) -> pd.DataFrame:
    # prova a prendere colonne PC*
    pc_cols = [c for c in df_feats.columns if PC_COL_RE.match(str(c))]
    if pc_cols:
        X = df_feats[pc_cols].apply(pd.to_numeric, errors="coerce")
        return X

    # altrimenti: usa tutte le colonne numeriche come PC (e rinomina)
    Xnum = df_feats.apply(pd.to_numeric, errors="coerce")
    Xnum = Xnum.loc[:, Xnum.notna().any(axis=0)]  # drop colonne tutte NaN

    if Xnum.shape[1] == 0:
        raise ValueError("No numeric PCA columns found in input features.")

    if require_pca:
        # Qui stiamo assumendo che il file contenga PCA scores anche senza header PC*
        Xnum.columns = [f"PC{i}" for i in range(1, Xnum.shape[1] + 1)]
        return Xnum

    # se allow non-pca: tieni numeriche e basta
    return Xnum


def load_features(path: str, require_pca: bool = True) -> tuple[np.ndarray, pd.Series, list[str]]:
    df = read_table_robust(path)
    sample_ids, df_feats = build_sample_id(df)
    Xdf = select_pc_matrix(df_feats, require_pca=require_pca)

    # PCA scores: niente NaN (fail fast)
    if Xdf.isna().any().any():
        n_missing = int(Xdf.isna().sum().sum())
        raise ValueError(
            f"Found {n_missing} missing values in PCA features. "
            "PCA scores should not contain NaN. Fix upstream."
        )

    return Xdf.values, sample_ids, list(Xdf.columns)


def main() -> None:
    ap = argparse.ArgumentParser(description="Cluster PCA-transformed features with KMeans or DBSCAN")

    ap.add_argument("--features", required=True, help="Input TSV of PCA scores (supports sample_id OR FID/IID)")
    ap.add_argument("--algorithm", choices=["kmeans", "dbscan"], default="kmeans")
    ap.add_argument("--k", type=int, default=3, help="Number of clusters for kmeans")
    ap.add_argument("--dbscan-eps", type=float, default=0.5)
    ap.add_argument("--dbscan-min-samples", type=int, default=5)
    ap.add_argument("--out-clusters", required=True)
    ap.add_argument("--out-info", required=True)

    # default: richiede PCA (PC-only). puoi disabilitare se vuoi
    ap.add_argument(
        "--allow-non-pca",
        dest="require_pca",
        action="store_false",
        help="Allow non-PCA columns (disables PCA-only enforcement). Default: PCA-only.",
    )
    ap.set_defaults(require_pca=True)

    args = ap.parse_args()

    X, sample_ids, feat_names = load_features(args.features, require_pca=args.require_pca)

    info: dict = {
        "algorithm": args.algorithm,
        "n_samples": int(X.shape[0]),
        "n_features": int(X.shape[1]),
        "features_type": "pca_scores" if args.require_pca else "numeric_features",
        "feature_names_head": feat_names[:10],
    }

    if args.algorithm == "kmeans":
        if args.k < 1:
            raise ValueError("k must be >= 1")
        try:
            model = KMeans(n_clusters=args.k, n_init=100, init="random", random_state=42)
        except TypeError:
            model = KMeans(n_clusters=args.k, n_init=10, random_state=42)

        labels = model.fit_predict(X)
        info.update({"k": int(args.k), "inertia": float(model.inertia_)})
    else:
        model = DBSCAN(eps=args.dbscan_eps, min_samples=args.dbscan_min_samples)
        labels = model.fit_predict(X)
        uniq = set(labels.tolist())
        n_clusters = len(uniq) - (1 if -1 in uniq else 0)
        info.update({
            "eps": float(args.dbscan_eps),
            "min_samples": int(args.dbscan_min_samples),
            "n_clusters_found": int(n_clusters),
            "n_noise": int(np.sum(labels == -1)),
        })

    out_df = pd.DataFrame({"sample_id": sample_ids.astype(str), "cluster": labels.astype(int)})
    out_df.to_csv(args.out_clusters, sep="\t", index=False)
    Path(args.out_info).write_text(json.dumps(info, indent=2))


if __name__ == "__main__":
    main()
