#!/usr/bin/env python3
import argparse
import pandas as pd

def read_robust(path: str) -> pd.DataFrame:
    try:
        df = pd.read_csv(path, sep="\t")
        if df.shape[1] == 1:
            raise ValueError("single col -> maybe whitespace")
        return df
    except Exception:
        return pd.read_csv(path, sep=r"\s+", engine="python")

def is_num(s: pd.Series) -> bool:
    return pd.to_numeric(s, errors="coerce").notna().mean() > 0.95

def main():
    ap = argparse.ArgumentParser(description="Convert pca_scores to features TSV (sample_id + numeric PCs).")
    ap.add_argument("--in", dest="inp", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    df = read_robust(args.inp)

    # Case A: already has sample_id
    if "sample_id" in df.columns:
        sid = df["sample_id"].astype(str)
        X = df.drop(columns=["sample_id"]).copy()

    # Case B: flashpca-like: FID IID PC1...
    elif df.shape[1] >= 3 and (not is_num(df.iloc[:, 0])) and (not is_num(df.iloc[:, 1])):
        sid = df.iloc[:, 0].astype(str) + ":" + df.iloc[:, 1].astype(str)
        X = df.iloc[:, 2:].copy()

    # Case C: id + numeric columns
    elif df.shape[1] >= 2 and (not is_num(df.iloc[:, 0])) and is_num(df.iloc[:, 1]):
        sid = df.iloc[:, 0].astype(str)
        X = df.iloc[:, 1:].copy()

    else:
        raise SystemExit("ERROR: could not infer sample_id from pca_scores")

    X = X.apply(pd.to_numeric, errors="coerce")
    X = X.loc[:, X.notna().any(axis=0)].fillna(0.0)

    # ensure PC-like column names
    X.columns = [
        c if str(c).lower().startswith("pc") else f"PC{i}"
        for i, c in enumerate(X.columns, start=1)
    ]

    out_df = pd.concat([sid.rename("sample_id"), X], axis=1)
    out_df.to_csv(args.out, sep="\t", index=False)

if __name__ == "__main__":
    main()
