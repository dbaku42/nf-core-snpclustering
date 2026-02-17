#!/usr/bin/env python3
import argparse
import json
import re
from pathlib import Path

def main():
    ap = argparse.ArgumentParser(description="Convert flashpca2 .pcs/.outpc to TSV with header.")

    # nomi canonici
    ap.add_argument("--in", dest="inp", required=False,
                    help="Input flashpca2 PCs file (e.g. .pcs or .outpc)")
    ap.add_argument("--out", dest="out", required=False,
                    help="Output TSV path")

    # alias compatibili con la pipeline Nextflow
    ap.add_argument("--outpc", dest="inp", required=False,
                    help="Alias of --in (flashpca outpc.raw)")
    ap.add_argument("--out-pca", dest="out", required=False,
                    help="Alias of --out (output TSV)")
    ap.add_argument("--n-pcs", dest="n_pcs", type=int, default=None,
                    help="Optional (compat)")
    ap.add_argument("--out-info", dest="out_info", default=None,
                    help="Optional JSON output (compat)")
    ap.add_argument("--id-mode", dest="id_mode", default=None,
                    help="Ignored/optional (compat)")

    ap.add_argument("--id_col", default="IID", help="Name for sample ID column")
    args = ap.parse_args()

    # valida input/output (accetta sia --in/--out sia --outpc/--out-pca)
    if not args.inp or not args.out:
        ap.error("Missing required args: provide --in/--out or --outpc/--out-pca")

    inp = Path(args.inp)
    out = Path(args.out)

    if not inp.exists():
        raise SystemExit(f"Input file not found: {inp}")

    lines = inp.read_text().splitlines()
    # drop empty/comment lines
    rows = [ln for ln in lines if ln.strip() and not ln.lstrip().startswith("#")]
    if not rows:
        raise SystemExit(f"No data in {inp}")

    # flashpca2 output is typically: <sample_id> <PC1> <PC2> ...
    first = re.split(r"\s+", rows[0].strip())
    if len(first) < 2:
        raise SystemExit(f"Unexpected format in {inp}: first row has <2 columns")

    n_pc_found = len(first) - 1
    # se la pipeline passa --n-pcs, rispettiamolo (senza rompere se differisce)
    n_pc = min(args.n_pcs, n_pc_found) if args.n_pcs else n_pc_found

    header = [args.id_col] + [f"PC{i}" for i in range(1, n_pc + 1)]

    out.parent.mkdir(parents=True, exist_ok=True)
    n_written = 0
    with out.open("w") as f:
        f.write("\t".join(header) + "\n")
        for ln in rows:
            parts = re.split(r"\s+", ln.strip())
            if len(parts) < n_pc + 1:
                continue
            # scrivi solo ID + primi n_pc
            f.write("\t".join([parts[0]] + parts[1:1+n_pc]) + "\n")
            n_written += 1

    # crea pca_info.json se richiesto dalla pipeline
    if args.out_info:
        info_path = Path(args.out_info)
        info = {
            "input": str(inp.name),
            "output": str(out.name),
            "n_pcs": int(n_pc),
            "n_samples_written": int(n_written),
            "id_col": str(args.id_col),
            "id_mode": args.id_mode,
        }
        info_path.write_text(json.dumps(info, indent=2) + "\n")

if __name__ == "__main__":
    main()

