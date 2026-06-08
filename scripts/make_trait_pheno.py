#!/usr/bin/env python3

import argparse
from pathlib import Path
import pandas as pd

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, help="Master phenotype TSV")
    parser.add_argument("--id-col", default="ID", help="Sample ID column")
    parser.add_argument("--trait-col", required=True, help="Trait column to extract")
    parser.add_argument("--output", required=True, help="Output 2-column TSV")
    args = parser.parse_args()

    df = pd.read_csv(args.input, sep="\t", dtype=str)

    if args.id_col not in df.columns:
        raise ValueError(f"ID column '{args.id_col}' not found. Columns: {list(df.columns)}")
    if args.trait_col not in df.columns:
        raise ValueError(f"Trait column '{args.trait_col}' not found. Columns: {list(df.columns)}")

    out = df[[args.id_col, args.trait_col]].copy()
    out = out.rename(columns={args.id_col: "sample", args.trait_col: "phenotype"})

    # drop blanks / missing values
    out["phenotype"] = pd.to_numeric(out["phenotype"], errors="coerce")
    out = out.dropna(subset=["sample", "phenotype"]).copy()

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.output, sep="\t", header=False, index=False)

if __name__ == "__main__":
    main()