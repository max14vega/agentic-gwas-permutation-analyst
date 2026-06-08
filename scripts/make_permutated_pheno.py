#!/usr/bin/env python3

import argparse
from pathlib import Path
import numpy as np
import pandas as pd

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--phenotype", required=True, help="2-column phenotype TSV")
    parser.add_argument("--outdir", required=True, help="Output dir")
    parser.add_argument("--n", type=int, required=True, help="Number of permutations")
    parser.add_argument("--seed", type=int, default=1)
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    pheno = pd.read_csv(args.phenotype, sep="\t", header=None, names=["sample", "phenotype"])
    rng = np.random.default_rng(args.seed)

    for i in range(1, args.n + 1):
        perm = pheno.copy()
        perm["phenotype"] = rng.permutation(perm["phenotype"].values)
        perm.to_csv(outdir / f"perm_{i:03d}.tsv", sep="\t", header=False, index=False)

if __name__ == "__main__":
    main()