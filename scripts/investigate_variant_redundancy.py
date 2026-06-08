#!/usr/bin/env python3
"""Check whether the `variant` column is reused across permutation files.

This answers the single biggest storage-architecture question: standard GWAS
permutation testing re-tests the same fixed variant set across all runs (only
the phenotype/label is shuffled). If that holds here, the variant sequences
can be normalized into a single reference table instead of being duplicated
in every one of the ~10,000 output files -- roughly halving (or more) the raw
data volume before any compression.

Usage:
    python scripts/investigate_variant_redundancy.py <run_dir> [--sample N] [--column variant]

`run_dir` should contain the permutation .txt files (or any directory that
`storage.txt_reader` can iterate -- see assumptions noted there about format).
"""

import argparse
import sys
from pathlib import Path

import polars as pl

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from gwas_analyst.storage.txt_reader import iter_result_files, read_column  # noqa: E402


def investigate(run_dir: Path, column: str = "variant", sample: int = 5) -> None:
    files = iter_result_files(run_dir)
    if not files:
        raise SystemExit(f"No result files found in {run_dir}")

    sample_files = files[:sample]
    print(f"Found {len(files)} files; inspecting {len(sample_files)} of them.\n")

    columns = [read_column(f, column) for f in sample_files]
    lengths = [len(c) for c in columns]
    print(f"Row counts per sampled file: {lengths}")

    base = columns[0]
    same_as_first = [c == base for c in columns]
    print(f"Identical to first file (same values, same order): {same_as_first}")

    base_set = set(base)
    same_set = [set(c) == base_set for c in columns]
    print(f"Same SET of values as first file (order-independent):  {same_set}")

    union = set().union(*columns)
    print(f"\nUnique '{column}' values across sampled files: {len(union)}")
    print(f"Sum of per-file unique counts:                  {sum(len(set(c)) for c in columns)}")

    if all(same_as_first):
        verdict = (
            "IDENTICAL across files (including order) -> normalize into a reference\n"
            "table keyed by a row index; fact tables can drop the column entirely."
        )
    elif all(same_set):
        verdict = (
            "Same SET of values, different order -> normalize into a reference table\n"
            "keyed by variant sequence (variant -> variant_id); fact tables store variant_id."
        )
    elif len(union) < sum(len(set(c)) for c in columns) * 0.8:
        verdict = (
            "Substantial overlap but not identical sets -> still worth a global\n"
            "reference table (dictionary-encode the union); fact tables store variant_id."
        )
    else:
        verdict = (
            "Little to no overlap across files -> normalization won't help much;\n"
            "rely on Parquet dictionary/RLE encoding (and consider 2-bit packing\n"
            "for long pure-ACGT sequences) instead."
        )

    print(f"\nVerdict: {verdict}")
    print(
        "\n(Sampled only the first "
        f"{len(sample_files)} of {len(files)} files -- re-run with --sample "
        "to widen the check before committing to a layout.)"
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run_dir", type=Path)
    parser.add_argument("--column", default="variant")
    parser.add_argument("--sample", type=int, default=5, help="number of files to compare")
    args = parser.parse_args()

    investigate(args.run_dir, column=args.column, sample=args.sample)
