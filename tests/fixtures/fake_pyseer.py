#!/usr/bin/env python3
"""Stand-in for the real `pyseer` binary, used to exercise the orchestrator
end to end without a conda env / variant dataset / similarity matrix.

Accepts the same CLI shape the pipeline invokes pyseer with (`--phenotypes`,
`--kmers`/`--vcf`/`--pres`, `--similarity`/`--save-lmm`/`--load-lmm`, `--lmm`,
`--min-af`, `--max-af`, `--cpu`, ...). Writes an 8-column header + N rows to
**stdout** in pyseer's real shape (tab-delimited, scientific-notation floats,
occasional empty-string nulls), exactly as the orchestrator's
`> "${name}.txt"` redirection expects.

Two extra checks exercise behaviour the orchestrator depends on but that a
synthetic-data-shaped stub wouldn't otherwise touch:

* **Phenotype header**: pyseer requires `--phenotypes` to point at a file
  whose first line is a header (PRD section 5.1: `sample\\tphenotype`); the
  stub fails loudly if that header is missing, the same way real pyseer would
  reject a headerless file from `make_permutated_pheno.py`.
* **LMM cache round-trip**: `--save-lmm <prefix>` writes a `<prefix>.cache`
  marker file; `--load-lmm <prefix>` requires that marker to already exist
  (and errors if it doesn't, like real pyseer would on a missing/stale
  decomposition) -- exercising the orchestrator's "permutation 1 saves, 2..N
  load" sequencing (PRD section 5.2).

The variant set is **fixed** (seeded from `--phenotypes`'s row count alone,
not its content) so every invocation produces the same variants in the same
order -- mirroring standard permutation testing (only the phenotype is
shuffled) and exercising the orchestrator's `normalized` layout path.
"""

import argparse
import random
import sys
from decimal import Decimal
from pathlib import Path

_BASES = "ACGT"
_NOTES_VOCAB = ["af-filter", "bad-chisq", "high-bse", "firth-fail"]
_HEADER = ["variant", "af", "filter-pvalue", "lrt-pvalue", "beta", "beta-std-err", "variant_h2", "notes"]


def _sci(x: float, rng: random.Random) -> str:
    if rng.random() < 0.01:
        return ""
    return "%.2E" % Decimal(x)


def _check_phenotype_header(pheno_path: Path) -> int:
    with pheno_path.open() as f:
        header = f.readline().rstrip("\n")
        if header != "sample\tphenotype":
            sys.stderr.write(
                f"fake_pyseer: ERROR -- {pheno_path} is missing the required "
                f"'sample\\tphenotype' header (got {header!r}); pyseer requires "
                "a header row in --phenotypes files\n"
            )
            sys.exit(1)
        return sum(1 for _ in f)


def _check_lmm_cache(args: argparse.Namespace) -> None:
    if args.save_lmm:
        marker = Path(f"{args.save_lmm}.cache")
        marker.parent.mkdir(parents=True, exist_ok=True)
        marker.write_text("fake-lmm-decomposition\n")
    if args.load_lmm:
        marker = Path(f"{args.load_lmm}.cache")
        if not marker.exists():
            sys.stderr.write(
                f"fake_pyseer: ERROR -- --load-lmm {args.load_lmm} has no cached "
                f"decomposition at {marker} (expected --save-lmm to have run first)\n"
            )
            sys.exit(1)
        if args.similarity:
            sys.stderr.write(
                "fake_pyseer: ERROR -- --load-lmm and --similarity are mutually "
                "exclusive (the cached decomposition replaces the similarity matrix)\n"
            )
            sys.exit(1)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--phenotypes", required=True)
    parser.add_argument("--similarity")
    parser.add_argument("--save-lmm")
    parser.add_argument("--load-lmm")
    parser.add_argument("--num-rows", type=int, default=200)
    args, _ignored = parser.parse_known_args()

    pheno_path = Path(args.phenotypes)
    n_samples = _check_phenotype_header(pheno_path)
    _check_lmm_cache(args)

    shared_rng = random.Random(1234)  # fixed seed -> identical variant set every invocation
    variants = ["".join(shared_rng.choice(_BASES) for _ in range(30)) for _ in range(args.num_rows)]

    row_rng = random.Random(hash(pheno_path.read_bytes()) & 0xFFFFFFFF)
    out = sys.stdout
    out.write("\t".join(_HEADER) + "\n")
    for variant in variants:
        notes = "" if row_rng.random() > 0.05 else row_rng.choice(_NOTES_VOCAB)
        row = [
            variant,
            _sci(row_rng.uniform(0.01, 0.5), row_rng),
            _sci(row_rng.uniform(0, 1), row_rng),
            _sci(row_rng.uniform(0, 1), row_rng),
            _sci(row_rng.uniform(-1, 1), row_rng),
            _sci(row_rng.uniform(0, 0.5), row_rng),
            _sci(row_rng.uniform(0, 0.05), row_rng),
            notes,
        ]
        out.write("\t".join(row) + "\n")
    sys.stderr.write(f"fake_pyseer: {n_samples} samples, {args.num_rows} variants -- done\n")


if __name__ == "__main__":
    main()
