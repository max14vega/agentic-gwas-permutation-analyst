#!/usr/bin/env python3
"""Run the permutation pipeline and talk to the analyst -- in one app.

Generates each permuted phenotype in memory (replacing `make_trait_pheno.py`
+ `make_permutated_pheno.py`'s file-based handoff -- PRD section 7), drives
`pyseer` once per permutation with LMM-cache reuse after the first run (PRD
section 5.2), converts each output to Parquet -- deleting the bloated `.txt`
the instant it's produced -- and, once the run finishes, drops straight into
the existing `Analyst` REPL pointed at the freshly-converted dataset. This is
"open one application, run your tests, watch them convert in real time, then
ask questions" end to end.

Usage:
    python scripts/run_and_analyze.py \\
        --master-pheno phenotypes.tsv --trait-col LD50 \\
        --kmers variant_files/unitigs_reordered.pyseer.gz \\
        --similarity recombination_masked_similarity.tsv \\
        --population LD50 --base-dir /data/gwas_run \\
        --n-perms 10000 --seed 84 \\
        [--layout inline|normalized] [--conda-env pyseer_env] [--cpu 10]
        [question]                                    # ask once and exit

Requires ANTHROPIC_API_KEY in the environment for the analyst phase.
"""

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from gwas_analyst.agent.analyst import Analyst  # noqa: E402
from gwas_analyst.orchestrator.runner import run_permutations  # noqa: E402
from gwas_analyst.storage import db  # noqa: E402


def _variant_arg(args: argparse.Namespace) -> tuple[str, str]:
    for flag, value in (("--kmers", args.kmers), ("--vcf", args.vcf), ("--pres", args.pres)):
        if value:
            return (flag, value)
    raise SystemExit("one of --kmers / --vcf / --pres is required")


def _load_trait_pheno(master_pheno: Path, id_col: str, trait_col: str) -> pd.DataFrame:
    """Extract a 2-column (sample, phenotype) frame for one trait, dropping
    blanks/NA -- replaces `make_trait_pheno.py`'s file-based step (its output
    is now generated in memory and re-shuffled per permutation, never written
    to disk as an intermediate)."""
    df = pd.read_csv(master_pheno, sep="\t", dtype=str)
    if id_col not in df.columns:
        raise SystemExit(f"ID column {id_col!r} not found. Columns: {list(df.columns)}")
    if trait_col not in df.columns:
        raise SystemExit(f"trait column {trait_col!r} not found. Columns: {list(df.columns)}")

    out = df[[id_col, trait_col]].rename(columns={id_col: "sample", trait_col: "phenotype"})
    out["phenotype"] = pd.to_numeric(out["phenotype"], errors="coerce")
    return out.dropna(subset=["sample", "phenotype"]).reset_index(drop=True)


def _print_event(i: int, total: int, event) -> None:
    progress = f"[{i:>5}/{total}]"
    if event.status == "converted":
        size = event.parquet_path.stat().st_size
        print(f"{progress} {event.name} -> {event.parquet_path.name} "
              f"({size / 1e6:.1f} MB, {event.elapsed_seconds:.0f}s), source deleted")
    else:
        print(f"{progress} {event.name}: {event.status} -- {event.error}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--master-pheno", required=True, type=Path, help="master phenotype TSV (one column per trait)")
    parser.add_argument("--id-col", default="ID", help="sample ID column in the master phenotype TSV [default: ID]")
    parser.add_argument("--trait-col", required=True, help="trait column to run permutations for")
    variant = parser.add_mutually_exclusive_group(required=True)
    variant.add_argument("--kmers", help="path to a pyseer k-mers/unitigs file")
    variant.add_argument("--vcf", help="path to a pyseer VCF variant file")
    variant.add_argument("--pres", help="path to a pyseer presence/absence Rtab file")
    parser.add_argument("--similarity", required=True, help="path to the (recombination-masked) similarity table -- used for permutation 1 only, then cached via --save-lmm/--load-lmm")
    parser.add_argument("--population", required=True, help="population label, e.g. 'ALL_CA' -- namespaces output dirs and the LMM cache")
    parser.add_argument("--base-dir", required=True, type=Path, help="base directory results/logs/cache are written under")
    parser.add_argument("--n-perms", type=int, default=10000, help="number of permutations to run [default: 10000]")
    parser.add_argument("--seed", type=int, default=84, help="random seed for phenotype shuffling [default: 84]")
    parser.add_argument("--layout", choices=["inline", "normalized"], default="normalized", help="Parquet layout -- run investigate_variant_redundancy.py first if unsure [default: normalized]")
    parser.add_argument("--conda-env", default="pyseer_env", help="conda environment pyseer runs in [default: pyseer_env]")
    parser.add_argument("--cpu", type=int, default=10, help="--cpu passed to pyseer [default: 10]")
    parser.add_argument("--min-af", default="0.05", help="--min-af passed to pyseer [default: 0.05]")
    parser.add_argument("--max-af", default="0.95", help="--max-af passed to pyseer [default: 0.95]")
    parser.add_argument("--keep-txt", action="store_true", help="keep the raw .txt output instead of deleting it after conversion")
    parser.add_argument("--pyseer-cmd", help="override the pyseer executable (space-separated), bypassing 'conda run -n <env> pyseer' -- e.g. for a stand-in/test binary or a directly-resolved env path")
    parser.add_argument("question", nargs="?", help="ask the analyst this one question after the run finishes, then exit")
    args = parser.parse_args()

    variant_arg = _variant_arg(args)
    extra_args = ["--lmm", "--min-af", args.min_af, "--max-af", args.max_af, "--cpu", str(args.cpu)]
    pyseer_cmd = args.pyseer_cmd.split() if args.pyseer_cmd else None

    pheno_df = _load_trait_pheno(args.master_pheno, args.id_col, args.trait_col)
    variant_stem = Path(variant_arg[1]).name.split(".")[0]  # strip all suffixes, e.g. "unitigs.pyseer.gz" -> "unitigs"
    job_name = f"{args.trait_col}_{variant_stem}"
    txt_out = args.base_dir / "permutations" / args.population / job_name
    log_dir = args.base_dir / "perm_logs" / args.population / job_name
    parquet_out = args.base_dir / "parquet" / args.population / job_name / "permutations"
    lmm_cache_prefix = args.base_dir / "lmm_cache" / args.population / "lmm"

    print(f"Running {args.n_perms} permutations of {args.trait_col!r} "
          f"({len(pheno_df)} samples, seed={args.seed}) -- "
          f"converting and deleting each .txt as it lands.\n")

    try:
        for i, event in enumerate(run_permutations(
            pheno_df, args.n_perms, args.seed, variant_arg, args.similarity,
            txt_out, log_dir, parquet_out,
            lmm_cache_prefix=lmm_cache_prefix, conda_env=args.conda_env,
            pyseer_extra_args=extra_args, layout=args.layout,
            delete_source=not args.keep_txt, pyseer_cmd=pyseer_cmd,
        ), start=1):
            _print_event(i, args.n_perms, event)
    except KeyboardInterrupt:
        print("\nInterrupted -- analyzing what's converted so far.")

    print()
    reference = parquet_out.parent / "reference" / "variants.parquet"
    con = db.connect(parquet_out, reference_path=reference if args.layout == "normalized" else None)
    analyst = Analyst(con)

    if args.question:
        print(analyst.ask(args.question))
        return

    print("Agentic GWAS analyst -- type a question, or Ctrl-D to quit.")
    try:
        while True:
            question = input("\n> ").strip()
            if not question:
                continue
            print(analyst.ask(question))
    except (EOFError, KeyboardInterrupt):
        print()


if __name__ == "__main__":
    main()
