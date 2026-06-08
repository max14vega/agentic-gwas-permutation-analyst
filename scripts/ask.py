#!/usr/bin/env python3
"""CLI for the agentic analyst -- the baseline interface before any UI exists.

Usage:
    python scripts/ask.py <permutations_glob_or_dir> [--reference path/to/variants.parquet] "question"
    python scripts/ask.py <permutations_glob_or_dir> [--reference ...]            # interactive REPL

Requires ANTHROPIC_API_KEY in the environment.
"""

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from gwas_analyst.agent.analyst import Analyst  # noqa: E402
from gwas_analyst.storage import db  # noqa: E402


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("permutations", help="directory or glob of permutation Parquet files")
    parser.add_argument("--reference", help="path to variants.parquet (normalized layouts)")
    parser.add_argument("question", nargs="?", help="ask a single question and exit")
    args = parser.parse_args()

    con = db.connect(args.permutations, reference_path=args.reference)
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
