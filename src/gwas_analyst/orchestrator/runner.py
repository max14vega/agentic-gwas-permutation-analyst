"""Drives the permutation pipeline per the PRD (`GWAS_Pipeline_PRD.md`) and
streams each output straight to Parquet -- the "one app, real-time
conversion" baseline.

This supersedes the original bash-loop replica: rather than globbing
pre-generated `perm_*.tsv` files (`make_permutated_pheno.py`'s output), each
permuted phenotype is generated **in memory** -- a numpy shuffle keyed by
`seed + perm_index`, written to a single rotating `NamedTemporaryFile` with
the header row pyseer's `--phenotypes` requires (`sample\\tphenotype` -- the
original scripts omit it), fed to pyseer, and deleted in a `finally` the
instant pyseer exits. At most one intermediate phenotype file exists on disk
at a time, instead of pre-generating all N up front (PRD section 7).

It also implements pyseer's LMM decomposition cache (PRD section 5.2): the
similarity-matrix decomposition is identical across every permutation of the
same population, so permutation 1 runs with `--similarity <path> --save-lmm
<prefix>` and every subsequent permutation runs with `--load-lmm <prefix>`
(omitting `--similarity` entirely) -- eliminating the most expensive
repeated computation across thousands of runs.

Because this code is itself the thing invoking `pyseer`, it knows with
certainty exactly when each `.txt` is complete -- the moment the subprocess
returns -- so no generic filesystem watcher (watchdog, polling, "has this
file stopped growing" heuristics) is needed. Convert-and-delete triggers
deterministically, in-process, right after each subprocess exits, keeping
peak disk usage bounded to roughly one extra file at a time, exactly like
`parquet_writer`'s streaming batch converters but live, alongside the run.
"""

from __future__ import annotations

import subprocess
import tempfile
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator, Literal

import numpy as np
import pandas as pd

from gwas_analyst.storage.parquet_writer import NormalizedConverter, _convert_one_inline

DEFAULT_PYSEER_ARGS = ["--lmm", "--min-af", "0.05", "--max-af", "0.95", "--cpu", "10"]


@dataclass
class RunEvent:
    """Reports the outcome of running pyseer on one phenotype permutation."""

    name: str  # e.g. "perm_0042" -- becomes the source_file
    status: Literal["converted", "pyseer_failed", "convert_failed"]
    elapsed_seconds: float
    txt_path: Path
    parquet_path: Path | None = None
    error: str | None = None


def _write_temp_phenotype(pheno_df: pd.DataFrame, perm_index: int, seed: int) -> Path:
    """Shuffle the phenotype column (seeded by `seed + perm_index`, mirroring
    the PRD's `run_single_permutation`) and write it to a single temp TSV
    with the header pyseer's `--phenotypes` requires -- a file path, not
    stdin, and the first line must be `sample\\tphenotype` (PRD section 5.1;
    `make_permutated_pheno.py` writes headerless files, which pyseer rejects).

    Caller is responsible for deleting the returned path in a `finally`.
    """
    rng = np.random.default_rng(seed + perm_index)
    perm_df = pheno_df.copy()
    perm_df["phenotype"] = rng.permutation(perm_df["phenotype"].to_numpy())

    with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as f:
        f.write("sample\tphenotype\n")
        perm_df.to_csv(f, sep="\t", header=False, index=False)
        return Path(f.name)


def _build_command(
    pyseer_cmd: list[str] | None,
    conda_env: str,
    phenotype_path: Path,
    variant_arg: tuple[str, str],
    similarity_path: Path,
    lmm_cache_prefix: Path,
    is_first: bool,
    extra_args: list[str],
) -> list[str]:
    """Build the subprocess argv for one pyseer invocation.

    `is_first` selects between the two LMM-cache modes (PRD section 5.2):
    the first permutation computes the similarity decomposition and saves it
    (`--similarity ... --save-lmm <prefix>`); every later one loads the saved
    decomposition directly (`--load-lmm <prefix>`, no `--similarity` at all).

    `pyseer_cmd`, when given, replaces only the *executable* (normally
    `conda run -n <env> pyseer`) -- e.g. a stand-in script for testing, or a
    directly-resolved `<env>/bin/pyseer` path if `conda run` doesn't suit a
    given machine. The per-file arguments are always appended identically, so
    swapping the executable doesn't change what the loop actually runs with.
    """
    base_cmd = pyseer_cmd if pyseer_cmd is not None else ["conda", "run", "-n", conda_env, "pyseer"]
    cmd = [*base_cmd, "--phenotypes", str(phenotype_path), *variant_arg]
    if is_first:
        cmd += ["--similarity", str(similarity_path), "--save-lmm", str(lmm_cache_prefix)]
    else:
        cmd += ["--load-lmm", str(lmm_cache_prefix)]
    cmd += extra_args
    return cmd


def run_permutations(
    pheno_df: pd.DataFrame,
    n_perms: int,
    seed: int,
    variant_arg: tuple[str, str],
    similarity_path: str | Path,
    txt_out_dir: str | Path,
    log_dir: str | Path,
    parquet_out_dir: str | Path,
    *,
    lmm_cache_prefix: str | Path,
    conda_env: str = "pyseer_env",
    pyseer_extra_args: list[str] | None = None,
    layout: Literal["inline", "normalized"] = "normalized",
    reference_path: str | Path | None = None,
    delete_source: bool = True,
    pyseer_cmd: list[str] | None = None,
    name_prefix: str = "perm",
) -> Iterator[RunEvent]:
    """Run pyseer once per permutation (1..n_perms, mirroring the bash loop's
    sequential `perm_NNNN` naming), generating each permuted phenotype
    in-memory and converting -- by default deleting -- each `.txt` output the
    instant its subprocess completes.

    `pheno_df` is the pre-loaded, pre-filtered 2-column (`sample`,
    `phenotype`) DataFrame for the trait being run (PRD's
    `run_permutation_batch` contract) -- replaces `make_trait_pheno.py` +
    `make_permutated_pheno.py`'s file-based handoff with an in-memory one.

    Yields one `RunEvent` per permutation as it finishes, so a CLI (or later
    a UI) can show live progress without re-deriving it. A `pyseer` failure
    or a conversion failure is reported as an event rather than raised, so
    one bad permutation doesn't abort the rest -- except a `ValueError` from
    the normalized-layout reference check, which signals a structural
    assumption violation and is allowed to propagate (mirrors
    `parquet_writer.write_normalized`'s fail-loudly behaviour).

    `variant_arg` is the pyseer flag/path pair selecting the variant input,
    e.g. `("--kmers", "/path/unitigs_reordered.pyseer.gz")`,
    `("--vcf", "/path/snpeff_annotated.vcf.gz")`, or
    `("--pres", "/path/gene_presence_absence.Rtab")`.
    """
    pyseer_extra_args = DEFAULT_PYSEER_ARGS if pyseer_extra_args is None else pyseer_extra_args
    similarity_path = Path(similarity_path)
    lmm_cache_prefix = Path(lmm_cache_prefix)
    txt_out_dir = Path(txt_out_dir)
    log_dir = Path(log_dir)
    parquet_out_dir = Path(parquet_out_dir)
    lmm_cache_prefix.parent.mkdir(parents=True, exist_ok=True)
    for d in (txt_out_dir, log_dir, parquet_out_dir):
        d.mkdir(parents=True, exist_ok=True)

    converter = NormalizedConverter(parquet_out_dir, reference_path) if layout == "normalized" else None

    for perm_index in range(1, n_perms + 1):
        name = f"{name_prefix}_{perm_index:04d}"
        txt_path = txt_out_dir / f"{name}.txt"
        log_path = log_dir / f"{name}.log"

        tmp_path = _write_temp_phenotype(pheno_df, perm_index, seed)
        try:
            cmd = _build_command(
                pyseer_cmd, conda_env, tmp_path, variant_arg, similarity_path,
                lmm_cache_prefix, perm_index == 1, pyseer_extra_args,
            )
            start = time.monotonic()
            with txt_path.open("wb") as out_f, log_path.open("wb") as err_f:
                result = subprocess.run(cmd, stdout=out_f, stderr=err_f)
            elapsed = time.monotonic() - start
        finally:
            tmp_path.unlink(missing_ok=True)  # always delete, even on error (PRD section 13)

        if result.returncode != 0:
            yield RunEvent(
                name=name, status="pyseer_failed", elapsed_seconds=elapsed, txt_path=txt_path,
                error=f"pyseer exited {result.returncode} -- see {log_path}",
            )
            continue

        try:
            if layout == "normalized":
                parquet_path = converter.convert_one(txt_path)
            else:
                parquet_path = _convert_one_inline(txt_path, parquet_out_dir)
        except ValueError:
            raise  # structural assumption violated -- fail loudly, don't keep converting
        except Exception as exc:  # noqa: BLE001 -- surfaced to the caller as an event, not raised
            yield RunEvent(
                name=name, status="convert_failed", elapsed_seconds=elapsed, txt_path=txt_path,
                error=str(exc),
            )
            continue

        if delete_source:
            txt_path.unlink()

        yield RunEvent(
            name=name, status="converted", elapsed_seconds=elapsed,
            txt_path=txt_path, parquet_path=parquet_path,
        )
