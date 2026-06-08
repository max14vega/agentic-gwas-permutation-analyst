"""Reads the permutation algorithm's raw .txt output.

Format assumption (update once the real format is confirmed): one header row
followed by whitespace-delimited rows matching `schema.COLUMN_NAMES`. The
synthetic generator in `gwas_analyst.data.synthetic` writes this same format,
so the rest of the pipeline can be developed and tested against it.
"""

from pathlib import Path

import polars as pl

from gwas_analyst.storage.schema import COLUMN_NAMES


def iter_result_files(run_dir: Path) -> list[Path]:
    """Return permutation result files in `run_dir`, sorted for determinism."""
    return sorted(Path(run_dir).glob("*.txt"))


def read_column(path: Path, column: str) -> list:
    """Read a single column from one result file without loading the rest."""
    if column not in COLUMN_NAMES:
        raise ValueError(f"Unknown column {column!r}; expected one of {COLUMN_NAMES}")
    series = pl.read_csv(
        path,
        separator="\t",
        columns=[column],
    ).to_series()
    return series.to_list()


def read_file(path: Path) -> pl.DataFrame:
    """Read a full result file as a typed DataFrame."""
    return pl.read_csv(path, separator="\t")
