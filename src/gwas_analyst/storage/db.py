"""DuckDB query layer over the Parquet dataset -- the "database baseline".

DuckDB queries Parquet files directly (out-of-core, no server process), which
fits a single-user setup: point it at a directory of files and run SQL.

Two dataset shapes are supported, mirroring `parquet_writer`:
  * normalized -- `permutations` (fact rows + variant_id) joined with
    `variants` (variant_id -> sequence) via the `results` view.
  * inline -- `permutations` already contains the full typed rows; `results`
    is just an alias for it.
"""

import re
from pathlib import Path
from typing import Any

import duckdb

_SELECT_ONLY = re.compile(r"^\s*(with\b|select\b)", re.IGNORECASE)


def connect(
    permutations_glob: str | Path,
    reference_path: str | Path | None = None,
) -> duckdb.DuckDBPyConnection:
    """Open an in-memory DuckDB connection with views over the Parquet dataset.

    `permutations_glob` is a glob pattern (or directory) matching the
    per-permutation Parquet files, e.g. "data/parquet_normalized/permutations/*.parquet".
    `reference_path`, if given, points at the variant reference table -- pass
    this for normalized layouts to get a `results` view with `variant` joined back in.
    """
    con = duckdb.connect(":memory:")
    glob = str(permutations_glob)
    if Path(glob).is_dir():
        glob = str(Path(glob) / "*.parquet")

    con.execute(
        "CREATE VIEW permutations AS "
        "SELECT *, regexp_extract(filename, '([^/\\\\]+)\\.parquet$', 1) AS source_file "
        f"FROM read_parquet('{glob}', filename = true)"
    )

    if reference_path is not None:
        con.execute(f"CREATE VIEW variants AS SELECT * FROM read_parquet('{reference_path}')")
        con.execute(
            "CREATE VIEW results AS "
            "SELECT p.*, v.variant "
            "FROM permutations p JOIN variants v USING (variant_id)"
        )
    else:
        con.execute("CREATE VIEW results AS SELECT * FROM permutations")

    return con


def run_readonly_query(con: duckdb.DuckDBPyConnection, sql: str) -> list[dict[str, Any]]:
    """Execute a read-only SQL query and return rows as dicts.

    Only SELECT/WITH statements are permitted -- this is the guard the agent's
    SQL tool relies on so a model-generated query can't mutate the dataset.
    """
    if not _SELECT_ONLY.match(sql):
        raise ValueError("Only SELECT/WITH statements are allowed")
    if ";" in sql.strip().rstrip(";"):
        raise ValueError("Multiple statements are not allowed")

    result = con.execute(sql)
    columns = [d[0] for d in result.description]
    return [dict(zip(columns, row)) for row in result.fetchall()]


# -- Canonical analytical queries -------------------------------------------------

def pooled_null_distribution(
    con: duckdb.DuckDBPyConnection, column: str = "lct_pvalue"
) -> dict[str, Any]:
    """Summary stats of `column` pooled across every permutation -- the basis
    for an empirical null distribution / permutation-derived significance threshold."""
    sql = f"""
        SELECT
            count(*)                          AS n,
            min({column})                     AS min,
            quantile_cont({column}, 0.01)     AS p01,
            quantile_cont({column}, 0.50)     AS median,
            quantile_cont({column}, 0.99)     AS p99,
            max({column})                     AS max,
            avg({column})                     AS mean,
            stddev({column})                  AS stddev
        FROM permutations
    """
    rows = run_readonly_query(con, sql)
    return rows[0]


def top_hits(
    con: duckdb.DuckDBPyConnection,
    column: str = "lct_pvalue",
    limit: int = 20,
    ascending: bool = True,
) -> list[dict[str, Any]]:
    """Top `limit` rows by `column` (smallest p-values by default)."""
    order = "ASC" if ascending else "DESC"
    sql = f"SELECT * FROM results ORDER BY {column} {order} LIMIT {int(limit)}"
    return run_readonly_query(con, sql)


def per_run_summary(con: duckdb.DuckDBPyConnection, column: str = "lct_pvalue") -> list[dict[str, Any]]:
    """Aggregate `column` stats grouped by source permutation file."""
    sql = f"""
        SELECT
            source_file,
            count(*)            AS n,
            min({column})       AS min,
            avg({column})       AS mean,
            max({column})       AS max
        FROM permutations
        GROUP BY source_file
        ORDER BY source_file
    """
    return run_readonly_query(con, sql)
