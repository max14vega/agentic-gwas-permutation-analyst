"""Convert raw permutation .txt output into compact, typed Parquet.

Two layouts are supported -- pick based on what
`scripts/investigate_variant_redundancy.py` reports for the real data:

* `write_inline`: each file becomes one Parquet file with the full 12-column
  typed schema (zstd-compressed, dictionary-encoded strings). Use this when
  variant sets differ across permutations.

* `write_normalized`: splits the `variant` column out into a single
  `reference/variants.parquet` (written once) and writes fact tables that
  carry only `variant_id` plus the per-permutation statistics. Use this when
  the same variant set (in the same order) is reused across every file --
  the expected case for standard permutation testing -- since it removes the
  largest column's ~10,000x duplication entirely.

Both functions stream one source file at a time (read -> cast -> write -> drop)
so peak memory and, if `delete_source=True`, peak disk usage stay bounded --
critical given the source data alone can be 300-400 GB.
"""

from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq

from gwas_analyst.storage.schema import (
    ARROW_SCHEMA,
    NORMALIZED_FACT_SCHEMA,
    VARIANT_REFERENCE_SCHEMA,
)
from gwas_analyst.storage.txt_reader import iter_result_files

COMPRESSION = "zstd"
COMPRESSION_LEVEL = 9


def _read_typed_table(path: Path, schema: pa.Schema) -> pa.Table:
    """Read a tab-delimited file and cast columns to the target Arrow schema.

    pyseer's real header uses hyphens for some fields (`filter-pvalue`,
    `lrt-pvalue`, `beta-std-err` -- see __main__.py:488-513) while the schema
    uses underscores throughout (so SQL/agent code can reference columns
    without quoting). Read with automatic type inference first -- the raw
    column names won't match an explicit `column_types` map keyed by schema
    names -- then normalize names and cast to the typed schema.
    """
    import pyarrow.csv as pcsv

    table = pcsv.read_csv(path, parse_options=pcsv.ParseOptions(delimiter="\t"))
    table = table.rename_columns([name.replace("-", "_") for name in table.column_names])
    return table.select(schema.names).cast(schema)


def _convert_one_inline(src: Path, out_dir: Path) -> Path:
    """Convert a single result file to a standalone typed, compressed Parquet
    file. Factored out so both the batch `write_inline` and a streaming caller
    (e.g. the orchestrator, converting each file the moment pyseer produces it)
    share one read -> cast -> write code path."""
    table = _read_typed_table(src, ARROW_SCHEMA)
    dst = out_dir / (src.stem + ".parquet")
    pq.write_table(table, dst, compression=COMPRESSION, compression_level=COMPRESSION_LEVEL)
    return dst


def write_inline(run_dir: Path, out_dir: Path, delete_source: bool = False) -> list[Path]:
    """Convert each result file to a standalone typed, compressed Parquet file."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    written = []
    for src in iter_result_files(run_dir):
        written.append(_convert_one_inline(src, out_dir))
        if delete_source:
            src.unlink()
    return written


class NormalizedConverter:
    """Stateful per-file converter for the normalized layout.

    Establishes the variant reference table from the *first* file it converts
    (variant_id = row index) and validates every subsequent file's `variant`
    column against it, raising ValueError on mismatch -- a violated "identical
    order" assumption fails loudly rather than silently corrupting the dataset.

    Factored out as a class (rather than a loop inside `write_normalized`) so
    the batch writer and a streaming caller (the orchestrator, converting each
    file right after the `pyseer` subprocess that produced it exits) can share
    the same incremental reference-building logic one file at a time.
    """

    def __init__(self, out_dir: Path, reference_path: Path | None = None):
        self.out_dir = Path(out_dir)
        self.out_dir.mkdir(parents=True, exist_ok=True)
        self.reference_path = reference_path or (self.out_dir.parent / "reference" / "variants.parquet")
        self.reference_path.parent.mkdir(parents=True, exist_ok=True)
        self._reference_variants: pa.Array | None = None

    def convert_one(self, src: Path) -> Path:
        table = _read_typed_table(src, ARROW_SCHEMA)
        variants = table.column("variant")

        if self._reference_variants is None:
            self._reference_variants = variants
            ref_table = pa.table(
                {
                    "variant_id": pa.array(range(len(variants)), type=pa.int32()),
                    "variant": variants,
                },
                schema=VARIANT_REFERENCE_SCHEMA,
            )
            pq.write_table(ref_table, self.reference_path, compression=COMPRESSION,
                           compression_level=COMPRESSION_LEVEL)
        elif not variants.equals(self._reference_variants):
            raise ValueError(
                f"{src.name}: variant column differs from the established reference -- "
                "the 'identical order' assumption does not hold for this dataset. "
                "Re-run investigate_variant_redundancy.py and use write_inline "
                "(or a set-based normalization) instead."
            )

        fact_table = table.drop_columns(["variant"])
        fact_table = fact_table.add_column(
            0, "variant_id", pa.array(range(len(variants)), type=pa.int32())
        )
        fact_table = fact_table.select(NORMALIZED_FACT_SCHEMA.names)

        dst = self.out_dir / (src.stem + ".parquet")
        pq.write_table(fact_table, dst, compression=COMPRESSION,
                       compression_level=COMPRESSION_LEVEL)
        return dst


def write_normalized(
    run_dir: Path,
    out_dir: Path,
    reference_path: Path | None = None,
    delete_source: bool = False,
) -> list[Path]:
    """Convert results to a normalized layout: one shared reference table of
    unique variants plus per-file fact tables keyed by `variant_id`.

    Assumes every file has the *same variants in the same row order* -- the
    expected layout for standard permutation testing (only labels are
    shuffled). See `NormalizedConverter` for the per-file conversion logic.
    """
    files = iter_result_files(run_dir)
    if not files:
        raise ValueError(f"No result files found in {run_dir}")

    converter = NormalizedConverter(out_dir, reference_path)
    written = []
    for src in files:
        written.append(converter.convert_one(src))
        if delete_source:
            src.unlink()
    return written
