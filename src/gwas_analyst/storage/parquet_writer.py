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
    """Read a tab-delimited file and cast columns to the target Arrow schema."""
    import pyarrow.csv as pcsv

    table = pcsv.read_csv(
        path,
        parse_options=pcsv.ParseOptions(delimiter="\t"),
        convert_options=pcsv.ConvertOptions(column_types={f.name: f.type for f in schema}),
    )
    return table.select(schema.names)


def write_inline(run_dir: Path, out_dir: Path, delete_source: bool = False) -> list[Path]:
    """Convert each result file to a standalone typed, compressed Parquet file."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    written = []
    for src in iter_result_files(run_dir):
        table = _read_typed_table(src, ARROW_SCHEMA)
        dst = out_dir / (src.stem + ".parquet")
        pq.write_table(table, dst, compression=COMPRESSION, compression_level=COMPRESSION_LEVEL)
        written.append(dst)
        if delete_source:
            src.unlink()
    return written


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
    shuffled). The first file's `variant` column becomes the reference table
    (variant_id = row index); each subsequent file is validated against it
    and raises ValueError on mismatch, so a violated assumption fails loudly
    rather than silently corrupting the dataset.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    reference_path = reference_path or (out_dir.parent / "reference" / "variants.parquet")
    reference_path.parent.mkdir(parents=True, exist_ok=True)

    files = iter_result_files(run_dir)
    if not files:
        raise ValueError(f"No result files found in {run_dir}")

    written = []
    reference_variants: pa.Array | None = None

    for i, src in enumerate(files):
        table = _read_typed_table(src, ARROW_SCHEMA)
        variants = table.column("variant")

        if reference_variants is None:
            reference_variants = variants
            ref_table = pa.table(
                {
                    "variant_id": pa.array(range(len(variants)), type=pa.int32()),
                    "variant": variants,
                },
                schema=VARIANT_REFERENCE_SCHEMA,
            )
            pq.write_table(ref_table, reference_path, compression=COMPRESSION,
                           compression_level=COMPRESSION_LEVEL)
        elif not variants.equals(reference_variants):
            raise ValueError(
                f"{src.name}: variant column differs from the reference established by "
                f"{files[0].name} -- the 'identical order' assumption does not hold for "
                "this dataset. Re-run investigate_variant_redundancy.py and use "
                "write_inline (or a set-based normalization) instead."
            )

        fact_table = table.drop_columns(["variant"])
        fact_table = fact_table.add_column(
            0, "variant_id", pa.array(range(len(variants)), type=pa.int32())
        )
        fact_table = fact_table.select(NORMALIZED_FACT_SCHEMA.names)

        dst = out_dir / (src.stem + ".parquet")
        pq.write_table(fact_table, dst, compression=COMPRESSION,
                       compression_level=COMPRESSION_LEVEL)
        written.append(dst)
        if delete_source:
            src.unlink()

    return written
