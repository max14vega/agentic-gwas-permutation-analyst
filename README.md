# agentic-gwas-permutation-analyst

[WIP] Storage efficient GWAS permutation application with agentic analyst layer for rapid analysis, built with claude code

## Layout

- `src/gwas_analyst/storage/schema.py` -- the 12-column result schema (7 columns
  confirmed, 5 placeholders pending the algorithm spec -- see comments there).
- `src/gwas_analyst/storage/txt_reader.py` -- reads the algorithm's raw `.txt` output.
- `src/gwas_analyst/storage/parquet_writer.py` -- converts raw output to compact,
  typed, zstd-compressed Parquet, either inline or normalized (variant sequences
  split into a shared reference table -- see below).
- `src/gwas_analyst/storage/db.py` -- DuckDB query layer over the Parquet dataset
  (the "database baseline": canonical queries + a guarded read-only SQL executor).
- `src/gwas_analyst/agent/analyst.py` -- the agentic analyst (the "AI baseline"):
  Claude + a read-only SQL tool + a system prompt encoding GWAS/microbiology domain
  knowledge.
- `src/gwas_analyst/data/synthetic.py` -- generates synthetic permutation output
  shaped like the real algorithm's, for development without the real data.
- `scripts/investigate_variant_redundancy.py` -- run this against the **real**
  output first. It checks whether the same variant set is reused across
  permutation files (expected for standard permutation testing) and recommends
  `write_normalized` vs `write_inline` accordingly.
- `scripts/ask.py` -- CLI/REPL for the agent (requires `ANTHROPIC_API_KEY`).

## Quickstart

```bash
python3 -m venv .venv && .venv/bin/pip install -e .

# Generate synthetic data shaped like the real output (until it's available)
.venv/bin/python -m gwas_analyst.data.synthetic data/run_fixed --num-files 8 --num-rows 4000

# 1. Check whether variant sets repeat across files (decides the storage layout)
.venv/bin/python scripts/investigate_variant_redundancy.py data/run_fixed

# 2a. If variants repeat (normalized layout -- the expected case):
.venv/bin/python -c "
from pathlib import Path
from gwas_analyst.storage.parquet_writer import write_normalized
write_normalized(Path('data/run_fixed'), Path('data/parquet/permutations'))
"
# -> writes data/parquet/permutations/*.parquet and data/parquet/reference/variants.parquet

# 2b. If variants differ per file, use write_inline instead (same call shape).

# 3. Ask the analyst questions


## Storage approach

Raw `.txt` output (~300-400 GB for a full run: 10,000 files x ~400k rows x 12
columns) is converted to typed, zstd-compressed **Parquet**. Two layouts are
supported -- `investigate_variant_redundancy.py` tells you which applies to your
data:

- **Normalized** (expected for standard permutation testing, where the same
  variant set is re-tested under shuffled phenotypes): the `variant` DNA
  sequence column -- almost certainly the largest -- is stored exactly once in
  `reference/variants.parquet`, and each permutation's results carry only an
  integer `variant_id`. This removes ~10,000x duplication of the largest column.
- **Inline**: each file is converted to a standalone typed Parquet file, relying
  on Parquet's columnar typing, dictionary encoding, and zstd compression.

Either way, **DuckDB** queries the resulting Parquet dataset directly (no server,
no separate database to manage -- a good fit for a single-user setup), and the
agent layer issues guarded read-only SQL through that same connection.
</content>
