"""Column schema for permutation result rows.

Cross-checked against the real pyseer LMM-mode output header
(reference/pyseer/pyseer/__main__.py:488-513) and value formatting
(reference/pyseer/pyseer/utils.py:39-105), traced against the *exact* flags
the user's pipeline passes (`--lmm --min-af 0.05 --max-af 0.95 --cpu N`, with
no `--lineage`, `--print-samples`, or `--wg`):

    variant  af  filter-pvalue  lrt-pvalue  beta  beta-std-err  variant_h2  notes

That header-construction logic only appends `lineage` when `--lineage` (or
`--wg`+`--sequence-reweighting`+`--lineage-clusters`) is passed, and only
appends `k-samples`/`nk-samples` when `--print-samples` is passed -- neither
applies here, so those columns do not exist in the real output.

The permutation identity is *not* a data column at all: the pipeline's bash
loop names each output file after its phenotype file
(`base=$(basename "$pheno" .tsv)` -> `${base}.txt`, e.g. `perm_0001.txt`).
`storage.db.connect` already extracts this as `source_file` from the Parquet
filename, so no placeholder/bookkeeping column is needed in the schema.

Note pyseer formats numeric fields as scientific-notation strings ('%.2E') or
empty string when non-finite -- the txt reader/writer must treat empty strings
as null when casting to float.
"""

import pyarrow as pa

# `lct_pvalue` in the user's original description was a mishearing of
# `lrt-pvalue` (likelihood ratio test) -- pyseer's actual column name.
ALL_COLUMNS = [
    ("variant", pa.string(),
     "k-mer DNA sequence string (or VCF/Rtab variant ID, depending on input type) under test"),
    ("af", pa.float32(),
     "allele frequency of the variant"),
    ("filter_pvalue", pa.float32(),
     "p-value from the pre-filtering association test"),
    ("lrt_pvalue", pa.float32(),
     "likelihood-ratio-test p-value from the full association model"),
    ("beta", pa.float32(),
     "estimated effect size of the variant on the trait"),
    ("beta_std_err", pa.float32(),
     "standard error of the beta estimate"),
    ("variant_h2", pa.float32(),
     "fraction of model heritability (h2) explained by the variant (LMM mode)"),
    ("notes", pa.string(),
     "comma-separated flags from a fixed vocabulary, e.g. af-filter, bad-chisq, "
     "pre-filtering-failed, high-bse, perfectly-separable-data, matrix-inversion-error, "
     "firth-fail, missing-data-error, lrt-filtering-failed -- mostly empty"),
]

assert len(ALL_COLUMNS) == 8, "Schema must have exactly 8 columns"

ARROW_SCHEMA = pa.schema([(name, dtype) for name, dtype, _ in ALL_COLUMNS])

COLUMN_NAMES = [name for name, _, _ in ALL_COLUMNS]

# Reference-table schema used when variant identities are normalized out of the
# per-permutation fact tables (see investigate_variant_redundancy.py).
VARIANT_REFERENCE_SCHEMA = pa.schema([
    ("variant_id", pa.int32()),
    ("variant", pa.string()),
])

# Fact-table schema for the normalized layout: same as ALL_COLUMNS but with
# `variant` replaced by an integer foreign key into the reference table.
NORMALIZED_FACT_SCHEMA = pa.schema(
    [("variant_id", pa.int32())]
    + [(name, dtype) for name, dtype, _ in ALL_COLUMNS if name != "variant"]
)
