"""Column schema for permutation result rows.

Seven columns are confirmed by the user; the remaining five are placeholders
until the permutation algorithm spec is finalized. Update PLACEHOLDER_COLUMNS
in place when the real names/types are known -- nothing else needs to change.
"""

import pyarrow as pa

# Confirmed columns (name, arrow type, description)
CONFIRMED_COLUMNS = [
    ("variant", pa.string(), "DNA sequence string identifying the variant under test"),
    ("filter_pvalue", pa.float32(), "p-value from the pre-filtering association test"),
    ("lct_pvalue", pa.float32(), "p-value from the LCT-locus association test"),
    ("beta", pa.float32(), "estimated effect size of the variant on the trait"),
    ("beta_std_err", pa.float32(), "standard error of the beta estimate"),
    ("variant_h2", pa.float32(), "heritability contribution attributed to the variant"),
    ("notes", pa.string(), "free-text annotation, mostly empty"),
]

# Placeholder columns -- replace with real name/type/description once confirmed.
# Guesses below follow common GWAS permutation output conventions.
PLACEHOLDER_COLUMNS = [
    ("placeholder_perm_index", pa.int32(), "TBD: index of the permutation run"),
    ("placeholder_chromosome", pa.int16(), "TBD: chromosome number"),
    ("placeholder_position", pa.int64(), "TBD: base-pair position on the chromosome"),
    ("placeholder_sample_size", pa.int32(), "TBD: number of samples used in the test"),
    ("placeholder_maf", pa.float32(), "TBD: minor allele frequency"),
]

ALL_COLUMNS = CONFIRMED_COLUMNS + PLACEHOLDER_COLUMNS

assert len(ALL_COLUMNS) == 12, "Schema must have exactly 12 columns"

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
