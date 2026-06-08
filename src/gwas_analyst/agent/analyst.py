"""The agentic analyst -- the "AI baseline".

A thin loop around the Anthropic Messages API: the model is given one tool
(read-only SQL against the DuckDB view `results`/`permutations`/`variants`)
plus a system prompt encoding GWAS-permutation and microbiology domain
knowledge, and translates natural-language questions into queries + narrative
interpretation.

No UI here by design (per the build order: get this working via CLI first,
wrap it in FastAPI + React once it's solid).
"""

import json
from typing import Any

import anthropic
import duckdb

from gwas_analyst.storage import db
from gwas_analyst.storage.schema import ALL_COLUMNS

MODEL = "claude-sonnet-4-6"

_SCHEMA_DESCRIPTION = "\n".join(
    f"  - {name} ({dtype}): {desc}" + ("  [PLACEHOLDER -- name/type TBD]" if name.startswith("placeholder_") else "")
    for name, dtype, desc in ALL_COLUMNS
)

SYSTEM_PROMPT = f"""You are a domain-expert analyst for GWAS (genome-wide association
study) permutation testing results in a microbiology research context.

## Domain knowledge you should apply

- **Permutation testing**: the phenotype/label is repeatedly shuffled and the
  association test is re-run on each shuffle, building an *empirical null
  distribution* for the test statistic. A variant's observed p-value is only
  meaningful when compared against this null -- e.g. an empirical permutation
  p-value is the fraction of permuted runs whose statistic is at least as
  extreme as the observed one. Always frame "significance" claims this way
  rather than relying on nominal p-value cutoffs alone.
- **Multiple testing**: with ~400k variants tested per run, raw p-values must
  be corrected (FDR / Benjamini-Hochberg, or permutation-derived thresholds)
  before calling anything "significant". Mention this when a user asks about
  "significant" or "top" hits.
- **beta / beta_std_err**: `beta` is the estimated effect size of a variant on
  the trait; `beta_std_err` is its precision. A large beta with a large
  std_err is much weaker evidence than the same beta with a small std_err
  (consider beta / beta_std_err, an approximate test statistic, when comparing
  variants).
- **variant_h2**: the heritability contribution attributed to a variant --
  i.e. how much trait variance it explains. Useful for ranking biological
  relevance independent of statistical significance.
- **lct_pvalue / filter_pvalue**: two stages of association testing (an
  initial filter and a focused locus-level test, named for the lactase (LCT)
  gene region -- a classic locus for studying gene-environment interaction
  and recent human/microbiome adaptation, e.g. lactase persistence and gut
  microbiota composition). When a user asks about "the LCT result" they
  almost certainly mean `lct_pvalue`.
- **variant**: the raw DNA sequence under test.

## Data schema

You query a DuckDB view called `results` (per-row permutation output, with
the variant sequence joined back in) and `permutations` (same data plus a
`source_file` column identifying which permutation run produced each row;
in the normalized layout this excludes `variant` -- join `variants` on
`variant_id` if you need the sequence). Columns:

{_SCHEMA_DESCRIPTION}

Five of the twelve columns above are placeholders pending the final algorithm
spec -- if a query needs one of them and it doesn't exist in the live schema,
say so plainly rather than guessing at a name.

## How to answer

1. Use the `query_results` tool to run a read-only SQL query (SELECT/WITH
   only) against the dataset to ground your answer in the actual data --
   don't speculate about numbers you haven't queried.
2. Then explain the result in plain language using the domain framing above
   (null distributions, multiple-testing correction, effect size vs.
   precision, heritability) -- a microbiologist should come away understanding
   both *what* the numbers say and *how much weight* to put on them.
3. Be concise. Lead with the answer, support it with the numbers you queried.
"""

QUERY_TOOL = {
    "name": "query_results",
    "description": (
        "Run a read-only SQL query (SELECT or WITH only) against the permutation "
        "results dataset. Available views: `results` (full per-row data including "
        "the variant sequence), `permutations` (per-row stats + source_file, "
        "without the joined sequence), `variants` (variant_id -> sequence, "
        "normalized layouts only). Returns up to 200 rows as JSON."
    ),
    "input_schema": {
        "type": "object",
        "properties": {
            "sql": {"type": "string", "description": "A SELECT/WITH SQL query"},
        },
        "required": ["sql"],
    },
}

_MAX_ROWS = 200


class Analyst:
    """Wraps an Anthropic client + DuckDB connection in a tool-use loop."""

    def __init__(self, con: duckdb.DuckDBPyConnection, client: anthropic.Anthropic | None = None):
        self.con = con
        self.client = client or anthropic.Anthropic()

    def _run_tool(self, name: str, tool_input: dict[str, Any]) -> str:
        if name != "query_results":
            return json.dumps({"error": f"unknown tool {name!r}"})
        sql = tool_input.get("sql", "")
        try:
            rows = db.run_readonly_query(self.con, sql)
        except Exception as exc:  # noqa: BLE001 -- surfaced to the model, not the user
            return json.dumps({"error": str(exc)})
        truncated = rows[:_MAX_ROWS]
        payload: dict[str, Any] = {"row_count": len(rows), "rows": truncated}
        if len(rows) > _MAX_ROWS:
            payload["note"] = f"truncated to first {_MAX_ROWS} of {len(rows)} rows"
        return json.dumps(payload, default=str)

    def ask(self, question: str, max_turns: int = 6) -> str:
        """Send `question` through the tool-use loop and return the final text answer."""
        messages: list[dict[str, Any]] = [{"role": "user", "content": question}]

        for _ in range(max_turns):
            response = self.client.messages.create(
                model=MODEL,
                max_tokens=2048,
                system=SYSTEM_PROMPT,
                tools=[QUERY_TOOL],
                messages=messages,
            )

            if response.stop_reason != "tool_use":
                return "".join(
                    block.text for block in response.content if block.type == "text"
                )

            messages.append({"role": "assistant", "content": response.content})
            tool_results = []
            for block in response.content:
                if block.type != "tool_use":
                    continue
                result = self._run_tool(block.name, block.input)
                tool_results.append(
                    {"type": "tool_result", "tool_use_id": block.id, "content": result}
                )
            messages.append({"role": "user", "content": tool_results})

        return "Reached the tool-use turn limit without a final answer -- try a narrower question."
