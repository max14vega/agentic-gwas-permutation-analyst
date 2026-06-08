# PRD: GWAS Permutation Pipeline UI
**Version:** 1.1 — MVP (updated after pyseer docs review)  
**Author:** PM Draft for Claude Code  
**Platform:** Local web app (Python backend + browser UI)  
**Audience:** Microbiology PhD researchers running pyseer-based GWAS permutation analyses

---

## 1. Problem Statement

The current GWAS permutation pipeline is entirely CLI-driven and requires researchers to manually edit hardcoded file paths across multiple scripts, create directories by hand, and manage thousands of intermediate TSV files (up to 10,000 per trait) that consume significant disk space. There is no visibility into run progress, no error recovery, and no audit trail. A researcher running 3 traits × 4 variant types × 10,000 permutations generates up to 120,000 intermediate files before any results exist.

**Primary pain points (from user):**
- Storage: each permutation generates a separate TSV file on disk
- Time: no parallelism control or job batching from the UI; LMM similarity matrix is recomputed on every run even though it's identical across all permutations of the same population
- UX: every run requires manually editing shell scripts and Python calls

---

## 2. Goals

| Goal | Success Metric |
|------|----------------|
| Eliminate manual path editing | Zero hardcoded paths in any user-facing step |
| Reduce intermediate storage footprint | Permuted phenotypes generated in-memory, written to a single temp file per permutation, deleted immediately after pyseer exits |
| Cache LMM decomposition across permutations | `--save-lmm` on first permutation, `--load-lmm` on all subsequent ones — eliminating redundant similarity matrix processing |
| Surface run progress | Live progress bar per permutation batch |
| Make the tool usable by non-CLI researchers | A lab member with no terminal experience can complete a full run |

**Out of scope for MVP:**
- Cloud execution or HPC job scheduling
- Multi-user / shared lab accounts
- Results visualization or Manhattan plots (p-value threshold calculation, Manhattan plots)
- pyseer installation management

---

## 3. Users

**Primary user:** Michael (and similar PhD researchers in the lab) — comfortable with Python and the terminal but wants to stop babysitting shell loops and managing file sprawl.

**Secondary user:** Lab members without CLI experience who need to rerun analyses with different traits or variant types.

---

## 4. Architecture Overview

```
[Browser UI] ←→ [Local FastAPI backend (Python)] ←→ [Filesystem + pyseer subprocess]
```

- **Frontend:** Single-page app (React, single HTML file, CDN imports — no build step)
- **Backend:** FastAPI (Python), chosen because the existing scripts are already Python and pyseer runs as a subprocess
- **Storage strategy:** Each permuted phenotype is generated in-memory, written to a single `tempfile.NamedTemporaryFile`, passed to pyseer via `--phenotypes`, then deleted immediately when pyseer exits. At any point in time, exactly 1 intermediate file exists on disk instead of 10,000.
- **No database required for MVP:** job state lives in memory; a JSON summary is written per run on completion

---

## 5. Pyseer Integration — Critical Details

> These come directly from the pyseer documentation and must be followed exactly.

### 5.1 Phenotype File Format
pyseer's `--phenotypes` flag requires a **file path** — it does not support stdin. The phenotype file must be tab-delimited with a **header row as the first line**. The existing `make_permutated_pheno.py` writes files with no header, so the backend must inject a header (`sample\tphenotype\n`) before writing each temp file.

```
sample  phenotype
sample_1    1.23
sample_2    0.87
```

### 5.2 LMM Caching (Major Time Win)
pyseer supports `--save-lmm <prefix>` to save the LMM similarity matrix decomposition after the first permutation, and `--load-lmm <prefix>` to reuse it on all subsequent permutations. This decomposition is computationally expensive and **identical across all permutations** for a given population and similarity matrix. The backend must implement this automatically:

- Permutation 1: run with `--similarity <path> --save-lmm <cache_prefix>`
- Permutations 2–N: run with `--load-lmm <cache_prefix>` (omit `--similarity`)
- Cache stored at `BASEDIR/lmm_cache/{population}/` and reused across runs with the same config

### 5.3 Variant Type Flags
| Variant Type | pyseer Flag | Notes |
|---|---|---|
| Unitigs / k-mers | `--kmers` | Assumes gzipped; add `--uncompressed` if not |
| Panfeed clusters | `--pres` | hashes_to_patterns.tsv format |
| SNPs / INDELs | `--vcf` | Only PASS sites processed; must be bgzipped + indexed |
| Gene PA | `--pres` | gene_presence_absence.Rtab from roary/piggy |

### 5.4 Full pyseer Command Template
```bash
pyseer \
  --phenotypes <temp_pheno_file_with_header> \
  --lmm \
  --kmers <variant_file> \          # or --vcf / --pres
  --similarity <sim_matrix> \       # permutation 1 only
  --save-lmm <lmm_cache_prefix> \   # permutation 1 only
  --load-lmm <lmm_cache_prefix> \   # permutations 2-N
  --min-af 0.05 \
  --max-af 0.95 \
  --cpu <n> \
  > <output_file> \
  2> <log_file>
```

---

## 6. Features (MVP Scope)

### Feature 1 — Project Configuration
The user defines a project once. All subsequent runs reference it.

**Fields:**
- Base directory (path input)
- Master phenotype TSV path
- Sample ID column name (default: `ID`)
- Similarity table path (recombination-corrected, for `--similarity`)
- Population label (e.g. `ALL_CA`, `ST2_Only`) — used for output directory naming and LMM cache namespacing
- Conda environment name (default: `pyseer_env`)
- CPU count for pyseer `--cpu` flag (default: 10)
- Min allele frequency (default: 0.05)
- Max allele frequency (default: 0.95)

Config is saved as `gwas_config.json` in the base directory and auto-loaded on next open.

---

### Feature 2 — Trait Selector
Reads the master phenotype TSV header and presents available trait columns as a checklist. User selects one or more traits to run. Replaces `make_trait_pheno.py` manual invocation.

**Behavior:**
- On TSV load, auto-detect and display all non-ID columns
- User checks traits to include in the current run
- Preview shows N samples with valid (non-null, non-NA) numeric values per trait

---

### Feature 3 — Variant Type Configuration
User selects one or more variant types and provides the file path for each. Maps directly to pyseer flags.

The UI shows a simple table: one row per variant type, with a checkbox to include it and a path field. User can configure multiple variant types; the backend generates one pyseer job per trait × variant type combination.

---

### Feature 4 — Permutation Settings
- Number of permutations (default: 10,000)
- Random seed (default: 84)
- **LMM cache:** checkbox to reuse existing cache if present (default: on). Displays cache status (hit/miss) per population after config is saved.
- **Debug mode:** checkbox to write all permuted TSV files to disk instead of using temp files (default: off)

---

### Feature 5 — Run Dashboard
Replaces the manual shell `for` loop.

**Displays:**
- Job queue table: rows of trait × variant type, with status (queued / running / done / error)
- Active job: progress bar showing `Permutation N of 10,000` + estimated time remaining (rolling average)
- LMM cache indicator: shows whether cache was hit or created for the active job
- Live log tail: last 20 lines of pyseer stderr for the active permutation
- Stop button: sends SIGTERM to active pyseer subprocess and marks job as stopped

**On completion:**
- Results written to `BASEDIR/permutations/{population}/{trait}_{varianttype}/{perm_NNN}.txt`
- Logs written to `BASEDIR/perm_logs/{population}/{trait}_{varianttype}/{perm_NNN}.log`
- LMM cache written to `BASEDIR/lmm_cache/{population}/`
- Run summary JSON written to `BASEDIR/run_summary_{timestamp}.json`

---

### Feature 6 — Run History
Simple log view listing past runs with: timestamp, population, traits run, variant types, permutation count, completion status, and a link to the output directory. Sourced from run summary JSONs in the base directory.

---

## 7. Storage Strategy

The biggest storage win is eliminating intermediate permuted phenotype file accumulation. The backend implements this as follows:

```python
import tempfile, os, subprocess, numpy as np, pandas as pd

def run_single_permutation(pheno_df, perm_index, seed, variant_config, output_path, log_path, lmm_cache_prefix, is_first):
    rng = np.random.default_rng(seed + perm_index)
    perm_df = pheno_df.copy()
    perm_df["phenotype"] = rng.permutation(perm_df["phenotype"].values)

    with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as f:
        f.write("sample\tphenotype\n")  # pyseer requires header
        perm_df.to_csv(f, sep="\t", header=False, index=False)
        tmp_path = f.name

    try:
        cmd = build_pyseer_cmd(tmp_path, variant_config, lmm_cache_prefix, is_first)
        with open(output_path, "w") as out, open(log_path, "w") as log:
            subprocess.run(cmd, stdout=out, stderr=log, check=True)
    finally:
        os.unlink(tmp_path)  # always delete, even on error
```

The LMM cache (`--save-lmm` / `--load-lmm`) removes the most expensive repeated computation across the 10,000 permutations.

---

## 8. Tech Stack

| Layer | Choice | Reason |
|---|---|---|
| Backend | Python + FastAPI | Matches existing script language; easy subprocess management |
| Frontend | React (single HTML file, CDN imports) | No build step; runs anywhere with a browser |
| Real-time updates | Server-Sent Events (SSE) | Progress streaming without WebSocket complexity |
| Config persistence | JSON file in project base dir | No database dependency; portable with the project |
| Process management | Python `subprocess.run` / `Popen` | Direct control over pyseer child processes |

---

## 9. File & Folder Structure (App)

```
gwas-pipeline-ui/
├── backend/
│   ├── main.py              # FastAPI app, all API routes
│   ├── pipeline.py          # Core: permutation loop, pyseer subprocess, LMM cache logic
│   ├── config.py            # Config load/save (gwas_config.json)
│   └── requirements.txt     # fastapi, uvicorn, numpy, pandas
├── frontend/
│   └── index.html           # Single-file React app (CDN imports)
├── make_trait_pheno.py      # Original script (preserved, not called by backend — logic inlined)
├── make_permutated_pheno.py # Original script (preserved, not called by backend — logic inlined)
└── README.md
```

---

## 10. API Routes

| Method | Route | Description |
|---|---|---|
| GET | `/config` | Load saved project config |
| POST | `/config` | Save project config |
| POST | `/traits` | Parse master TSV → return trait columns + valid sample counts |
| GET | `/lmm-cache` | Return cache status per population (hit/miss/stale) |
| POST | `/run` | Start a pipeline run; returns `run_id` |
| GET | `/run/{id}/progress` | SSE stream: `{perm: N, total: M, job: "trait_vartype", eta_seconds: X}` |
| GET | `/run/{id}/logs` | SSE stream: last 20 lines of active pyseer stderr |
| POST | `/run/{id}/stop` | Gracefully stop active run |
| GET | `/history` | List past run summaries from JSON files in base dir |

---

## 11. MVP User Flow

1. Open browser to `localhost:8000`
2. Set base directory and config fields → save config
3. Load master phenotype TSV → trait columns appear as checkboxes; valid sample count shown per trait
4. Select traits, variant type(s) + file paths, permutation count and seed
5. Click "Run" → job queue table populates
6. Watch per-job progress bars + live pyseer stderr tail
7. LMM cache created on first permutation, reused for all subsequent ones automatically
8. Run completes → results in configured output directories
9. View run history for past jobs

---

## 12. Out of Scope (Post-MVP Backlog)

- HPC / SLURM job submission
- p-value threshold calculation via `count_patterns.py`
- Results visualization (Manhattan plots, Q-Q plots)
- Multi-user lab server deployment with auth
- Automatic pyseer/conda environment setup
- Email/Slack notifications on job completion

---

## 13. Instructions for Claude Code

**Build order:**

1. Scaffold `gwas-pipeline-ui/` with the folder structure above
2. Install: `pip install fastapi uvicorn numpy pandas`
3. **Implement `pipeline.py` first.** Key function signature:
   ```python
   def run_permutation_batch(
       pheno_df: pd.DataFrame,       # pre-loaded, pre-filtered trait DataFrame (2 cols: sample, phenotype)
       n_perms: int,
       seed: int,
       variant_config: dict,         # {type: "kmers"|"vcf"|"pres", path: str, uncompressed: bool}
       output_dir: Path,
       log_dir: Path,
       lmm_cache_prefix: Path,       # path prefix for --save-lmm / --load-lmm
       similarity_path: Path,        # used only for permutation 1
       pyseer_env: str,              # conda env name
       cpu: int,
       min_af: float,
       max_af: float,
       progress_callback: callable   # called with (current_perm, total_perms) after each perm
   )
   ```
4. Implement FastAPI routes in `main.py` with SSE for `/run/{id}/progress` and `/run/{id}/logs`
5. Build the React frontend in `frontend/index.html` — must include: config form, trait checklist with sample counts, variant type rows with path inputs, run button, job queue table, progress bars, log tail panel, history table
6. Test end-to-end with n=3 permutations and a dummy 5-row phenotype TSV before wiring real pyseer

**Critical constraints:**
- pyseer's `--phenotypes` requires a **file path, not stdin** — always write to a temp file with a header row
- Always inject `sample\tphenotype` as the header line in temp phenotype files (the original scripts omit this)
- LMM cache: use `--save-lmm` on permutation 1, `--load-lmm` on permutations 2–N for the same population + similarity matrix. Cache is keyed by `{population}` subdirectory under `BASEDIR/lmm_cache/`
- Temp files must be deleted in a `finally` block — never leave them on error
- Activate the conda environment before calling pyseer: `conda run -n {env} pyseer ...` or equivalent
- No hardcoded paths anywhere in the codebase
