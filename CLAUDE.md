# CLAUDE.md — GSE Router Project Guidance

## Project Context

This is a **smart routing proxy** for converting arbitrary GEO single-cell RNA-seq datasets into standardized `.h5ad` objects. The sole purpose: take a GEO accession ID, output a standardized AnnData file with raw counts preserved and full provenance tracked.

The core problem: GEO data is submitter-supplied and chaotically organized — formats vary wildly across datasets. There is no single tool that handles all cases. This project fills that gap with an intelligent routing layer.

## Current State (2026-07-07)

### v1→v2 Intelligence Port — COMPLETE

Six phases of biological intelligence were ported from Aging_OS v1 into the GSE Router v2, all additive (no existing logic replaced):

| Phase | Feature | File(s) | Tests | Status |
|-------|---------|---------|-------|--------|
| 1 | **Shape-Matching MEX Fallback** — reads MTX dimensions, matches barcode/feature files by line count | `router/classifier.py` | 6 new | ✅ Done |
| 2 | **Gene Symbol Washing** — strips genome prefixes, splits Ensembl concatenations, remaps numeric indices | `tools/converters.py` (all 10 converters) | 7 new | ✅ Done |
| 3 | **Recursive Archive Extraction** — re-extracts nested archives, exempts 10x compressed files | `probes/file_probe.py` | 6 new | ✅ Done |
| 4 | **Gene-Case Species Cross-Validation** — detects species from gene casing (UPPER=Human, Title=Mouse) | `probes/metadata.py`, `pipeline/orchestrator.py` | 8 new | ✅ Done |
| 5 | **Graceful Degradation Chains** — falls back to lower-priority routes on primary conversion failure | `pipeline/orchestrator.py` | — | ✅ Done |
| 6 | **LLM Condition Inference** — optional DeepSeek-powered control/treatment classification | `probes/metadata.py` | — | ✅ Done |

### Phase 6 — Upstream Polish — COMPLETE (2026-07-07)

Three HIGH-priority improvements implemented:

| # | Feature | What changed | Impact |
|---|---------|-------------|--------|
| 1 | **Archived 10x Bundle Recognition** | Added 8 more filename patterns to `file_probe.py` and `classifier.py` (`barcodes.tsv.tar.gz`, `cellranger_count.tar.gz`, `outs.tar.gz`, etc.). `ARCHIVE_10X` already existed — expanded coverage. | ~3-5% success rate gain |
| 2 | **Non-Count Supplementary File Filtering** | Reordered `_match_filename` to check BCR/TCR/peak/diff patterns BEFORE generic CSV/TSV fallback. Added GFF/GTF detection in `_classify_gzip`. Added spatial transcriptomics detection. | Prevents ~2% misrouting |
| 3 | **Technology Detection from Supplementary Filenames** | Orchestrator now collects supplementary filenames from GEO FTP listing and passes them to `probe_all()` → `_infer_technology()` Layer 3. Manifest now also propagates `tissue_or_organ`, `disease`, `condition_groups`, `publication`, `metadata_confidence` from study context. | Closes ~2-3% detection gap |
| — | **Test Suite** | Created `tests/test_file_probe.py` (32 tests), `tests/test_classifier.py` (22 tests), `tests/test_route_planner.py` (37 tests). All 91 tests passing. | Regression safety |

### Success Rate (post-improvements, projected)

| Outcome | Before | After | Description |
|---------|--------|-------|-------------|
| **Clean success** | ~35% | **~42-45%** | h5ad with raw counts in `layers["counts"]` |
| **Partial success** | ~48% | **~44-47%** | h5ad with caveats |
| **Hard failure** | ~17% | **~10-12%** | No path to h5ad |
| **No data** | ~5% | ~5% | Supplementary not downloadable |

**Projected usable rate: ~87-89%** (up from ~83%)

## Key Design Rules

1. **Manifest-driven, not folder-driven.** The agent must never reason directly from directory listings. Every file gets probed and documented in a JSON sidecar first. All downstream decisions are based on these structured manifests.

2. **Three intelligence layers (upstream only):**
   - **Database metadata probe** — queries GEO/GSM/SRA/publication metadata (species, tissue, disease, platform, sample groups)
   - **File probe** — inspects local files, writes per-file JSON sidecars with detected type, confidence, paired files
   - **Route planner** — combines both evidence sources, applies biological judgment, selects processing method

3. **Separation of concerns — the MCP+Skills decoupling:**
   - **Skills contain judgment, policy, and biological knowledge** (which route to choose, what QC thresholds to apply, how to analyze)
   - **MCP tools execute deterministic operations** (convert formats, validate schemas, merge metadata, run analysis)
   - The planner decides; the tools execute. This decoupling ensures:
     - **Change isolation**: updating a routing rule doesn't touch conversion code
     - **Smartness**: skills can be upgraded (better heuristics, better LLM prompts) without touching tools
     - **Testability**: tools are pure functions; skills are policy-tested separately

4. **Provenance is first-class.** Every processing step is recorded. Raw counts are preserved in `AnnData.layers["counts"]`. Source file hashes are tracked.

5. **Fail safely, not silently.** When confidence is low or files are ambiguous, the system stops and reports for manual review. It must never guess.

6. **All logic lives in MCP tools or skills — nothing ad-hoc.** The orchestrator, classifier, converters, and probes are all exposed through the MCP boundary. The skill encodes policy; tools execute deterministically.

## Routing Rules (biological judgment)

### Cardinal rules (from GEO examples research)

1. **RAW.tar ≠ FASTQ.** Many GEO `RAW.tar` files contain processed matrices, `.loom`, `.h5ad`, or metadata — not necessarily FASTQ. Always check `filelist.txt` or FTP index before downloading.

2. **Skip FASTQ if processed data exists.** If ANY usable processed count/object file is available (priority 1-6 below), do NOT run FASTQ/SRA processing by default. Only fall back to FASTQ if ALL higher-priority routes fail.

3. **Never treat normalized as raw counts.** TPM, FPKM, log-normalized, and other normalized data go in `layers["normalized"]`, NOT `layers["counts"]`. Mark output as `limited_use = True`.

### Explicit 8-level routing priority

When multiple data formats exist in a dataset, prefer in this exact order:

| Priority | Format | Route | Feasibility |
|----------|--------|-------|-------------|
| 1 | `.h5ad` | validate/standardize | High |
| 2 | 10x `.h5` | `read_10x_h5` | High |
| 3 | `.loom` / `.loom.gz` | `read_loom` | High |
| 4 | 10x MEX triplet | `read_10x_mtx` | High |
| 5 | raw UMI TXT/CSV/RDS matrix | parse/convert (confirmed raw counts) | High |
| 6 | Seurat/SCE `.rds` / `.h5seurat` | R bridge conversion | Medium |
| 7 | normalized matrix only | mark limited-use; do NOT treat as counts | High |
| 8 | FASTQ/SRA | **avoid** unless no processed data exists | Medium-Hard |

### Graceful degradation (v1→v2 Phase 5)

When a primary conversion method fails, the orchestrator tries fallback routes automatically:

| Primary | Fallback(s) |
|---------|-------------|
| `READ_10X_MTX` | `READ_CSV_MATRIX` |
| `READ_10X_H5` | `READ_H5AD` |
| `READ_LOOM` | `READ_CSV_MATRIX` |
| `READ_RDS_SEURAT` ↔ `READ_RDS_SCE` | (cross-fallback) |
| `READ_RAW_UMI_MATRIX` | `READ_CSV_MATRIX` |

### Standard rules

- Species → reference genome: GRCh38 (human), GRCm39 (mouse), custom otherwise
- Protocol determines pipeline: 10x → Cell Ranger/STARsolo, Smart-seq2 → different route
- Control vs experiment → obs condition labels
- Preserve raw counts in `AnnData.layers["counts"]` when available
- Attach sample/cell metadata to `obs`
- Attach gene metadata to `var`
- If confidence < 0.7, stop and report instead of guessing
- **Gene symbols are washed** after every conversion (strip genome prefixes, split Ensembl concatenations)
- **Species is cross-validated** post-conversion (gene casing vs GEO metadata)

### MCP tool boundary

The agent (skill) decides routes; MCP tools execute deterministically:

```
fetch_geo_metadata(gse_id)       -> StudyContext JSON
list_geo_supplementary(gse_id)   -> File listing with sizes/types
read_geo_filelist(gse_id)        -> Parsed filelist.txt contents
probe_local_directory(path)      -> FileSidecar[] JSON
write_file_sidecars(path)        -> Per-file JSON sidecars written
classify_scrna_files(manifest)   -> file classification
group_related_files(manifest)    -> FileGroup[]
build_manifest(sidecars, ctx)    -> DatasetManifest JSON
plan_route(manifest)             -> RoutePlan JSON
convert_to_h5ad(group, method)   -> .h5ad + provenance
validate_h5ad(h5ad)              -> QCResult
write_processing_report(...)     -> ProcessingReport JSON
```

### Environment management (layered, not monolithic)

Use **separate micromamba/conda environments** — do NOT install everything into one:

| Environment | Contents | When needed |
|-------------|----------|-------------|
| `router-core` | pydantic, httpx, pandas, click, structlog, jsonschema | Always |
| `scrna-python` | scanpy, anndata, scipy, numpy, h5py, loompy | Most datasets |
| `r-converters` | R, Seurat, SeuratDisk, SingleCellExperiment, zellkonverter | Only when `.rds` appears |
| `raw-pipeline` | nextflow, nf-core/fetchngs, nf-core/scrnaseq, STARsolo | Last resort; disabled by default |

The MCP server reads `config/tools.yaml` to dispatch each tool call to the right environment. On Windows, run bioinformatics environments inside **WSL2 Ubuntu** — native Windows is fine for router-core and scrna-python, but R/Seurat, Nextflow, and Cell Ranger need Linux.

### Logging architecture (two layers)

**Layer 1 — MCP tools:** `data/GSE123456/logs/<run_id>/`
- `tool_calls.jsonl` — every MCP invocation: inputs, outputs, timing, errors
- `probe_events.jsonl` — file probing results with confidence and evidence
- `route_decisions.jsonl` — route selection rationale, rules triggered
- `conversion.log` — format conversion output and warnings
- `validation.json` — QC metrics and pass/fail

**Layer 2 — Claude Code hooks:** agent-level events via hooks on `SessionStart`, `PreToolUse`, `PostToolUse`, `Stop`, `SessionEnd`. Hook input includes `transcript_path` for full traceability.

### Development phases

```
Phase 1: project skill + simple MCP server                  ✅ DONE
Phase 2: schema + sidecar + manifest + structured logging   ✅ DONE
Phase 3: processed-file conversion (h5ad, 10x, loom, MEX, CSV)  ✅ DONE
Phase 3.5: v1 intelligence port (6 sub-phases)              ✅ DONE
Phase 4: RDS bridge (Seurat/SCE → h5ad)                     ⚠️ R scripts exist; env setup incomplete
Phase 5: optional FASTQ/SRA route (disabled by default)     ⚠️ Stubbed; not operational
Phase 6: upstream polish (HIGH-priority gaps)               📋 NEXT
```

### Development phases

```
Phase 1: project skill + simple MCP server                  ✅ DONE
Phase 2: schema + sidecar + manifest + structured logging   ✅ DONE
Phase 3: processed-file conversion (h5ad, 10x, loom, MEX, CSV, .mat)  ✅ DONE
Phase 3.5: v1 intelligence port (6 sub-phases)              ✅ DONE
Phase 4: RDS bridge (Seurat/SCE → h5ad)                     ⚠️ R scripts exist; env setup incomplete
Phase 5: optional FASTQ/SRA route (disabled by default)     ⚠️ Stubbed; not operational
Phase 6: upstream polish (archived 10x, non-count filtering, tech detection)  ✅ DONE
```

## Test Suite

```bash
# Run all tests (91 tests, ~0.6s)
pytest tests/ -v

# Test categories:
#   tests/test_file_probe.py      — format detection, archive handling, non-count filtering
#   tests/test_classifier.py      — grouping logic, ARCHIVE_10X mapping, strategy selection
#   tests/test_route_planner.py   — priority ordering, reference genome, rule application
```

## Technology Stack

- **Python 3.10+** with type hints throughout
- **AnnData / Scanpy** — target format, basic analysis, and downstream backbone
- **GEOparse** — GEO metadata parsing (`router-core` env)
- **Pydantic** — schema validation for all JSON structures (`router-core` env)
- **Click** — CLI interface (`router-core` env)
- **structlog** — structured logging (`router-core` env)
- **httpx** — async HTTP for GEO FTP supplementary file probing
- **micromamba** — environment management for layered tool envs
- **R / Seurat / SeuratDisk / zellkonverter** — RDS conversion (`r-converters` env, Linux/WSL2)

## Code Conventions

- All JSON-serializable data uses Pydantic models defined in `schemas/`
- Probing is read-only and side-effect-free
- Tools are idempotent where possible
- Every conversion tool returns both the object and a provenance record
- Use type hints on all public functions
- Docstrings in Google style
- Tool scripts log environment + package versions on startup
- MCP tools declare their required environment in `config/tools.yaml`
- `ProcessingStep.parameters` is typed `dict[str, str | int | float | bool]` — NEVER put `None` or `list` values in it. Join lists with `"; "`.
- Gene washing is applied after every conversion; original names preserved in `adata.var["_original_name"]`

## Project Skill & MCP

The agent's routing policy lives in `.claude/skills/geo-scrna-router/SKILL.md` with references for routing rules, file signatures, biology rules, and output schemas. The MCP server at `gse_router_mcp/server.py` exposes deterministic tools configured in `.mcp.json`.

## Key Files

| File | Purpose |
|------|---------|
| `src/gse_router/probes/metadata.py` | GEO metadata probe, organism detection (3-layer fallback), technology detection (2-layer), species-from-genes, LLM condition inference |
| `src/gse_router/probes/file_probe.py` | Local file inspection, per-file sidecar generation, recursive archive extraction with format-aware exemption |
| `src/gse_router/schemas/sidecar.py` | `FileDetectedType` enum, `FileSidecar`, `ProcessingStage` |
| `src/gse_router/schemas/manifest.py` | `DatasetManifest`, `FileGroup`, `FileGroupType`, `Strategy` |
| `src/gse_router/schemas/route_plan.py` | `RoutePlan`, `ProcessingStep`, `RouteMethod`, `ReferenceGenome` |
| `src/gse_router/router/classifier.py` | `FileClassifier` — groups files into MEX triplets (with shape-matching fallback), H5, FASTQ pairs, etc. |
| `src/gse_router/router/rules.py` | Routing rules, `decode_route()`, `select_reference_genome()` |
| `src/gse_router/router/route_planner.py` | `RoutePlanner` — builds ordered `RoutePlan` from manifest + context |
| `src/gse_router/tools/converters.py` | Format converters (10 subcommands + gene symbol washing in all converters) |
| `src/gse_router/tools/geo_downloader.py` | GEO supplementary file downloader with filelist.txt parsing |
| `src/gse_router/tools/validators.py` | h5ad validation and QC checks |
| `src/gse_router/tools/metadata_tools.py` | Cell/gene metadata merging into AnnData |
| `src/gse_router/pipeline/orchestrator.py` | Full pipeline orchestrator (probe → plan → execute → validate → report) with graceful degradation and species cross-validation |
| `gse_router_mcp/server.py` | MCP server — JSON-RPC stdio, 27 tools |

## Reference

- `.claude/skills/geo-scrna-router/` — skill definition and reference docs
