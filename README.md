# GSE Router — GEO scRNA-seq to .h5ad Converter

Takes a GEO single-cell RNA-seq accession ID (e.g., GSE123456) and outputs a standardized `.h5ad` (AnnData) file with preserved raw counts, cell/gene metadata, and full provenance tracking.

## Quick Start

```bash
# 1. Install (router-core only — lightweight)
pip install -e .

# 2. Set up environment variables
cp .env.example .env
# Edit .env and add your NCBI API key (register at https://www.ncbi.nlm.nih.gov/account/)
# Optional: add DEEPSEEK_API_KEY for LLM-powered control/treatment classification

# 3. Process a dataset (auto-download + convert)
gse-router process GSE115469 --auto-download

# 4. Check the output
gse-router validate output/h5ad/GSE115469_final.h5ad
```

## Commands

| Command | What it does |
|---------|-------------|
| `gse-router process GSE123456 --auto-download` | Full pipeline: download, probe, convert, validate |
| `gse-router process GSE123456 -d ./my_data/` | Process from local data directory |
| `gse-router list-supplementary GSE123456` | List downloadable files with sizes |
| `gse-router download GSE123456` | Download supplementary files only |
| `gse-router probe GSE123456 -d ./data/` | Probe files, write manifest (no conversion) |
| `gse-router validate output/file.h5ad` | Validate an existing .h5ad |

### Process flags

```
--auto-download / --no-auto-download   Download GEO files before processing
--dry-run                              Plan route only, no execution
--skip-metadata                        Skip GEO metadata queries
-d, --data-dir PATH                    Local data directory
-o, --output-dir PATH                  Output directory (default: ./output)
-e, --email EMAIL                      NCBI contact email
--api-key KEY                          NCBI API key
```

## Supported Formats (8-level priority)

| Priority | Format | Automatic? |
|----------|--------|-----------|
| 1 | `.h5ad` | ✅ Validate & standardize |
| 2 | 10x `.h5` | ✅ Convert via scanpy |
| 3 | `.loom` | ✅ Convert via scanpy |
| 4 | 10x MEX triplet (matrix.mtx + barcodes + features) | ✅ Convert via scanpy |
| 5 | Raw UMI TXT/CSV/RDS matrix | ✅ Convert if confirmed counts |
| 6 | Seurat/SCE `.rds` / `.h5seurat` | ⚠️ Needs R environment (WSL2) |
| 7 | Normalized matrix only | ⚠️ Converted, marked limited-use |
| 8 | FASTQ/SRA raw reads | ⚠️ Disabled by default |

The pipeline auto-detects format and picks the best route. If a conversion fails, it falls back to alternative methods automatically.

## Output

```
output/
├── GSE123456_report.json          # Full provenance report
└── h5ad/
    └── GSE123456_final.h5ad       # Standardized AnnData
        ├── .X                      # Expression matrix
        ├── .layers["counts"]       # Raw counts (preserved)
        ├── .obs                    # Cell metadata
        └── .var                    # Gene metadata
```

## Requirements

- **Python 3.10+** with packages: `pydantic`, `click`, `structlog`, `httpx`, `pandas`, `pyyaml`, `jsonschema`
- **scrna-python** env (micromamba): `scanpy`, `anndata`, `scipy`, `numpy`, `h5py`, `loompy` — for format conversion
- **r-converters** env (optional, WSL2): R, Seurat, zellkonverter — only needed for `.rds` files

Install via micromamba:
```bash
micromamba create -n scrna-python python=3.11 scanpy anndata scipy numpy h5py loompy pandas
```

## Development

```bash
# Install dev dependencies
pip install -e ".[dev]"

# Run tests
pytest tests/ -v

# Run specific test file
pytest tests/test_file_probe.py -v
```

## Architecture

```
GSE accession
  → GEO metadata probe (species, platform, conditions)
  → File probe (detect format, write per-file sidecars)
  → Route planner (pick conversion method by priority)
  → Execute conversion (with graceful degradation)
  → Gene symbol washing + species cross-validation
  → Validate + report
  → .h5ad output
```

Design: **Skills contain judgment; MCP tools execute deterministically.** See `CLAUDE.md` for full architecture.
