"""GSE Router MCP Server — stdio JSON-RPC entry point.

Exposes deterministic tools for probing, classifying, converting, and
validating single-cell RNA-seq datasets. The agent (skill) decides routes;
the MCP server executes operations.

Usage:
    python gse_router_mcp/server.py
"""

from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Any


# ---------------------------------------------------------------------------
# Tool registry — maps tool names to handler functions
# ---------------------------------------------------------------------------

# Add src to path so we can import gse_router
_src = Path(__file__).resolve().parent.parent / "src"
if str(_src) not in sys.path:
    sys.path.insert(0, str(_src))

# Load .env if present (for NCBI_API_KEY, NCBI_EMAIL, etc.)
try:
    from gse_router.config import _load_dotenv
    _load_dotenv()
except Exception:
    pass


def handle_probe_local_directory(params: dict) -> dict:
    """Probe a local directory and return file sidecars."""
    from gse_router.probes.file_probe import FileProbe

    path = Path(params["path"])
    probe = FileProbe()
    sidecars = probe.probe_directory(path)
    return {
        "sidecars": [sc.model_dump() for sc in sidecars],
        "count": len(sidecars),
    }


def handle_write_file_sidecars(params: dict) -> dict:
    """Write per-file JSON sidecars to disk."""
    from gse_router.probes.file_probe import FileProbe

    path = Path(params["path"])
    output_dir = Path(params.get("output_dir", path / "sidecars"))
    output_dir.mkdir(parents=True, exist_ok=True)

    probe = FileProbe()
    sidecars = probe.probe_directory(path)
    for sc in sidecars:
        sidecar_path = output_dir / f"{sc.file_id.replace(':', '_')}.json"
        sidecar_path.write_text(sc.model_dump_json(indent=2))

    return {"written": len(sidecars), "output_dir": str(output_dir)}


def handle_classify_scrna_files(params: dict) -> dict:
    """Classify probed files into scRNA-seq types."""
    from gse_router.router.classifier import FileClassifier
    from gse_router.schemas.sidecar import FileSidecar

    sidecars_raw = params.get("sidecars", [])
    sidecars: dict[str, FileSidecar] = {}
    for sc_data in sidecars_raw:
        sc = FileSidecar.model_validate(sc_data)
        sidecars[sc.file_id] = sc

    classifier = FileClassifier()
    groups = classifier.group_files(sidecars)
    strategy = classifier.select_strategy(groups)

    return {
        "groups": [g.model_dump() for g in groups],
        "strategy": strategy.value,
        "group_count": len(groups),
    }


def handle_plan_route(params: dict) -> dict:
    """Apply routing rules and produce a route plan."""
    from gse_router.router.route_planner import RoutePlanner
    from gse_router.schemas.manifest import DatasetManifest

    manifest_data = params["manifest"]
    manifest = DatasetManifest.model_validate(manifest_data)

    planner = RoutePlanner()
    route_plan = planner.plan(manifest)

    return json.loads(route_plan.model_dump_json())


def handle_validate_h5ad(params: dict) -> dict:
    """Validate an h5ad file against QC standards."""
    from gse_router.tools.validators import AnnDataValidator

    h5ad_path = Path(params["h5ad_path"])
    validator = AnnDataValidator()
    result = validator.validate(h5ad_path)

    return {
        "passed": result.passed,
        "metrics": {k: str(v) if not isinstance(v, (int, float, str, bool)) else v
                    for k, v in result.metrics.items()},
        "failures": result.failures,
        "warnings": result.warnings,
    }


def handle_build_manifest(params: dict) -> dict:
    """Build a dataset manifest from sidecars and optional study context."""
    from gse_router.router.classifier import FileClassifier
    from gse_router.schemas.manifest import DatasetManifest
    from gse_router.schemas.sidecar import FileSidecar

    sidecars_raw = params.get("sidecars", [])
    dataset_id = params["dataset_id"]
    title = params.get("title")
    organism = params.get("organism", [])

    sidecars: dict[str, FileSidecar] = {}
    for sc_data in sidecars_raw:
        sc = FileSidecar.model_validate(sc_data)
        sidecars[sc.file_id] = sc

    classifier = FileClassifier()
    file_groups = classifier.group_files(sidecars)
    strategy = classifier.select_strategy(file_groups)

    manifest = DatasetManifest(
        dataset_id=dataset_id,
        title=title,
        organism=organism,
        files=sidecars,
        file_groups=file_groups,
        selected_strategy=strategy,
    )

    return json.loads(manifest.model_dump_json())


# ---------------------------------------------------------------------------
# Conversion tool handlers
# ---------------------------------------------------------------------------


def _conversion_result_to_dict(result: Any) -> dict:
    """Serialize a ConversionResult to a JSON-safe dict."""
    return {
        "success": result.success,
        "output_path": str(result.output_path) if result.output_path else None,
        "provenance": result.provenance,
        "error": result.error,
        "warnings": result.warnings,
    }


def _metadata_result_to_dict(result: Any) -> dict:
    """Serialize a MetadataMergeResult to a JSON-safe dict."""
    return {
        "success": result.success,
        "output_path": str(result.output_path) if result.output_path else None,
        "n_cells_matched": result.n_cells_matched,
        "n_cells_total": result.n_cells_total,
        "n_columns_added": result.n_columns_added,
        "error": result.error,
        "provenance": result.provenance,
    }


def handle_convert_10x_mtx_to_h5ad(params: dict) -> dict:
    """Convert a 10x MEX triplet (matrix.mtx + barcodes.tsv + features.tsv) to h5ad."""
    from gse_router.tools.converters import FormatConverters

    mtx_path = Path(params["mtx_path"])
    barcodes_path = Path(params["barcodes_path"])
    features_path = Path(params["features_path"])
    sample_id = params.get("sample_id")
    var_names = params.get("var_names", "gene_symbols")
    output_dir = Path(params["output_dir"]) if "output_dir" in params else None

    converter = FormatConverters(output_dir=output_dir)
    result = converter.convert_10x_mtx_to_h5ad(
        mtx_path=mtx_path,
        barcodes_path=barcodes_path,
        features_path=features_path,
        sample_id=sample_id,
        var_names=var_names,
    )
    return _conversion_result_to_dict(result)


def handle_convert_10x_h5_to_h5ad(params: dict) -> dict:
    """Convert a 10x HDF5 file to h5ad."""
    from gse_router.tools.converters import FormatConverters

    h5_path = Path(params["h5_path"])
    sample_id = params.get("sample_id")
    genome = params.get("genome")
    output_dir = Path(params["output_dir"]) if "output_dir" in params else None

    converter = FormatConverters(output_dir=output_dir)
    result = converter.convert_10x_h5_to_h5ad(
        h5_path=h5_path,
        sample_id=sample_id,
        genome=genome,
    )
    return _conversion_result_to_dict(result)


def handle_convert_loom_to_h5ad(params: dict) -> dict:
    """Convert a Loom file to h5ad."""
    from gse_router.tools.converters import FormatConverters

    loom_path = Path(params["loom_path"])
    sample_id = params.get("sample_id")
    output_dir = Path(params["output_dir"]) if "output_dir" in params else None

    converter = FormatConverters(output_dir=output_dir)
    result = converter.convert_loom_to_h5ad(
        loom_path=loom_path,
        sample_id=sample_id,
    )
    return _conversion_result_to_dict(result)


def handle_read_h5ad_validate(params: dict) -> dict:
    """Validate and standardize an existing h5ad file."""
    from gse_router.tools.converters import FormatConverters

    h5ad_path = Path(params["h5ad_path"])
    sample_id = params.get("sample_id")
    output_dir = Path(params["output_dir"]) if "output_dir" in params else None

    converter = FormatConverters(output_dir=output_dir)
    result = converter.validate_and_standardize_h5ad(
        h5ad_path=h5ad_path,
        sample_id=sample_id,
    )
    return _conversion_result_to_dict(result)


def handle_convert_raw_umi_matrix_to_h5ad(params: dict) -> dict:
    """Convert a raw UMI count matrix (TXT/CSV/RDS) to h5ad."""
    from gse_router.tools.converters import FormatConverters

    matrix_path = Path(params["matrix_path"])
    sample_id = params.get("sample_id")
    delimiter = params.get("delimiter")
    file_format = params.get("file_format", "txt")
    output_dir = Path(params["output_dir"]) if "output_dir" in params else None

    converter = FormatConverters(output_dir=output_dir)
    result = converter.convert_raw_umi_matrix_to_h5ad(
        matrix_path=matrix_path,
        sample_id=sample_id,
        delimiter=delimiter,
        file_format=file_format,
    )
    return _conversion_result_to_dict(result)


def handle_convert_normalized_matrix_to_h5ad(params: dict) -> dict:
    """Convert a normalized expression matrix to h5ad (limited use)."""
    from gse_router.tools.converters import FormatConverters

    matrix_path = Path(params["matrix_path"])
    sample_id = params.get("sample_id")
    delimiter = params.get("delimiter")
    normalization_method = params.get("normalization_method", "unknown")
    output_dir = Path(params["output_dir"]) if "output_dir" in params else None

    converter = FormatConverters(output_dir=output_dir)
    result = converter.convert_normalized_matrix_to_h5ad(
        matrix_path=matrix_path,
        sample_id=sample_id,
        delimiter=delimiter,
        normalization_method=normalization_method,
    )
    return _conversion_result_to_dict(result)


def handle_convert_csv_matrix_to_h5ad(params: dict) -> dict:
    """Convert a CSV/TSV count matrix to h5ad."""
    from gse_router.tools.converters import FormatConverters

    csv_path = Path(params["csv_path"])
    sample_id = params.get("sample_id")
    delimiter = params.get("delimiter")
    transpose = params.get("transpose", False)
    is_raw_counts = params.get("is_raw_counts")
    output_dir = Path(params["output_dir"]) if "output_dir" in params else None

    converter = FormatConverters(output_dir=output_dir)
    result = converter.convert_csv_matrix_to_h5ad(
        csv_path=csv_path,
        sample_id=sample_id,
        delimiter=delimiter,
        transpose=transpose,
        is_raw_counts=is_raw_counts,
    )
    return _conversion_result_to_dict(result)


def handle_merge_cell_metadata(params: dict) -> dict:
    """Merge cell-level metadata from a CSV/TSV file into adata.obs."""
    from gse_router.tools.metadata_tools import MetadataTools

    h5ad_path = Path(params["h5ad_path"])
    metadata_path = Path(params["metadata_path"])
    join_column = params.get("join_column")  # Agent provides explicitly after inspecting metadata
    delimiter = params.get("delimiter")
    output_path = Path(params["output_path"]) if "output_path" in params else None

    tools = MetadataTools()
    result = tools.merge_cell_metadata(
        h5ad_path=h5ad_path,
        metadata_path=metadata_path,
        join_column=join_column,
        delimiter=delimiter,
        output_path=output_path,
    )
    return _metadata_result_to_dict(result)


def handle_merge_gene_metadata(params: dict) -> dict:
    """Merge gene-level annotations into adata.var."""
    from gse_router.tools.metadata_tools import MetadataTools

    h5ad_path = Path(params["h5ad_path"])
    annotation_path = Path(params["annotation_path"])
    join_column = params.get("join_column", "gene_id")
    output_path = Path(params["output_path"]) if "output_path" in params else None

    tools = MetadataTools()
    result = tools.merge_gene_annotation(
        h5ad_path=h5ad_path,
        annotation_path=annotation_path,
        join_column=join_column,
        output_path=output_path,
    )
    return _metadata_result_to_dict(result)


# ---------------------------------------------------------------------------
# RDS conversion tool handlers
# ---------------------------------------------------------------------------


def handle_convert_rds_seurat_to_h5ad(params: dict) -> dict:
    """Convert a Seurat RDS file to h5ad."""
    from gse_router.tools.converters import FormatConverters

    rds_path = Path(params["rds_path"])
    sample_id = params.get("sample_id")
    assay = params.get("assay", "RNA")
    output_dir = Path(params["output_dir"]) if "output_dir" in params else None

    converter = FormatConverters(output_dir=output_dir)
    result = converter.convert_rds_to_h5ad(
        rds_path,
        sample_id=sample_id,
        object_type="seurat",
    )
    return _conversion_result_to_dict(result)


def handle_convert_rds_sce_to_h5ad(params: dict) -> dict:
    """Convert a SingleCellExperiment RDS file to h5ad."""
    from gse_router.tools.converters import FormatConverters

    rds_path = Path(params["rds_path"])
    sample_id = params.get("sample_id")
    output_dir = Path(params["output_dir"]) if "output_dir" in params else None

    converter = FormatConverters(output_dir=output_dir)
    result = converter.convert_rds_to_h5ad(
        rds_path,
        sample_id=sample_id,
        object_type="sce",
    )
    return _conversion_result_to_dict(result)


def handle_convert_h5seurat_to_h5ad(params: dict) -> dict:
    """Convert an H5Seurat file to h5ad."""
    from gse_router.tools.converters import FormatConverters

    h5seurat_path = Path(params["h5seurat_path"])
    sample_id = params.get("sample_id")
    output_dir = Path(params["output_dir"]) if "output_dir" in params else None

    converter = FormatConverters(output_dir=output_dir)
    result = converter.convert_rds_to_h5ad(
        h5seurat_path,
        sample_id=sample_id,
        object_type="h5seurat",
    )
    return _conversion_result_to_dict(result)


# ---------------------------------------------------------------------------
# FASTQ/SRA pipeline handlers — disabled by default (cardinal rule: skip
# FASTQ if any processed data exists; priority 8 of 8 routing levels)
# ---------------------------------------------------------------------------


def _fastq_disabled_response(tool_name: str) -> dict:
    """Return a consistent disabled-by-default response for FASTQ/SRA tools."""
    return {
        "success": False,
        "error": (
            f"Raw FASTQ/SRA processing ('{tool_name}') is disabled by default. "
            "The project treats FASTQ/SRA as the lowest priority route (level 8 of 8). "
            "Per the cardinal rule: 'Skip FASTQ if processed data exists.' If any usable "
            "processed count/object file is available (priority 1–7), do NOT process raw "
            "reads. Only fall back to FASTQ/SRA if ALL higher-priority routes fail."
        ),
        "warnings": [
            "FASTQ/SRA processing is resource-intensive and should be a last resort.",
            "Many GEO 'RAW.tar' files contain processed matrices, not FASTQ — always check file contents first.",
        ],
        "note": (
            "To enable this tool, all of the following are required:\n"
            "  * Explicit opt-in via skip_raw_fastq=False in the route config\n"
            "  * WSL2 Ubuntu (on Windows) — raw read pipelines require Linux\n"
            "  * raw-pipeline micromamba environment (nextflow, nf-core/fetchngs,\n"
            "    nf-core/scrnaseq, STARsolo)\n"
            "  * Pre-built reference genome index for the target species\n"
            "  * Sufficient disk space for FASTQ download and alignment\n"
            "These requirements are intentional: the project is optimized for\n"
            "processed-data reuse, not raw-data reprocessing."
        ),
    }


def handle_fetchngs_download(params: dict) -> dict:
    """Download FASTQ files via nf-core/fetchngs (disabled by default)."""
    return _fastq_disabled_response("fetchngs_download")


def handle_fasterq_dump(params: dict) -> dict:
    """Convert SRA to FASTQ via fasterq-dump (disabled by default)."""
    return _fastq_disabled_response("fasterq_dump")


def handle_cellranger_count(params: dict) -> dict:
    """Run Cell Ranger count on FASTQ data (disabled by default)."""
    return _fastq_disabled_response("cellranger_count")


def handle_starsolo(params: dict) -> dict:
    """Run STARsolo alignment on FASTQ data (disabled by default)."""
    return _fastq_disabled_response("starsolo")


def handle_scrnaseq_pipeline(params: dict) -> dict:
    """Run the full nf-core/scrnaseq pipeline (disabled by default)."""
    return _fastq_disabled_response("scrnaseq_pipeline")


# ---------------------------------------------------------------------------
# GEO metadata tool handlers
# ---------------------------------------------------------------------------


def handle_fetch_geo_metadata(params: dict) -> dict:
    """Query GEO metadata for a GSE accession: organism, platform, samples, conditions."""
    import asyncio

    from gse_router_mcp.tools.geo_metadata import fetch_geo_metadata

    gse_id = params["gse_id"]
    email = params.get("email")
    return asyncio.run(fetch_geo_metadata(gse_id, email=email))


def handle_list_geo_supplementary(params: dict) -> dict:
    """List supplementary files for a GSE accession from GEO FTP."""
    import asyncio

    from gse_router_mcp.tools.geo_metadata import list_geo_supplementary

    gse_id = params["gse_id"]
    return asyncio.run(list_geo_supplementary(gse_id))


def handle_read_geo_filelist(params: dict) -> dict:
    """Read and parse GEO supplementary filelist.txt with checksums."""
    import asyncio

    from gse_router_mcp.tools.geo_metadata import read_geo_filelist

    gse_id = params["gse_id"]
    return asyncio.run(read_geo_filelist(gse_id))


def handle_download_geo_supplementary(params: dict) -> dict:
    """Download GEO supplementary files for a GSE dataset."""
    import asyncio

    from gse_router_mcp.tools.geo_metadata import download_geo_supplementary

    gse_id = params["gse_id"]
    data_dir = params.get("data_dir", f"data/{gse_id}/supplementary")
    file_types = params.get("file_types")
    max_size_mb = params.get("max_size_mb", 500)
    download_method = params.get("download_method", "auto")
    return asyncio.run(
        download_geo_supplementary(
            gse_id=gse_id,
            data_dir=data_dir,
            file_types=file_types,
            max_size_mb=max_size_mb,
            download_method=download_method,
        )
    )


def handle_check_mirrors(params: dict) -> dict:
    """Check reachability of GEO/SRA mirror sites."""
    import asyncio

    from gse_router.tools.mirrors import MirrorResolver

    gse_id = params.get("gse_id")
    resolver = MirrorResolver()
    health = asyncio.run(resolver.check_mirror_health(gse_id=gse_id))
    return {"mirror_health": health}


# Tool dispatch table
TOOLS: dict[str, dict] = {
    "probe_local_directory": {
        "handler": handle_probe_local_directory,
        "description": "Recursively probe a directory, detect file formats via magic bytes and headers",
        "inputSchema": {
            "type": "object",
            "properties": {
                "path": {"type": "string", "description": "Path to the dataset directory"},
            },
            "required": ["path"],
        },
    },
    "write_file_sidecars": {
        "handler": handle_write_file_sidecars,
        "description": "Write per-file JSON sidecars for all files in a directory",
        "inputSchema": {
            "type": "object",
            "properties": {
                "path": {"type": "string", "description": "Path to the dataset directory"},
                "output_dir": {"type": "string", "description": "Directory for sidecar output"},
            },
            "required": ["path"],
        },
    },
    "classify_scrna_files": {
        "handler": handle_classify_scrna_files,
        "description": "Classify probed files into scRNA-seq types and group related files",
        "inputSchema": {
            "type": "object",
            "properties": {
                "sidecars": {
                    "type": "array",
                    "description": "List of FileSidecar JSON objects",
                    "items": {"type": "object"},
                },
            },
            "required": ["sidecars"],
        },
    },
    "build_dataset_manifest": {
        "handler": handle_build_manifest,
        "description": "Combine sidecars + study context into a DatasetManifest JSON",
        "inputSchema": {
            "type": "object",
            "properties": {
                "dataset_id": {"type": "string", "description": "GEO accession (e.g., GSE123456)"},
                "sidecars": {"type": "array", "items": {"type": "object"}},
                "title": {"type": "string"},
                "organism": {"type": "array", "items": {"type": "string"}},
            },
            "required": ["dataset_id", "sidecars"],
        },
    },
    "plan_route": {
        "handler": handle_plan_route,
        "description": "Apply routing rules to manifest, produce an ordered RoutePlan",
        "inputSchema": {
            "type": "object",
            "properties": {
                "manifest": {"type": "object", "description": "DatasetManifest JSON"},
            },
            "required": ["manifest"],
        },
    },
    "validate_h5ad": {
        "handler": handle_validate_h5ad,
        "description": "Run QC checks on h5ad: cell/gene counts, MT%, raw counts verification",
        "inputSchema": {
            "type": "object",
            "properties": {
                "h5ad_path": {"type": "string", "description": "Path to the .h5ad file"},
            },
            "required": ["h5ad_path"],
        },
    },
    "convert_10x_mtx_to_h5ad": {
        "handler": handle_convert_10x_mtx_to_h5ad,
        "description": "Convert a 10x MEX triplet (matrix.mtx, barcodes.tsv, features.tsv) to a standardized .h5ad file",
        "inputSchema": {
            "type": "object",
            "properties": {
                "mtx_path": {"type": "string", "description": "Path to the MatrixMarket .mtx file"},
                "barcodes_path": {"type": "string", "description": "Path to barcodes.tsv(.gz)"},
                "features_path": {"type": "string", "description": "Path to features/genes.tsv(.gz)"},
                "sample_id": {"type": "string", "description": "Optional sample identifier"},
                "var_names": {"type": "string", "description": "How to name variables: 'gene_symbols' or 'gene_ids'"},
                "output_dir": {"type": "string", "description": "Output directory for the .h5ad file"},
            },
            "required": ["mtx_path", "barcodes_path", "features_path"],
        },
    },
    "convert_10x_h5_to_h5ad": {
        "handler": handle_convert_10x_h5_to_h5ad,
        "description": "Convert a 10x HDF5 (.h5) file to a standardized .h5ad file",
        "inputSchema": {
            "type": "object",
            "properties": {
                "h5_path": {"type": "string", "description": "Path to the 10x .h5 file"},
                "sample_id": {"type": "string", "description": "Optional sample identifier"},
                "genome": {"type": "string", "description": "Optional genome string for multi-genome H5 files"},
                "output_dir": {"type": "string", "description": "Output directory for the .h5ad file"},
            },
            "required": ["h5_path"],
        },
    },
    "convert_loom_to_h5ad": {
        "handler": handle_convert_loom_to_h5ad,
        "description": "Convert a Loom (.loom) file to a standardized .h5ad file",
        "inputSchema": {
            "type": "object",
            "properties": {
                "loom_path": {"type": "string", "description": "Path to the .loom file"},
                "sample_id": {"type": "string", "description": "Optional sample identifier"},
                "output_dir": {"type": "string", "description": "Output directory for the .h5ad file"},
            },
            "required": ["loom_path"],
        },
    },
    "read_h5ad_validate": {
        "handler": handle_read_h5ad_validate,
        "description": "Validate an existing .h5ad and standardize its schema (preserve raw counts, normalize metadata columns)",
        "inputSchema": {
            "type": "object",
            "properties": {
                "h5ad_path": {"type": "string", "description": "Path to the .h5ad file"},
                "sample_id": {"type": "string", "description": "Optional sample identifier to add to obs"},
                "output_dir": {"type": "string", "description": "Output directory for the standardized .h5ad file"},
            },
            "required": ["h5ad_path"],
        },
    },
    "convert_raw_umi_matrix_to_h5ad": {
        "handler": handle_convert_raw_umi_matrix_to_h5ad,
        "description": "Convert a raw UMI count matrix (TXT/CSV/RDS) to .h5ad — priority 5 route: confirmed raw counts",
        "inputSchema": {
            "type": "object",
            "properties": {
                "matrix_path": {"type": "string", "description": "Path to the count matrix file"},
                "sample_id": {"type": "string", "description": "Optional sample identifier"},
                "delimiter": {"type": "string", "description": "Column delimiter (auto-detected if omitted)"},
                "file_format": {"type": "string", "description": "File format: 'txt', 'csv', 'tsv', or 'rds'"},
                "output_dir": {"type": "string", "description": "Output directory for the .h5ad file"},
            },
            "required": ["matrix_path"],
        },
    },
    "convert_normalized_matrix_to_h5ad": {
        "handler": handle_convert_normalized_matrix_to_h5ad,
        "description": "Convert a normalized expression matrix to .h5ad — LIMITED USE: raw counts NOT preserved",
        "inputSchema": {
            "type": "object",
            "properties": {
                "matrix_path": {"type": "string", "description": "Path to the normalized matrix file"},
                "sample_id": {"type": "string", "description": "Optional sample identifier"},
                "delimiter": {"type": "string", "description": "Column delimiter (auto-detected if omitted)"},
                "normalization_method": {"type": "string", "description": "e.g. 'log2TPM', 'FPKM', 'log-normalized'"},
                "output_dir": {"type": "string", "description": "Output directory for the .h5ad file"},
            },
            "required": ["matrix_path"],
        },
    },
    "convert_csv_matrix_to_h5ad": {
        "handler": handle_convert_csv_matrix_to_h5ad,
        "description": "Convert a generic CSV/TSV count matrix to .h5ad — auto-detects gene vs cell orientation",
        "inputSchema": {
            "type": "object",
            "properties": {
                "csv_path": {"type": "string", "description": "Path to the CSV/TSV file"},
                "sample_id": {"type": "string", "description": "Optional sample identifier"},
                "delimiter": {"type": "string", "description": "Column delimiter (auto-detected if omitted)"},
                "transpose": {"type": "boolean", "description": "If true, transpose the matrix before building AnnData"},
                "is_raw_counts": {"type": "boolean", "description": "Whether these are known raw counts (None = unknown)"},
                "output_dir": {"type": "string", "description": "Output directory for the .h5ad file"},
            },
            "required": ["csv_path"],
        },
    },
    "merge_cell_metadata": {
        "handler": handle_merge_cell_metadata,
        "description": "Merge cell-level metadata from a CSV/TSV file into adata.obs (joins on barcode by default)",
        "inputSchema": {
            "type": "object",
            "properties": {
                "h5ad_path": {"type": "string", "description": "Path to the .h5ad file to enrich"},
                "metadata_path": {"type": "string", "description": "Path to the CSV/TSV metadata file"},
                "join_column": {"type": "string", "description": "Column in metadata to join on (default: 'barcode')"},
                "delimiter": {"type": "string", "description": "Column delimiter in metadata file (auto-detected if omitted)"},
                "output_path": {"type": "string", "description": "Path for the enriched .h5ad (overwrites input if omitted)"},
            },
            "required": ["h5ad_path", "metadata_path"],
        },
    },
    "merge_gene_metadata": {
        "handler": handle_merge_gene_metadata,
        "description": "Merge gene-level annotations (CSV/TSV/GTF) into adata.var",
        "inputSchema": {
            "type": "object",
            "properties": {
                "h5ad_path": {"type": "string", "description": "Path to the .h5ad file to enrich"},
                "annotation_path": {"type": "string", "description": "Path to the gene annotation file (CSV/TSV/GTF)"},
                "join_column": {"type": "string", "description": "Column in annotation to join on (default: 'gene_id')"},
                "output_path": {"type": "string", "description": "Path for the enriched .h5ad (overwrites input if omitted)"},
            },
            "required": ["h5ad_path", "annotation_path"],
        },
    },
    "convert_rds_seurat_to_h5ad": {
        "handler": handle_convert_rds_seurat_to_h5ad,
        "description": "Convert a Seurat RDS file to a standardized .h5ad file — priority 6 route: R bridge conversion",
        "inputSchema": {
            "type": "object",
            "properties": {
                "rds_path": {"type": "string", "description": "Path to the Seurat .rds file"},
                "sample_id": {"type": "string", "description": "Optional sample identifier to add to obs"},
                "assay": {"type": "string", "description": "Assay name to extract (default: 'RNA')"},
                "output_dir": {"type": "string", "description": "Output directory for the .h5ad file"},
            },
            "required": ["rds_path"],
        },
    },
    "convert_rds_sce_to_h5ad": {
        "handler": handle_convert_rds_sce_to_h5ad,
        "description": "Convert a SingleCellExperiment RDS file to a standardized .h5ad file — priority 6 route: R bridge conversion",
        "inputSchema": {
            "type": "object",
            "properties": {
                "rds_path": {"type": "string", "description": "Path to the SingleCellExperiment .rds file"},
                "sample_id": {"type": "string", "description": "Optional sample identifier to add to obs"},
                "output_dir": {"type": "string", "description": "Output directory for the .h5ad file"},
            },
            "required": ["rds_path"],
        },
    },
    "convert_h5seurat_to_h5ad": {
        "handler": handle_convert_h5seurat_to_h5ad,
        "description": "Convert an H5Seurat file to a standardized .h5ad file — priority 6 route: R bridge conversion",
        "inputSchema": {
            "type": "object",
            "properties": {
                "h5seurat_path": {"type": "string", "description": "Path to the .h5seurat file"},
                "sample_id": {"type": "string", "description": "Optional sample identifier to add to obs"},
                "output_dir": {"type": "string", "description": "Output directory for the .h5ad file"},
            },
            "required": ["h5seurat_path"],
        },
    },
    "fetchngs_download": {
        "handler": handle_fetchngs_download,
        "description": "Download FASTQ from SRA via nf-core/fetchngs (priority 8 — DISABLED by default; use only if no processed data exists)",
        "inputSchema": {
            "type": "object",
            "properties": {
                "gse_id": {"type": "string", "description": "GEO accession ID (e.g., GSE123456)"},
                "output_dir": {"type": "string", "description": "Output directory for downloaded FASTQ files"},
                "download_method": {"type": "string", "description": "Download method: 'sratools', 'ftp', or 'ena'"},
            },
            "required": ["gse_id"],
        },
    },
    "fasterq_dump": {
        "handler": handle_fasterq_dump,
        "description": "Convert SRA accessions to FASTQ via fasterq-dump (priority 8 — DISABLED by default; use only if no processed data exists)",
        "inputSchema": {
            "type": "object",
            "properties": {
                "sra_accession": {"type": "string", "description": "SRA run accession (e.g., SRR12345678)"},
                "output_dir": {"type": "string", "description": "Output directory for FASTQ files"},
                "threads": {"type": "integer", "description": "Number of threads for parallel processing"},
            },
            "required": ["sra_accession"],
        },
    },
    "cellranger_count": {
        "handler": handle_cellranger_count,
        "description": "Run Cell Ranger count on FASTQ data (priority 8 — DISABLED by default; use only if no processed data exists)",
        "inputSchema": {
            "type": "object",
            "properties": {
                "fastq_dir": {"type": "string", "description": "Directory containing FASTQ files"},
                "sample_id": {"type": "string", "description": "Sample identifier"},
                "transcriptome": {"type": "string", "description": "Path to Cell Ranger reference transcriptome"},
                "output_dir": {"type": "string", "description": "Output directory for Cell Ranger results"},
            },
            "required": ["fastq_dir", "sample_id", "transcriptome"],
        },
    },
    "starsolo": {
        "handler": handle_starsolo,
        "description": "Run STARsolo alignment on FASTQ data (priority 8 — DISABLED by default; use only if no processed data exists)",
        "inputSchema": {
            "type": "object",
            "properties": {
                "fastq_path": {"type": "string", "description": "Path to FASTQ file(s) — can be a directory or glob pattern"},
                "sample_id": {"type": "string", "description": "Sample identifier"},
                "genome_dir": {"type": "string", "description": "Path to STAR genome index directory"},
                "output_dir": {"type": "string", "description": "Output directory for STARsolo results"},
                "paired": {"type": "boolean", "description": "Whether reads are paired-end (default: auto-detect)"},
            },
            "required": ["fastq_path", "sample_id", "genome_dir"],
        },
    },
    "scrnaseq_pipeline": {
        "handler": handle_scrnaseq_pipeline,
        "description": "Run the full nf-core/scrnaseq pipeline (priority 8 — DISABLED by default; use only if no processed data exists)",
        "inputSchema": {
            "type": "object",
            "properties": {
                "gse_id": {"type": "string", "description": "GEO accession ID to process (e.g., GSE123456)"},
                "output_dir": {"type": "string", "description": "Output directory for pipeline results"},
                "profile": {"type": "string", "description": "Nextflow execution profile (e.g., 'docker', 'singularity', 'conda')"},
            },
            "required": ["gse_id"],
        },
    },
    "fetch_geo_metadata": {
        "handler": handle_fetch_geo_metadata,
        "description": "Query GEO metadata for a GSE accession: organism, platform, samples, conditions",
        "inputSchema": {
            "type": "object",
            "properties": {
                "gse_id": {"type": "string", "description": "GEO accession (e.g., GSE123456)"},
                "email": {"type": "string", "description": "Optional NCBI contact email for E-utilities rate-limiting"},
            },
            "required": ["gse_id"],
        },
    },
    "list_geo_supplementary": {
        "handler": handle_list_geo_supplementary,
        "description": "List supplementary files for a GSE accession from GEO FTP",
        "inputSchema": {
            "type": "object",
            "properties": {
                "gse_id": {"type": "string", "description": "GEO accession (e.g., GSE123456)"},
            },
            "required": ["gse_id"],
        },
    },
    "read_geo_filelist": {
        "handler": handle_read_geo_filelist,
        "description": "Read and parse GEO supplementary filelist.txt with checksums",
        "inputSchema": {
            "type": "object",
            "properties": {
                "gse_id": {"type": "string", "description": "GEO accession (e.g., GSE123456)"},
            },
            "required": ["gse_id"],
        },
    },
    "download_geo_supplementary": {
        "handler": handle_download_geo_supplementary,
        "description": "Download GEO supplementary files for a GSE dataset",
        "inputSchema": {
            "type": "object",
            "properties": {
                "gse_id": {"type": "string", "description": "GEO accession (e.g., GSE123456)"},
                "data_dir": {"type": "string", "description": "Local directory for downloaded files (default: data/<gse_id>/supplementary)"},
                "file_types": {
                    "type": "array",
                    "items": {"type": "string"},
                    "description": "Optional list of file extensions to filter (e.g., ['h5ad', 'loom', 'mtx', 'gz'])",
                },
                "max_size_mb": {"type": "number", "description": "Skip files larger than this in megabytes (default: 500)"},
                "download_method": {
                    "type": "string",
                    "description": "Download backend: 'auto' (smart), 'aspera', 'aria2', or 'httpx'. Default: 'auto'.",
                    "enum": ["auto", "aspera", "aria2", "httpx"],
                },
            },
            "required": ["gse_id"],
        },
    },
    "check_mirrors": {
        "handler": handle_check_mirrors,
        "description": "Check reachability of GEO/SRA mirror sites (NCBI, EBI, CNCB, DDBJ). Use before downloading to verify mirror availability.",
        "inputSchema": {
            "type": "object",
            "properties": {
                "gse_id": {"type": "string", "description": "Optional GEO accession for targeted health check on specific dataset directory"},
            },
        },
    },
}


# ---------------------------------------------------------------------------
# JSON-RPC stdio server
# ---------------------------------------------------------------------------

def _send_response(response_id: int | str | None, result: Any) -> None:
    """Send a JSON-RPC success response."""
    msg = json.dumps({"jsonrpc": "2.0", "id": response_id, "result": result})
    sys.stdout.write(msg + "\n")
    sys.stdout.flush()


def _send_error(response_id: int | str | None, code: int, message: str) -> None:
    """Send a JSON-RPC error response."""
    msg = json.dumps({
        "jsonrpc": "2.0",
        "id": response_id,
        "error": {"code": code, "message": message},
    })
    sys.stdout.write(msg + "\n")
    sys.stdout.flush()


def _handle_request(request: dict) -> None:
    """Handle a single JSON-RPC request."""
    req_id = request.get("id")
    method = request.get("method", "")
    params = request.get("params", {})

    if method == "initialize":
        _send_response(req_id, {
            "protocolVersion": "2024-11-05",
            "serverInfo": {"name": "gse-router", "version": "0.1.0"},
            "capabilities": {"tools": {}},
        })
        return

    if method == "notifications/initialized":
        # No response needed for notifications
        return

    if method == "tools/list":
        tools_list = [
            {
                "name": name,
                "description": info["description"],
                "inputSchema": info["inputSchema"],
            }
            for name, info in TOOLS.items()
        ]
        _send_response(req_id, {"tools": tools_list})
        return

    if method == "tools/call":
        tool_name = params.get("name", "")
        tool_args = params.get("arguments", {})
        tool_info = TOOLS.get(tool_name)

        if tool_info is None:
            _send_error(req_id, -32601, f"Unknown tool: {tool_name}")
            return

        try:
            result = tool_info["handler"](tool_args)
            _send_response(req_id, {
                "content": [{"type": "text", "text": json.dumps(result, indent=2, default=str)}],
            })
        except Exception as exc:
            _send_error(req_id, -32000, f"Tool execution failed: {exc}")
        return

    _send_error(req_id, -32601, f"Unknown method: {method}")


def main() -> None:
    """Run the MCP server on stdio."""
    for line in sys.stdin:
        line = line.strip()
        if not line:
            continue
        try:
            request = json.loads(line)
            _handle_request(request)
        except json.JSONDecodeError:
            _send_error(None, -32700, "Parse error")
        except Exception as exc:
            _send_error(None, -32603, f"Internal error: {exc}")


if __name__ == "__main__":
    main()
