"""Per-file JSON sidecar schema.

Each downloaded file gets one sidecar describing its detected type,
confidence, related files, and recommended processing routes.
"""

from __future__ import annotations

from datetime import datetime, timezone
from enum import Enum

from pydantic import BaseModel, Field


class FileDetectedType(str, Enum):
    """Detected single-cell file format.

    Types are ordered approximately by routing priority (highest first)
    to guide the classifier's preference when multiple formats are found.
    """

    # High-priority processed objects
    H5AD = "h5ad"
    H5_10X = "10x_h5"
    LOOM = "loom"
    # Count matrices
    MEX_MATRIX_10X = "10x_mtx_matrix"
    MEX_BARCODES_10X = "10x_mtx_barcodes"
    MEX_FEATURES_10X = "10x_mtx_features"
    # Raw UMI matrices (explicitly raw counts in delimited format)
    RAW_UMI_MATRIX = "raw_umi_matrix"
    # R objects
    RDS_SEURAT = "seurat_rds"
    RDS_SCE = "sce_rds"
    H5SEURAT = "h5seurat"
    # Generic delimited matrices (ambiguous: could be counts or normalized)
    CSV_MATRIX = "csv_matrix"
    TSV_MATRIX = "tsv_matrix"
    TXT_MATRIX = "txt_matrix"
    # Normalized matrix — do NOT treat as raw counts
    NORMALIZED_MATRIX = "normalized_matrix"
    # Raw sequencing
    FASTQ = "fastq"
    SRA_ARCHIVE = "sra_archive"
    # Alignment
    BAM = "bam"
    CRAM = "cram"
    # Metadata and annotation
    CELL_METADATA = "cell_metadata"
    SAMPLE_METADATA = "sample_metadata"
    METADATA = "metadata"  # generic metadata (e.g., GEO manifest, checksum files)
    GENE_ANNOTATION = "gene_annotation"
    # Non-count supplementary data (must NOT be routed as count matrices)
    IMMUNE_RECEPTOR_DATA = "immune_receptor_data"  # BCR/TCR contig annotations
    PEAK_DATA = "peak_data"  # BED peaks (ATAC/ChIP-seq)
    DIFF_RESULTS = "diff_results"  # DiffBind/DESeq2 differential expression results
    # Specialized matrix formats
    MATRIX_MAT = "matlab_matrix"  # MATLAB .mat matrix file (SPLiT-seq, etc.)
    # Special
    FILELIST_TXT = "filelist_txt"  # GEO supplementary file listing
    ARCHIVE = "archive"
    RAW_TAR = "raw_tar"  # GEO RAW.tar — may contain processed data, not necessarily FASTQ
    UNKNOWN = "unknown"


class ProcessingStage(str, Enum):
    """Processing stage of the detected file."""

    RAW_SEQUENCING = "raw_sequencing"
    RAW_COUNTS = "processed_counts"
    NORMALIZED = "normalized"
    PROCESSED_OBJECT = "processed_object"
    METADATA = "metadata"
    ANNOTATION = "annotation"
    UNKNOWN = "unknown"


class RouteRecommendation(BaseModel):
    """A recommended processing route for this file."""

    route: str = Field(..., description="Route identifier, e.g. 'read_10x_mtx_to_h5ad'")
    priority: int = Field(..., description="Priority (1 = highest)")
    description: str | None = Field(None, description="Human-readable description of this route")


class FileSidecar(BaseModel):
    """Complete per-file sidecar describing one downloaded file."""

    file_id: str = Field(
        ...,
        description="Content-hash based ID, e.g. 'sha256:abcdef...'",
        pattern=r"^sha256:[a-f0-9]{64}$",
    )
    path: str = Field(..., description="Relative path to the file within the dataset directory")
    filename: str = Field(..., description="Original filename")
    size_bytes: int = Field(..., ge=0, description="File size in bytes")
    extension: str = Field(default="", description="File extension, lowercase")
    mime_type: str | None = Field(None, description="Detected MIME type")

    detected_type: FileDetectedType = Field(..., description="Detected single-cell file format")
    processing_stage: ProcessingStage = Field(..., description="Inferred processing stage")
    modality: str = Field(default="scRNA-seq", description="Modality, e.g. 'scRNA-seq', 'scATAC-seq'")
    organism: str | None = Field(None, description="Detected or inferred organism")

    sample_id_guess: str | None = Field(None, description="Best guess at sample identifier")
    gse: str | None = Field(None, description="Parent GSE accession")
    gsm: str | None = Field(None, description="Parent GSM accession")
    srr: str | None = Field(None, description="Parent SRR accession")

    paired_files: list[str] = Field(
        default_factory=list,
        description="File IDs of companion files (e.g., barcodes + features for 10x matrix)",
    )
    confidence: float = Field(..., ge=0.0, le=1.0, description="Detection confidence score")
    evidence: list[str] = Field(
        default_factory=list,
        description="Specific observations supporting the classification",
    )

    # Data quality flags
    is_raw_counts: bool | None = Field(
        None,
        description="Whether this file contains raw (un-normalized) counts. "
                    "None = unknown, True = confirmed raw, False = confirmed normalized. "
                    "Critical for routing: never treat normalized data as raw counts.",
    )
    limited_use: bool = Field(
        default=False,
        description="If True, data can be converted but should be marked as limited-use "
                    "(e.g., normalized-only matrices without raw counts available).",
    )

    recommended_routes: list[RouteRecommendation] = Field(
        default_factory=list,
        description="Ordered list of recommended processing routes",
    )

    parent_file_ids: list[str] = Field(
        default_factory=list,
        description="File IDs this file was extracted from (e.g., archive)",
    )
    extracted_files: list[str] = Field(
        default_factory=list,
        description="File IDs of files extracted from this archive",
    )

    probe_version: str = Field(default="v1", description="Version of the probe tool")
    probed_at: datetime = Field(default_factory=lambda: datetime.now(timezone.utc), description="When this file was probed")
    checksum_sha256: str | None = Field(None, description="SHA-256 checksum of the file")
