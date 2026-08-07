"""Dataset manifest schema.

Aggregates all per-file sidecars into a structured dataset-level view
with file groups, sample assignments, and combined evidence.
"""

from __future__ import annotations

from datetime import datetime, timezone
from enum import Enum

from pydantic import BaseModel, Field

from .sidecar import FileSidecar


class FileGroupType(str, Enum):
    """Type of a related file group.

    Ordered by routing priority — the classifier prefers higher-priority
    types when deciding which data to process first.
    """

    H5AD = "h5ad"                    # Priority 1: existing h5ad
    H5_10X = "10x_h5"                # Priority 2: 10x HDF5
    LOOM = "loom"                    # Priority 3: Loom file
    MEX_10X = "10x_mex"              # Priority 4: 10x MEX triplet
    RAW_UMI_MATRIX = "raw_umi_matrix"  # Priority 5: raw UMI count matrix (TXT/CSV/RDS)
    SEURAT_RDS = "seurat_rds"        # Priority 6: Seurat RDS
    SCE_RDS = "sce_rds"              # Priority 6: SCE RDS
    H5SEURAT = "h5seurat"            # Priority 6: h5Seurat
    CSV_MATRIX = "csv_matrix"        # Generic delimited matrix (ambiguous)
    NORMALIZED_MATRIX = "normalized_matrix"  # Priority 7: normalized only — limited use
    FASTQ_PAIRED = "fastq_paired"    # Priority 8: raw FASTQ (avoid if processed exists)
    FASTQ_SINGLE = "fastq_single"    # Priority 8: raw FASTQ single-end
    ARCHIVE_10X = "archive_10x"     # Archive known to contain 10x MEX data (extract → probe → route)
    METADATA = "metadata"
    ANNOTATION = "annotation"
    BAM = "bam"
    ARCHIVE = "archive"
    UNKNOWN = "unknown"


class Strategy(str, Enum):
    """Overall processing strategy for a dataset."""

    PROCESSED_COUNTS_FIRST = "processed_counts_first"
    RAW_FASTQ_ONLY = "raw_fastq_only"
    HYBRID = "hybrid"
    METADATA_ONLY = "metadata_only"
    NORMALIZED_ONLY = "normalized_only"  # Only normalized data available — limited use
    UNSOLVABLE = "unsolvable"


class FileGroup(BaseModel):
    """A group of related files that form one logical input."""

    group_id: str = Field(..., description="Unique group identifier within the dataset")
    group_type: FileGroupType = Field(..., description="Type of this file group")
    sample_id: str | None = Field(None, description="Sample this group belongs to")
    gsm: str | None = Field(None, description="GSM accession for this sample")

    files: dict[str, str] = Field(
        ...,
        description="Role -> file_id mapping (e.g. {'matrix': 'sha256:...', 'barcodes': 'sha256:...'})",
    )
    route: str = Field(..., description="Selected processing route for this group")
    route_description: str | None = Field(None, description="Human-readable route description")
    confidence: float = Field(..., ge=0.0, le=1.0, description="Group confidence score")

    notes: list[str] = Field(default_factory=list, description="Notes about this group")


class DatasetManifest(BaseModel):
    """Complete dataset manifest aggregating all file evidence and metadata."""

    dataset_id: str = Field(..., description="Primary accession, e.g. 'GSE123456'")
    title: str | None = Field(None, description="Dataset title from GEO")
    organism: list[str] = Field(default_factory=list, description="Organisms detected")
    platform: list[str] = Field(default_factory=list, description="Platforms used (e.g. '10x Genomics')")
    technology: str = Field(default="scRNA-seq", description="Sequencing technology")

    # Study context from web/database metadata probe
    tissue_or_organ: list[str] = Field(
        default_factory=list,
        description="Tissue or organ(s) studied (e.g. ['lung', 'liver'])",
    )
    disease: list[str] = Field(
        default_factory=list,
        description="Disease or condition studied (e.g. ['lung cancer'])",
    )
    condition_groups: list[str] = Field(
        default_factory=list,
        description="Experimental condition groups (e.g. ['tumor', 'normal', 'treated'])",
    )
    publication: dict[str, str] = Field(
        default_factory=dict,
        description="Publication metadata: {'pmid': '...', 'doi': '...'}",
    )
    metadata_confidence: dict[str, float] = Field(
        default_factory=dict,
        description="Confidence scores for metadata inferences: "
                    "{'organism': 0.98, 'condition_groups': 0.82, ...}",
    )

    files: dict[str, FileSidecar] = Field(
        default_factory=dict,
        description="File ID -> FileSidecar mapping for all files",
    )
    file_groups: list[FileGroup] = Field(
        default_factory=list,
        description="Logical groups of related files",
    )

    # Routing priority: explicit ordered list of group IDs to process first
    routing_priority: list[str] = Field(
        default_factory=list,
        description="Ordered list of group_ids in processing priority order. "
                    "Follows the rule: h5ad > 10x H5 > loom > MEX > raw UMI > "
                    "RDS > normalized-only > FASTQ.",
    )

    selected_strategy: Strategy = Field(..., description="Overall processing strategy")
    final_target: str = Field(default="h5ad", description="Target output format")

    unclassified_files: list[str] = Field(
        default_factory=list,
        description="File IDs that could not be classified",
    )
    manual_review_needed: bool = Field(
        default=False,
        description="Whether this dataset requires manual review before processing",
    )
    review_reasons: list[str] = Field(
        default_factory=list,
        description="Reasons manual review is needed",
    )

    manifest_version: str = Field(default="v1", description="Manifest schema version")
    created_at: datetime = Field(default_factory=lambda: datetime.now(timezone.utc), description="When this manifest was created")
    source_accessions: list[str] = Field(
        default_factory=list,
        description="All accessions (GSE, GSM, SRR) referenced in this dataset",
    )
