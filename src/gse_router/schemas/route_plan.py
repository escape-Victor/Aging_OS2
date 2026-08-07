"""Route plan schema.

Describes the selected processing plan: which groups to process,
in what order, with what methods and parameters.
"""

from __future__ import annotations

from datetime import datetime, timezone
from enum import Enum

from pydantic import BaseModel, Field


class RouteMethod(str, Enum):
    """Available processing methods."""

    # Raw FASTQ routes
    CELL_RANGER_COUNT = "cellranger_count"
    STAR_SOLO = "starsolo"
    KALLISTO_BUSTOOLS = "kallisto_bustools"
    ALEVIN_FRY = "alevin_fry"

    # Processed matrix routes
    READ_10X_MTX = "read_10x_mtx_to_h5ad"
    READ_10X_H5 = "read_10x_h5_to_h5ad"
    READ_H5AD = "read_h5ad_validate"
    READ_LOOM = "read_loom_to_h5ad"
    READ_CSV_MATRIX = "read_csv_matrix_to_h5ad"
    READ_RAW_UMI_MATRIX = "read_raw_umi_matrix_to_h5ad"
    READ_NORMALIZED_MATRIX = "read_normalized_matrix_to_h5ad"
    READ_RDS_SEURAT = "read_rds_seurat_to_h5ad"
    READ_RDS_SCE = "read_rds_sce_to_h5ad"

    # Metadata routes
    MERGE_CELL_METADATA = "merge_cell_metadata"
    MERGE_GENE_ANNOTATION = "merge_gene_annotation"

    # Special
    MANUAL_REVIEW = "manual_review"
    SKIP = "skip"


class ReferenceGenome(str, Enum):
    """Supported reference genomes."""

    GRCH38 = "GRCh38"
    GRCH37 = "GRCh37"
    GRCM39 = "GRCm39"
    GRCM38 = "GRCm38"
    CUSTOM = "custom"


class ProcessingStep(BaseModel):
    """A single processing step in the route plan."""

    step_id: str = Field(..., description="Unique step identifier")
    method: RouteMethod = Field(..., description="Processing method to use")
    input_group_id: str = Field(..., description="File group ID to process")
    output_label: str = Field(..., description="Label for the output of this step")

    parameters: dict[str, str | int | float | bool] = Field(
        default_factory=dict,
        description="Method-specific parameters",
    )
    reference_genome: ReferenceGenome | None = Field(
        None, description="Reference genome (for alignment steps)"
    )
    expected_runtime_minutes: int | None = Field(None, description="Estimated runtime")

    depends_on: list[str] = Field(
        default_factory=list,
        description="Step IDs that must complete before this step",
    )
    is_critical: bool = Field(default=True, description="Whether failure blocks the overall pipeline")
    retry_count: int = Field(default=0, description="Number of retries attempted")


class RoutePlan(BaseModel):
    """Complete route plan for processing a dataset."""

    dataset_id: str = Field(..., description="Dataset accession")
    manifest_version: str = Field(..., description="Version of the source manifest")

    strategy: str = Field(..., description="Overall processing strategy")
    reference_genome: ReferenceGenome | None = Field(None, description="Selected reference genome")
    organism_consensus: str | None = Field(
        None, description="Consensus organism determined from metadata + file evidence"
    )

    steps: list[ProcessingStep] = Field(
        default_factory=list,
        description="Ordered processing steps",
    )

    expected_outputs: list[str] = Field(
        default_factory=list,
        description="Expected output artifact labels",
    )
    quality_checks: list[str] = Field(
        default_factory=list,
        description="QC steps to run after processing",
    )

    low_confidence_steps: list[str] = Field(
        default_factory=list,
        description="Step IDs with confidence < 0.7 that need review",
    )
    warnings: list[str] = Field(
        default_factory=list,
        description="Warnings and caveats about this plan",
    )

    plan_version: str = Field(default="v1", description="Route plan schema version")
    created_at: datetime = Field(default_factory=lambda: datetime.now(timezone.utc), description="When this plan was generated")
