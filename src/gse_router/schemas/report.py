"""Processing report schema.

Final report recording provenance, decisions, QC results, and output metadata.
"""

from __future__ import annotations

from datetime import datetime, timezone
from enum import Enum

from pydantic import BaseModel, Field


class StepStatus(str, Enum):
    """Status of a processing step."""

    PENDING = "pending"
    RUNNING = "running"
    COMPLETED = "completed"
    FAILED = "failed"
    SKIPPED = "skipped"
    MANUAL_REVIEW = "manual_review"


class QCMetric(BaseModel):
    """A single QC metric."""

    name: str = Field(..., description="Metric name, e.g. 'n_genes_by_counts'")
    value: float | int | str = Field(..., description="Metric value")
    threshold: str | None = Field(None, description="Expected threshold, e.g. '>200'")
    passed: bool | None = Field(None, description="Whether the metric passed QC")
    notes: str | None = Field(None, description="Additional notes")


class StepResult(BaseModel):
    """Result of a single processing step."""

    step_id: str = Field(..., description="Step ID from the route plan")
    method: str = Field(..., description="Method used")
    status: StepStatus = Field(..., description="Execution status")
    started_at: datetime | None = Field(None, description="When this step started")
    finished_at: datetime | None = Field(None, description="When this step finished")
    duration_seconds: float | None = Field(None, description="Wall-clock duration")

    input_artifacts: list[str] = Field(default_factory=list, description="Input file IDs / labels")
    output_artifacts: list[str] = Field(default_factory=list, description="Output file labels / paths")

    error_message: str | None = Field(None, description="Error message if failed")
    warnings: list[str] = Field(default_factory=list, description="Warnings from this step")
    logs_path: str | None = Field(None, description="Path to detailed logs for this step")


class AnnDataSummary(BaseModel):
    """Summary statistics for the output AnnData object."""

    n_obs: int = Field(..., description="Number of cells")
    n_vars: int = Field(..., description="Number of genes/features")
    n_layers: int = Field(default=0, description="Number of additional layers")
    layer_names: list[str] = Field(default_factory=list, description="Names of layers")
    obs_columns: list[str] = Field(default_factory=list, description="Column names in obs")
    var_columns: list[str] = Field(default_factory=list, description="Column names in var")
    has_raw_counts: bool = Field(default=False, description="Whether raw counts are preserved")
    raw_counts_layer: str | None = Field(None, description="Layer name for raw counts")
    file_path: str | None = Field(None, description="Path to the output .h5ad file")
    file_size_mb: float | None = Field(None, description="Output file size in MB")


class ProcessingReport(BaseModel):
    """Complete processing report with full provenance."""

    dataset_id: str = Field(..., description="Dataset accession")
    report_id: str = Field(..., description="Unique report identifier")

    status: StepStatus = Field(..., description="Overall processing status")
    strategy: str = Field(..., description="Strategy used")
    reference_genome: str | None = Field(None, description="Reference genome used")

    steps: list[StepResult] = Field(default_factory=list, description="Results of each processing step")
    qc_metrics: list[QCMetric] = Field(default_factory=list, description="QC metrics and their results")

    output: AnnDataSummary | None = Field(None, description="Output object summary")
    provenance: dict[str, object] = Field(
        default_factory=dict,
        description="Full provenance chain: inputs -> tools -> outputs",
    )

    unprocessed_groups: list[str] = Field(
        default_factory=list,
        description="File groups that could not be processed",
    )
    manual_review_items: list[str] = Field(
        default_factory=list,
        description="Items requiring manual review",
    )

    total_duration_seconds: float | None = Field(None, description="Total wall-clock time")
    peak_memory_mb: float | None = Field(None, description="Peak memory usage")

    warnings: list[str] = Field(default_factory=list, description="Global warnings")
    errors: list[str] = Field(default_factory=list, description="Global errors")

    report_version: str = Field(default="v1", description="Report schema version")
    created_at: datetime = Field(default_factory=lambda: datetime.now(timezone.utc), description="When this report was generated")

    source_manifest_id: str | None = Field(None, description="ID of the source manifest")
    source_route_plan_id: str | None = Field(None, description="ID of the source route plan")
