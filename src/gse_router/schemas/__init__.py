"""Schema definitions for GSE Router data structures."""

from .sidecar import (
    FileDetectedType,
    FileSidecar,
    ProcessingStage,
    RouteRecommendation,
)
from .manifest import (
    DatasetManifest,
    FileGroup,
    FileGroupType,
    Strategy,
)
from .route_plan import (
    ProcessingStep,
    ReferenceGenome,
    RouteMethod,
    RoutePlan,
)
from .report import (
    AnnDataSummary,
    ProcessingReport,
    QCMetric,
    StepResult,
    StepStatus,
)

__all__ = [
    # Sidecar
    "FileDetectedType",
    "FileSidecar",
    "ProcessingStage",
    "RouteRecommendation",
    # Manifest
    "DatasetManifest",
    "FileGroup",
    "FileGroupType",
    "Strategy",
    # Route plan
    "ProcessingStep",
    "ReferenceGenome",
    "RouteMethod",
    "RoutePlan",
    # Report
    "AnnDataSummary",
    "ProcessingReport",
    "QCMetric",
    "StepResult",
    "StepStatus",
]
