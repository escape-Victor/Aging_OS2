"""Routing rules — biological judgment encoded as policy.

These rules embody the biological knowledge that determines which
processing route to choose, what parameters to use, and when to stop
and ask for manual review.

Rules are functions that take a manifest and return decisions —
they do NOT execute any processing themselves.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Callable

from ..schemas.manifest import DatasetManifest, FileGroup, FileGroupType
from ..schemas.route_plan import ReferenceGenome, RouteMethod


class RulePriority(str, Enum):
    """Priority level of a routing rule."""

    CRITICAL = "critical"  # Must pass — violation stops pipeline
    HIGH = "high"          # Strong recommendation
    MEDIUM = "medium"      # Default guidance
    LOW = "low"            # Soft preference


@dataclass
class RuleResult:
    """Result of applying a routing rule."""

    rule_name: str
    passed: bool
    priority: RulePriority
    message: str
    action: str = ""  # What to do about it


@dataclass
class RouteDecision:
    """A decoded routing decision for one file group."""

    group_id: str
    method: RouteMethod
    reference_genome: ReferenceGenome | None = None
    parameters: dict[str, str | int | float | bool] = field(default_factory=dict)
    warnings: list[str] = field(default_factory=list)


# ---------------------------------------------------------------------------
# Rule registry
# ---------------------------------------------------------------------------

# Each rule is a callable: (manifest) -> RuleResult
_rules: list[Callable[[DatasetManifest], RuleResult]] = []


def register_rule(priority: RulePriority):
    """Decorator to register a routing rule function."""
    def decorator(fn: Callable[[DatasetManifest], RuleResult]) -> Callable[[DatasetManifest], RuleResult]:
        _rules.append(fn)
        return fn
    return decorator


def apply_all_rules(manifest: DatasetManifest) -> list[RuleResult]:
    """Apply all registered rules to a manifest and return results."""
    return [rule(manifest) for rule in _rules]


def get_failures(results: list[RuleResult]) -> list[RuleResult]:
    """Return only failed critical/high rules that should block processing."""
    return [r for r in results if not r.passed and r.priority in (RulePriority.CRITICAL, RulePriority.HIGH)]


# ---------------------------------------------------------------------------
# Core routing rules
# ---------------------------------------------------------------------------

@register_rule(RulePriority.CRITICAL)
def rule_must_have_usable_data(manifest: DatasetManifest) -> RuleResult:
    """There must be at least one processable file group."""
    processable = {
        FileGroupType.H5AD, FileGroupType.H5_10X, FileGroupType.LOOM,
        FileGroupType.MEX_10X, FileGroupType.ARCHIVE_10X,
        FileGroupType.RAW_UMI_MATRIX,
        FileGroupType.SEURAT_RDS, FileGroupType.SCE_RDS, FileGroupType.H5SEURAT,
        FileGroupType.CSV_MATRIX, FileGroupType.NORMALIZED_MATRIX,
        FileGroupType.FASTQ_PAIRED, FileGroupType.FASTQ_SINGLE,
    }
    has_data = any(g.group_type in processable for g in manifest.file_groups)
    # Check for limited-use only (must have at least one processable group)
    processable_groups = [
        g for g in manifest.file_groups if g.group_type in processable
    ]
    only_normalized = (
        len(processable_groups) > 0
        and all(g.group_type == FileGroupType.NORMALIZED_MATRIX for g in processable_groups)
    )
    if only_normalized:
        return RuleResult(
            rule_name="must_have_usable_data",
            passed=True,
            priority=RulePriority.CRITICAL,
            message="Only normalized data available — will be marked as limited-use",
            action="Process but mark output as 'normalized_only, no raw counts'",
        )
    return RuleResult(
        rule_name="must_have_usable_data",
        passed=has_data,
        priority=RulePriority.CRITICAL,
        message="Dataset has processable data" if has_data else "No processable data found in dataset",
        action="" if has_data else "Stop — set strategy to UNSOLVABLE and report for manual review",
    )


@register_rule(RulePriority.CRITICAL)
def rule_organism_must_be_determined(manifest: DatasetManifest) -> RuleResult:
    """An organism must be determined before processing."""
    has_organism = len(manifest.organism) > 0
    return RuleResult(
        rule_name="organism_must_be_determined",
        passed=has_organism,
        priority=RulePriority.CRITICAL,
        message=f"Organism: {manifest.organism}" if has_organism else "Organism not determined",
        action="" if has_organism else "Run organism resolution before processing",
    )


@register_rule(RulePriority.CRITICAL)
def rule_skip_fastq_if_processed_exists(manifest: DatasetManifest) -> RuleResult:
    """If any usable processed count/object file exists, do NOT run FASTQ/SRA processing.

    This is the cardinal rule from the GEO examples research:
    'RAW.tar' does NOT mean FASTQ — many GEO RAW.tar files contain processed matrices.
    Always prefer the highest-priority processed data available.
    """
    processed_types = {
        FileGroupType.H5AD, FileGroupType.H5_10X, FileGroupType.LOOM,
        FileGroupType.MEX_10X, FileGroupType.ARCHIVE_10X,
        FileGroupType.RAW_UMI_MATRIX,
        FileGroupType.SEURAT_RDS, FileGroupType.SCE_RDS, FileGroupType.H5SEURAT,
        FileGroupType.CSV_MATRIX,
    }
    has_processed = any(g.group_type in processed_types for g in manifest.file_groups)
    has_raw = any(
        g.group_type in (FileGroupType.FASTQ_PAIRED, FileGroupType.FASTQ_SINGLE)
        for g in manifest.file_groups
    )
    if has_processed and has_raw:
        raw_groups = [g.group_id for g in manifest.file_groups
                      if g.group_type in (FileGroupType.FASTQ_PAIRED, FileGroupType.FASTQ_SINGLE)]
        return RuleResult(
            rule_name="skip_fastq_if_processed_exists",
            passed=True,
            priority=RulePriority.CRITICAL,
            message=f"Processed data available — skipping FASTQ/SRA groups: {raw_groups}",
            action="Exclude FASTQ groups from processing. Only use FASTQ if ALL processed routes fail.",
        )
    return RuleResult(
        rule_name="skip_fastq_if_processed_exists",
        passed=True,
        priority=RulePriority.CRITICAL,
        message="No FASTQ to skip or no processed alternative available",
        action="",
    )


@register_rule(RulePriority.CRITICAL)
def rule_raw_tar_inspection(manifest: DatasetManifest) -> RuleResult:
    """Before downloading RAW.tar, check filelist.txt or FTP index first.

    GEO RAW.tar archives often contain processed matrices, .loom, .h5ad,
    or metadata — not necessarily FASTQ. Checking the file listing first
    prevents unnecessary downloads and misrouting.
    """
    has_raw_tar = any(
        g.group_type == FileGroupType.ARCHIVE
        and any(
            sc.detected_type.value == "raw_tar"
            for sc in manifest.files.values()
            if sc.file_id in g.files.values()
        )
        for g in manifest.file_groups
    )
    has_filelist = any(
        sc.detected_type.value == "filelist_txt"
        for sc in manifest.files.values()
    )
    if has_raw_tar and not has_filelist:
        return RuleResult(
            rule_name="raw_tar_inspection",
            passed=False,
            priority=RulePriority.CRITICAL,
            message="RAW.tar detected but no filelist.txt found — contents unknown",
            action="Fetch filelist.txt from GEO FTP before downloading RAW.tar. "
                    "RAW.tar may contain processed data, not FASTQ.",
        )
    return RuleResult(
        rule_name="raw_tar_inspection",
        passed=True,
        priority=RulePriority.CRITICAL,
        message="RAW.tar inspection OK or no RAW.tar present",
        action="",
    )


@register_rule(RulePriority.HIGH)
def rule_prefer_processed_over_raw(manifest: DatasetManifest) -> RuleResult:
    """Prefer processed count matrices over raw FASTQ when both are available."""
    processed_types = {
        FileGroupType.H5AD, FileGroupType.H5_10X, FileGroupType.LOOM,
        FileGroupType.MEX_10X, FileGroupType.ARCHIVE_10X,
        FileGroupType.RAW_UMI_MATRIX,
        FileGroupType.SEURAT_RDS, FileGroupType.SCE_RDS, FileGroupType.H5SEURAT,
        FileGroupType.CSV_MATRIX,
    }
    has_processed = any(g.group_type in processed_types for g in manifest.file_groups)
    has_raw = any(
        g.group_type in (FileGroupType.FASTQ_PAIRED, FileGroupType.FASTQ_SINGLE)
        for g in manifest.file_groups
    )
    if has_processed and has_raw:
        return RuleResult(
            rule_name="prefer_processed_over_raw",
            passed=True,
            priority=RulePriority.HIGH,
            message="Both processed and raw data available — preferring processed",
            action="Route processed groups first, skip raw FASTQ unless processed fails",
        )
    return RuleResult(
        rule_name="prefer_processed_over_raw",
        passed=True,
        priority=RulePriority.HIGH,
        message="Single data type available",
        action="Use whichever data is available",
    )


@register_rule(RulePriority.HIGH)
def rule_strict_routing_priority(manifest: DatasetManifest) -> RuleResult:
    """Apply the explicit 8-level routing priority.

    Priority order (highest to lowest):
      1. .h5ad -> validate/standardize
      2. 10x .h5 -> read_10x_h5
      3. .loom / .loom.gz -> read_loom
      4. 10x MEX triplet -> read_10x_mtx
      5. raw UMI TXT/CSV/RDS matrix -> parse/convert
      6. Seurat/SCE .rds -> R bridge conversion
      7. normalized matrix only -> mark limited-use
      8. FASTQ/SRA -> avoid unless no processed data exists
    """
    priority_map = {
        FileGroupType.H5AD: 1,
        FileGroupType.H5_10X: 2,
        FileGroupType.LOOM: 3,
        FileGroupType.MEX_10X: 4,
        FileGroupType.ARCHIVE_10X: 4,  # Contains 10x MEX — same priority as MEX_10X
        FileGroupType.RAW_UMI_MATRIX: 5,
        FileGroupType.SEURAT_RDS: 6,
        FileGroupType.SCE_RDS: 6,
        FileGroupType.H5SEURAT: 6,
        FileGroupType.CSV_MATRIX: 6,
        FileGroupType.NORMALIZED_MATRIX: 7,
        FileGroupType.FASTQ_PAIRED: 8,
        FileGroupType.FASTQ_SINGLE: 8,
    }
    data_groups = [
        g for g in manifest.file_groups
        if g.group_type not in (FileGroupType.METADATA, FileGroupType.ANNOTATION,
                                FileGroupType.ARCHIVE, FileGroupType.UNKNOWN)
    ]
    if not data_groups:
        return RuleResult(
            rule_name="strict_routing_priority",
            passed=True,
            priority=RulePriority.HIGH,
            message="No data groups to prioritize",
            action="",
        )
    best = min(data_groups, key=lambda g: priority_map.get(g.group_type, 99))
    best_priority = priority_map.get(best.group_type, 99)
    has_fastq = any(g.group_type in (FileGroupType.FASTQ_PAIRED, FileGroupType.FASTQ_SINGLE) for g in data_groups)
    if has_fastq and best_priority < 8:
        return RuleResult(
            rule_name="strict_routing_priority",
            passed=True,
            priority=RulePriority.HIGH,
            message=f"Best available: {best.group_type.value} (priority {best_priority}). "
                    f"FASTQ present but will be skipped.",
            action="Process priority 1-7 groups; skip FASTQ unless all higher-priority routes fail.",
        )
    return RuleResult(
        rule_name="strict_routing_priority",
        passed=True,
        priority=RulePriority.HIGH,
        message=f"Processing with best available: {best.group_type.value} (priority {best_priority})",
        action="",
    )


@register_rule(RulePriority.HIGH)
def rule_warn_normalized_as_counts(manifest: DatasetManifest) -> RuleResult:
    """Warn if the only available data appears to be normalized, not raw counts."""
    # Check for normalized-only matrices
    normalized_groups = [
        g for g in manifest.file_groups
        if g.group_type == FileGroupType.NORMALIZED_MATRIX
    ]
    # Also check for suspected normalized data in h5ad or CSV
    suspect_groups = [
        g for g in manifest.file_groups
        if g.group_type in (FileGroupType.H5AD, FileGroupType.CSV_MATRIX)
        and g.confidence < 0.85
    ]
    # Check per-file is_raw_counts flag
    normalized_files = [
        fid for fid, sc in manifest.files.items()
        if sc.is_raw_counts is False
    ]

    if normalized_groups or normalized_files:
        return RuleResult(
            rule_name="warn_normalized_as_counts",
            passed=True,
            priority=RulePriority.HIGH,
            message=f"Normalized data detected: {len(normalized_groups)} groups, "
                    f"{len(normalized_files)} files flagged as non-raw. "
                    f"Do NOT treat normalized expression as raw counts.",
            action="If raw counts unavailable: mark output as 'limited-use: normalized only'. "
                    "Preserve in layers['normalized'], do NOT copy to layers['counts'].",
        )
    if suspect_groups:
        return RuleResult(
            rule_name="warn_normalized_as_counts",
            passed=True,
            priority=RulePriority.HIGH,
            message=f"{len(suspect_groups)} groups may contain normalized data, not raw counts",
            action="Flag for QC: verify these contain raw counts; do not treat normalized as raw",
        )
    return RuleResult(
        rule_name="warn_normalized_as_counts",
        passed=True,
        priority=RulePriority.HIGH,
        message="No suspected normalized-as-counts issues",
        action="",
    )


@register_rule(RulePriority.MEDIUM)
def rule_check_metadata_available(manifest: DatasetManifest) -> RuleResult:
    """Check whether cell/sample metadata is available to attach."""
    meta_groups = [g for g in manifest.file_groups if g.group_type == FileGroupType.METADATA]
    if meta_groups:
        return RuleResult(
            rule_name="check_metadata_available",
            passed=True,
            priority=RulePriority.MEDIUM,
            message=f"{len(meta_groups)} metadata groups found",
            action="Plan metadata merge steps for each matched sample group",
        )
    return RuleResult(
        rule_name="check_metadata_available",
        passed=True,
        priority=RulePriority.MEDIUM,
        message="No auxiliary metadata found — using GEO sample annotations only",
        action="",
    )


@register_rule(RulePriority.LOW)
def rule_confidence_threshold(manifest: DatasetManifest) -> RuleResult:
    """Flag low-confidence groups for review."""
    low_conf = [g for g in manifest.file_groups if g.confidence < 0.7]
    if low_conf:
        names = ", ".join(g.group_id for g in low_conf)
        return RuleResult(
            rule_name="confidence_threshold",
            passed=False,
            priority=RulePriority.LOW,
            message=f"{len(low_conf)} groups below confidence threshold: {names}",
            action="Mark for manual review before processing",
        )
    return RuleResult(
        rule_name="confidence_threshold",
        passed=True,
        priority=RulePriority.LOW,
        message="All groups above confidence threshold",
        action="",
    )


# ---------------------------------------------------------------------------
# Route planning helpers
# ---------------------------------------------------------------------------

def select_reference_genome(organism: str) -> ReferenceGenome | None:
    """Map organism name to reference genome."""
    organism_lower = organism.lower().strip()
    human_terms = {"homo sapiens", "human", "hsapiens"}
    mouse_terms = {"mus musculus", "mouse", "mmusculus"}

    if organism_lower in human_terms:
        return ReferenceGenome.GRCH38
    if organism_lower in mouse_terms:
        return ReferenceGenome.GRCM39
    if "homo" in organism_lower or "human" in organism_lower:
        return ReferenceGenome.GRCH38
    if "mus" in organism_lower or "mouse" in organism_lower:
        return ReferenceGenome.GRCM39
    return ReferenceGenome.CUSTOM


def decode_route(group: FileGroup) -> RouteDecision:
    """Map a file group's route string to a concrete RouteDecision.

    This is the bridge between the classifier's route strings and the
    route planner's typed RouteMethod + parameters.
    """
    route_map: dict[str, RouteDecision] = {
        # Priority 1: existing h5ad
        "read_h5ad_validate": RouteDecision(
            group_id=group.group_id,
            method=RouteMethod.READ_H5AD,
        ),
        # Priority 2: 10x H5
        "read_10x_h5_to_h5ad": RouteDecision(
            group_id=group.group_id,
            method=RouteMethod.READ_10X_H5,
        ),
        # Priority 3: Loom
        "read_loom_to_h5ad": RouteDecision(
            group_id=group.group_id,
            method=RouteMethod.READ_LOOM,
        ),
        # Priority 4: 10x MEX
        "read_10x_mtx_to_h5ad": RouteDecision(
            group_id=group.group_id,
            method=RouteMethod.READ_10X_MTX,
            parameters={"var_names": "gene_symbols"},
        ),
        # Priority 5: raw UMI matrix
        "read_raw_umi_matrix_to_h5ad": RouteDecision(
            group_id=group.group_id,
            method=RouteMethod.READ_RAW_UMI_MATRIX,
            parameters={"var_names": "gene_symbols"},
        ),
        # Priority 6: RDS (Seurat/SCE)
        "read_rds_seurat_to_h5ad": RouteDecision(
            group_id=group.group_id,
            method=RouteMethod.READ_RDS_SEURAT,
        ),
        "read_rds_sce_to_h5ad": RouteDecision(
            group_id=group.group_id,
            method=RouteMethod.READ_RDS_SCE,
        ),
        # Generic CSV matrix (ambiguous)
        "read_csv_matrix_to_h5ad": RouteDecision(
            group_id=group.group_id,
            method=RouteMethod.READ_CSV_MATRIX,
        ),
        # Priority 7: normalized only — limited use
        "read_normalized_matrix_to_h5ad": RouteDecision(
            group_id=group.group_id,
            method=RouteMethod.READ_NORMALIZED_MATRIX,
            parameters={"normalization_method": "unknown"},
        ),
        # Priority 8: FASTQ/SRA (avoid if processed exists)
        "cellranger_count": RouteDecision(
            group_id=group.group_id,
            method=RouteMethod.CELL_RANGER_COUNT,
        ),
        "starsolo": RouteDecision(
            group_id=group.group_id,
            method=RouteMethod.STAR_SOLO,
        ),
        # Archive extraction (10x bundles)
        "extract_10x_archive_to_h5ad": RouteDecision(
            group_id=group.group_id,
            method=RouteMethod.READ_10X_MTX,  # Extract → probe → treat as MEX
            parameters={"source": "archive_10x"},
        ),
        # Metadata
        "merge_cell_metadata": RouteDecision(
            group_id=group.group_id,
            method=RouteMethod.MERGE_CELL_METADATA,
        ),
        "merge_gene_annotation": RouteDecision(
            group_id=group.group_id,
            method=RouteMethod.MERGE_GENE_ANNOTATION,
        ),
    }

    decision = route_map.get(group.route)
    if decision is None:
        return RouteDecision(
            group_id=group.group_id,
            method=RouteMethod.MANUAL_REVIEW,
            warnings=[f"No route mapping for: {group.route}"],
        )
    return decision
