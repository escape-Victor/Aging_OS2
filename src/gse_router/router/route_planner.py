"""Biology-aware route planner.

Combines web metadata + local file evidence to produce a detailed
route plan: which groups to process, in what order, with what methods
and parameters.
"""

from __future__ import annotations


import structlog

from ..probes.metadata import StudyContext
from ..schemas.manifest import DatasetManifest, FileGroup, FileGroupType, Strategy
from ..schemas.route_plan import ProcessingStep, ReferenceGenome, RouteMethod, RoutePlan
from .rules import (
    apply_all_rules,
    decode_route,
    get_failures,
    select_reference_genome,
)

logger = structlog.get_logger(__name__)


def _any_mex_in_extracted(
    groups: list[FileGroup],
    sidecars: dict | None,
) -> bool:
    """Return True if any MEX_10X group has files inside an ``_extracted_*`` dir.

    When this is True, the route planner should skip ARCHIVE_10X groups because
    the archive's contents are already available as extracted MEX triplets.
    Processing both would convert the same data twice.
    """
    if not sidecars:
        return False
    for group in groups:
        if group.group_type != FileGroupType.MEX_10X:
            continue
        for file_id in group.files.values():
            sc = sidecars.get(file_id)
            if sc is not None and "_extracted_" in str(sc.path):
                return True
    return False


class RoutePlanner:
    """The third intelligence layer: combines metadata + file evidence
    to produce a concrete, ordered, parameterized processing plan.

    This planner:
    1. Applies routing rules to the manifest (biological judgment)
    2. Decodes each file group into a specific RouteMethod + parameters
    3. Computes dependencies between steps
    4. Selects reference genome and QC checks
    5. Produces a RoutePlan ready for execution
    """

    def __init__(self) -> None:
        self._logger = logger.bind(component="RoutePlanner")

    def plan(
        self,
        manifest: DatasetManifest,
        study_context: StudyContext | None = None,
    ) -> RoutePlan:
        """Generate a route plan from a manifest and optional study context.

        Args:
            manifest: The complete dataset manifest.
            study_context: Optional web-derived study metadata.

        Returns:
            A complete RoutePlan ready for execution.
        """
        self._logger.info("Planning route", dataset_id=manifest.dataset_id)

        # 1. Apply routing rules
        rule_results = apply_all_rules(manifest)
        failures = get_failures(rule_results)

        # 2. Select reference genome
        ref_genome = self._resolve_reference_genome(manifest, study_context)

        # 3. Split groups by processing order
        data_groups = self._prioritize_groups(
            manifest.file_groups, manifest.selected_strategy,
            sidecars=manifest.files,
        )
        meta_groups = [g for g in manifest.file_groups if g.route == "merge_cell_metadata"]
        anno_groups = [g for g in manifest.file_groups if g.route == "merge_gene_annotation"]

        # 4. Build steps
        steps: list[ProcessingStep] = []
        step_counter = 0

        # Phase 1: Process data groups (primary matrices/objects)
        data_steps = self._build_data_steps(data_groups, ref_genome, start_counter=step_counter)
        steps.extend(data_steps)
        step_counter += len(data_steps)

        # Collect data step IDs so metadata/annotation steps can declare
        # their dependency on h5ad outputs from prior data conversion steps.
        data_step_ids = [s.step_id for s in data_steps]

        # Phase 2: Merge metadata
        meta_steps = self._build_metadata_steps(
            meta_groups, start_counter=step_counter, data_step_ids=data_step_ids,
        )
        steps.extend(meta_steps)
        step_counter += len(meta_steps)

        # Phase 3: Merge annotations
        anno_steps = self._build_annotation_steps(
            anno_groups, start_counter=step_counter, data_step_ids=data_step_ids,
        )
        steps.extend(anno_steps)
        step_counter += len(anno_steps)

        # 5. Collect warnings and low-confidence items
        warnings = [r.message for r in failures]
        low_conf_steps = [
            s.step_id for s in steps
            if any(g.confidence < 0.7 for g in manifest.file_groups if g.group_id == s.input_group_id)
        ]

        # 6. Build QC checklist
        qc_checks = self._build_qc_checklist(manifest, ref_genome)

        plan = RoutePlan(
            dataset_id=manifest.dataset_id,
            manifest_version=manifest.manifest_version,
            strategy=manifest.selected_strategy.value,
            reference_genome=ref_genome,
            organism_consensus=manifest.organism[0] if manifest.organism else None,
            steps=steps,
            expected_outputs=["final_h5ad", "provenance_report"],
            quality_checks=qc_checks,
            low_confidence_steps=low_conf_steps,
            warnings=warnings,
        )
        self._logger.info(
            "Route plan generated",
            dataset_id=manifest.dataset_id,
            step_count=len(steps),
            strategy=manifest.selected_strategy.value,
        )
        return plan

    # ------------------------------------------------------------------
    # Step builders
    # ------------------------------------------------------------------

    def _build_data_steps(
        self,
        groups: list[FileGroup],
        ref_genome: ReferenceGenome | None,
        start_counter: int = 0,
    ) -> list[ProcessingStep]:
        """Build processing steps for data groups (matrices, objects, FASTQ)."""
        steps: list[ProcessingStep] = []
        for i, group in enumerate(groups):
            decision = decode_route(group)
            step = ProcessingStep(
                step_id=f"step_{start_counter + i:03d}",
                method=decision.method,
                input_group_id=group.group_id,
                output_label=f"h5ad_{group.group_id}",
                parameters=decision.parameters,
                reference_genome=decision.reference_genome or ref_genome,
                expected_runtime_minutes=self._estimate_runtime(decision.method),
                depends_on=[],
                is_critical=decision.method != RouteMethod.SKIP,
            )
            if decision.warnings:
                step.parameters["_warnings"] = "; ".join(decision.warnings)
            steps.append(step)
        return steps

    def _build_metadata_steps(
        self,
        groups: list[FileGroup],
        start_counter: int = 0,
        data_step_ids: list[str] | None = None,
    ) -> list[ProcessingStep]:
        """Build steps for merging cell/sample metadata.

        Metadata steps depend on all prior data steps for their h5ad inputs.
        """
        steps: list[ProcessingStep] = []
        for i, group in enumerate(groups):
            step = ProcessingStep(
                step_id=f"step_{start_counter + i:03d}",
                method=RouteMethod.MERGE_CELL_METADATA,
                input_group_id=group.group_id,
                output_label=f"metadata_{group.group_id}",
                parameters={},
                depends_on=data_step_ids or [],
                is_critical=False,
            )
            steps.append(step)
        return steps

    def _build_annotation_steps(
        self,
        groups: list[FileGroup],
        start_counter: int = 0,
        data_step_ids: list[str] | None = None,
    ) -> list[ProcessingStep]:
        """Build steps for merging gene annotations.

        Annotation steps depend on all prior data steps for their h5ad inputs.
        """
        steps: list[ProcessingStep] = []
        for i, group in enumerate(groups):
            step = ProcessingStep(
                step_id=f"step_{start_counter + i:03d}",
                method=RouteMethod.MERGE_GENE_ANNOTATION,
                input_group_id=group.group_id,
                output_label=f"annotation_{group.group_id}",
                parameters={},
                depends_on=data_step_ids or [],
                is_critical=False,
            )
            steps.append(step)
        return steps

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------

    def _prioritize_groups(
        self,
        groups: list[FileGroup],
        strategy: Strategy,
        sidecars: dict[str, "FileSidecar"] | None = None,
    ) -> list[FileGroup]:
        """Order groups for processing using the explicit 8-level routing priority.

        Priority (from GEO examples research):
          1. h5ad -> validate/standardize
          2. 10x H5 -> read_10x_h5
          3. loom/loom.gz -> read_loom
          4. 10x MEX triplet -> read_10x_mtx
          5. raw UMI TXT/CSV/RDS matrix -> parse/convert (confirmed raw counts)
          6. Seurat/SCE .rds -> R bridge conversion
          7. CSV/TXT matrix (ambiguous: may be raw or normalized)
          8. normalized matrix only -> mark limited-use
          9. FASTQ/SRA -> avoid unless no processed data exists

        RAW_UMI_MATRIX is distinct from CSV_MATRIX: RAW_UMI_MATRIX is explicitly
        confirmed raw counts and gets higher priority. CSV_MATRIX is ambiguous —
        the probe could not confirm whether counts are raw or normalized — and is
        treated at a lower priority tier alongside Seurat/RDS objects.
        NORMALIZED_MATRIX is explicitly confirmed normalized (non-raw) data and
        is placed just above FASTQ (last resort).
        """
        priority_order: dict[FileGroupType, int] = {
            FileGroupType.H5AD: 0,              # 1
            FileGroupType.H5_10X: 1,             # 2
            FileGroupType.LOOM: 2,               # 3
            FileGroupType.MEX_10X: 3,            # 4
            FileGroupType.ARCHIVE_10X: 3,        # 4 (contains 10x MEX — extract then same as MEX_10X)
            FileGroupType.RAW_UMI_MATRIX: 4,     # 5 (confirmed raw counts)
            FileGroupType.SEURAT_RDS: 5,         # 6
            FileGroupType.SCE_RDS: 5,            # 6
            FileGroupType.H5SEURAT: 5,           # 6
            FileGroupType.CSV_MATRIX: 5,         # 6 (ambiguous: raw or normalized)
            FileGroupType.NORMALIZED_MATRIX: 6,   # 7 (confirmed normalized — limited use)
            FileGroupType.FASTQ_PAIRED: 7,       # 8
            FileGroupType.FASTQ_SINGLE: 7,       # 8
            FileGroupType.METADATA: 10,
            FileGroupType.ANNOTATION: 11,
            FileGroupType.BAM: 12,
            FileGroupType.UNKNOWN: 99,
        }
        data_groups = [
            g for g in groups
            if g.group_type
            not in (FileGroupType.METADATA, FileGroupType.ANNOTATION,
                    FileGroupType.ARCHIVE, FileGroupType.UNKNOWN)
            or g.group_type == FileGroupType.ARCHIVE_10X
        ]
        if strategy == Strategy.NORMALIZED_ONLY:
            data_groups = [
                g for g in data_groups
                if g.group_type == FileGroupType.NORMALIZED_MATRIX
            ]
        elif strategy == Strategy.RAW_FASTQ_ONLY:
            data_groups = [
                g for g in data_groups
                if g.group_type in (FileGroupType.FASTQ_PAIRED, FileGroupType.FASTQ_SINGLE)
            ]
        elif strategy == Strategy.PROCESSED_COUNTS_FIRST:
            # Exclude FASTQ if ANY processed data exists (priority 0-6)
            has_processed = any(
                priority_order.get(g.group_type, 99) <= 6
                for g in data_groups
            )
            if has_processed:
                data_groups = [
                    g for g in data_groups
                    if g.group_type not in (FileGroupType.FASTQ_PAIRED, FileGroupType.FASTQ_SINGLE)
                ]

        # Deduplicate: if MEX triplets already exist in _extracted_* directories,
        # skip ARCHIVE_10X groups for the corresponding archives. Without this,
        # the planner creates steps for both, causing duplicate conversion of
        # the same underlying data (once directly, once via archive extraction).
        mex_from_extracted = _any_mex_in_extracted(data_groups, sidecars)
        if mex_from_extracted:
            data_groups = [
                g for g in data_groups
                if g.group_type != FileGroupType.ARCHIVE_10X
            ]

        data_groups.sort(key=lambda g: priority_order.get(g.group_type, 50))
        return data_groups

    def _resolve_reference_genome(
        self,
        manifest: DatasetManifest,
        study_context: StudyContext | None = None,
    ) -> ReferenceGenome | None:
        """Resolve the reference genome from manifest and study context."""
        organisms = manifest.organism.copy()
        if study_context and study_context.organism:
            organisms.extend(study_context.organism)

        if not organisms:
            return None

        # Use the most confident organism
        return select_reference_genome(organisms[0])

    def _estimate_runtime(self, method: RouteMethod) -> int | None:
        """Estimate runtime in minutes for a method (rough heuristic)."""
        estimates: dict[RouteMethod, int] = {
            RouteMethod.CELL_RANGER_COUNT: 120,
            RouteMethod.STAR_SOLO: 90,
            RouteMethod.KALLISTO_BUSTOOLS: 30,
            RouteMethod.ALEVIN_FRY: 30,
            RouteMethod.READ_10X_MTX: 2,
            RouteMethod.READ_10X_H5: 1,
            RouteMethod.READ_H5AD: 1,
            RouteMethod.READ_LOOM: 1,
            RouteMethod.READ_CSV_MATRIX: 5,
            RouteMethod.READ_RAW_UMI_MATRIX: 5,
            RouteMethod.READ_NORMALIZED_MATRIX: 2,
            RouteMethod.READ_RDS_SEURAT: 3,
            RouteMethod.READ_RDS_SCE: 3,
            RouteMethod.MERGE_CELL_METADATA: 1,
            RouteMethod.MERGE_GENE_ANNOTATION: 1,
        }
        return estimates.get(method)

    def _build_qc_checklist(
        self, manifest: DatasetManifest, ref_genome: ReferenceGenome | None
    ) -> list[str]:
        """Build a QC checklist based on the dataset context."""
        checks: list[str] = [
            "n_genes_by_counts > 200",
            "total_counts > 500",
            "pct_mitochondrial < 20",
            "doublet_detection",
            "n_cells >= 100",
        ]
        # Add organism-specific checks
        if ref_genome == ReferenceGenome.GRCH38:
            checks.append("verify_human_gene_symbols")
        elif ref_genome == ReferenceGenome.GRCM39:
            checks.append("verify_mouse_gene_symbols")
        checks.append("preserve_raw_counts_in_layer")
        checks.append("validate_obs_metadata_attached")
        checks.append("validate_var_metadata_attached")

        # Add matrix-type-specific QC checks based on available group types
        group_types = {g.group_type for g in manifest.file_groups}

        if FileGroupType.RAW_UMI_MATRIX in group_types:
            checks.append("raw_umi_matrix: verify integer counts (no floats, no log-transform)")
            checks.append("raw_umi_matrix: verify gene identifiers in first column")
            checks.append("raw_umi_matrix: preserve raw counts in layers['counts']")

        if FileGroupType.CSV_MATRIX in group_types:
            checks.append("csv_matrix: verify matrix shape (genes x cells)")
            checks.append("csv_matrix: verify delimiter consistency (comma/tab)")
            checks.append("csv_matrix: verify gene identifiers in first column")
            checks.append("csv_matrix: flag for manual review if counts appear normalized (non-integer)")

        if FileGroupType.NORMALIZED_MATRIX in group_types:
            checks.append("normalized_matrix: mark output as limited-use (no raw counts available)")
            checks.append("normalized_matrix: preserve in layers['normalized'] — do NOT copy to layers['counts']")
            checks.append("normalized_matrix: verify normalization method if detectable from value distributions")

        return checks
