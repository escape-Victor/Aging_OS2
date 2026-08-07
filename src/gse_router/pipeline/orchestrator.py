"""Pipeline orchestrator — ties probes, router, and tools into one workflow.

The orchestrator runs the full loop:
  probe -> manifest -> route plan -> execute -> validate -> report
"""

from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path

import structlog

from ..probes.file_probe import FileProbe
from ..probes.metadata import MetadataProbe, StudyContext
from ..router.classifier import FileClassifier
from ..router.route_planner import RoutePlanner
from ..schemas.manifest import DatasetManifest, FileGroup, FileGroupType, Strategy
from ..schemas.sidecar import FileDetectedType
from ..schemas.report import (
    AnnDataSummary,
    ProcessingReport,
    QCMetric,
    StepResult,
    StepStatus,
)
from ..schemas.route_plan import RouteMethod, RoutePlan
from ..tools.converters import ConversionResult, FormatConverters
from ..tools.metadata_tools import MetadataTools
from ..tools.validators import AnnDataValidator, QCResult

logger = structlog.get_logger(__name__)


class Orchestrator:
    """Orchestrates the full probe -> plan -> execute -> validate pipeline.

    This is the top-level coordinator that:
    1. Runs metadata and file probes
    2. Builds the dataset manifest
    3. Plans the processing route
    4. Executes each step via the appropriate tool
    5. Validates outputs
    6. Produces the final report
    """

    # Fallback route chain: when a primary conversion method fails,
    # try the listed alternatives in order before giving up.
    # Ported from v1's consistent graceful-degradation pattern.
    FALLBACK_ROUTES: dict[RouteMethod, list[RouteMethod]] = {
        RouteMethod.READ_10X_MTX: [RouteMethod.READ_CSV_MATRIX],
        RouteMethod.READ_10X_H5: [RouteMethod.READ_H5AD],
        RouteMethod.READ_LOOM: [RouteMethod.READ_CSV_MATRIX],
        RouteMethod.READ_RDS_SEURAT: [RouteMethod.READ_RDS_SCE],
        RouteMethod.READ_RDS_SCE: [RouteMethod.READ_RDS_SEURAT],
        RouteMethod.READ_RAW_UMI_MATRIX: [RouteMethod.READ_CSV_MATRIX],
    }

    def __init__(
        self,
        output_dir: Path | None = None,
        data_dir: Path | None = None,
        email: str | None = None,
        api_key: str | None = None,
        download_backend: str = "auto",
        aspera_key_path: str | None = None,
    ) -> None:
        self.output_dir = output_dir or Path.cwd() / "output"
        self.data_dir = data_dir or Path.cwd() / "data"
        self.output_dir.mkdir(parents=True, exist_ok=True)

        self.download_backend = download_backend
        self.aspera_key_path = aspera_key_path

        self.metadata_probe = MetadataProbe(email=email, api_key=api_key)
        self.file_probe = FileProbe()
        self.classifier = FileClassifier()
        self.route_planner = RoutePlanner()
        self.converters = FormatConverters(output_dir=self.output_dir / "h5ad")
        self.metadata_tools = MetadataTools()
        self.validator = AnnDataValidator()

        self._logger = logger.bind(component="Orchestrator")

    # ------------------------------------------------------------------
    # Main pipeline
    # ------------------------------------------------------------------

    async def process(
        self,
        gse_id: str,
        data_path: Path | None = None,
        skip_metadata_probe: bool = False,
        skip_raw_fastq: bool = True,
        dry_run: bool = False,
        auto_download: bool = False,
        # Agent-controlled policy parameters
        skip_metadata_if_no_h5ad: bool = True,
        metadata_is_critical: bool = False,
    ) -> ProcessingReport:
        """Run the full pipeline for a GSE accession.

        Args:
            gse_id: GEO Series accession (e.g., 'GSE123456').
            data_path: Path to downloaded data. If None, expects data in self.data_dir/gse_id/.
            skip_metadata_probe: Skip GEO/SRA metadata queries.
            skip_raw_fastq: Skip raw FASTQ processing, prefer processed counts.
            dry_run: Plan only, no execution.
            auto_download: Download GEO supplementary files before the file probe.
            skip_metadata_if_no_h5ad: If True (default), metadata/annotation merge
                steps are SKIPPED when no h5ad was produced by prior data steps.
                Set to False to run metadata steps regardless.
            metadata_is_critical: If True, metadata merge failures cause overall
                FAILED status. If False (default), metadata failures produce
                SKIPPED status and don't block the pipeline.

        Returns:
            ProcessingReport with full provenance.
        """
        started_at = datetime.now(timezone.utc)
        self._logger.info("Starting pipeline", gse_id=gse_id)

        # Phase 1: Metadata probe
        study_context: StudyContext | None = None
        if not skip_metadata_probe:
            # Collect supplementary filenames for Layer 3 technology detection
            # (checks GEO FTP listing for 10x format markers in filenames)
            supplementary_filenames: list[str] | None = None
            try:
                from ..tools.geo_downloader import GeoSupplementaryLister
                lister = GeoSupplementaryLister()
                files_listing = await lister.list_supplementary(gse_id)
                if files_listing:
                    supplementary_filenames = [f.filename for f in files_listing if f.filename]
            except Exception:
                pass  # Best-effort; metadata probe works without filenames

            study_context = await self.metadata_probe.probe_all(
                gse_id, supplementary_filenames=supplementary_filenames
            )
        if auto_download:
            try:
                from ..tools.geo_downloader import GeoFileDownloader
                from ..tools.download_config import load_download_config, DownloadBackend

                dl_config = load_download_config()
                dl_config.backend = DownloadBackend(self.download_backend)
                if self.aspera_key_path:
                    dl_config.aspera_key_path = self.aspera_key_path

                downloader = GeoFileDownloader(
                    data_dir=data_path.parent if data_path else self.data_dir,
                    use_accelerator=True,
                    download_config=dl_config,
                )

                # Download ALL data files — don't pre-filter by type.
                # The file probe will inspect everything and the route planner
                # will pick the best processing route. Archives are downloaded
                # and auto-unpacked so their contents are visible to the probe.
                self._logger.info("Auto-downloading supplementary files", gse_id=gse_id)
                dataset_dir, _all_files, _downloaded = await downloader.download_dataset(
                    gse_id=gse_id,
                    file_types=None,       # greedy: download all data categories
                    max_size_mb=dl_config.max_download_size_mb,
                    skip_raw=True,         # still skip FASTQ/BAM
                )
                data_path = dataset_dir  # Unpack from (dir, files, downloads) tuple
                self._logger.info(
                    "Download complete",
                    gse_id=gse_id,
                    path=str(data_path),
                )
            except Exception as exc:
                self._logger.warning(
                    "Auto-download failed, falling back to local data",
                    gse_id=gse_id,
                    error=str(exc),
                )
                # Fall back to local data path
                data_path = data_path or self.data_dir / gse_id

        # Phase 2: File probe
        data_path = data_path or self.data_dir / gse_id
        sidecars = self.file_probe.probe_directory(data_path)
        sidecar_dict = {sc.file_id: sc for sc in sidecars}

        # Phase 3: Build manifest
        manifest = self._build_manifest(gse_id, sidecar_dict, study_context)

        # Phase 4: Plan route
        route_plan = self.route_planner.plan(manifest, study_context)

        # Dry run: return report with plan but no execution
        if dry_run:
            report = self._build_dry_run_report(gse_id, manifest, route_plan, started_at)
            report_path = self.output_dir / f"{gse_id}_report.json"
            report_path.write_text(report.model_dump_json(indent=2))
            return report

        # Phase 5: Execute steps
        step_results = self._execute_plan(
            route_plan, manifest, data_path,
            skip_metadata_if_no_h5ad=skip_metadata_if_no_h5ad,
            metadata_is_critical=metadata_is_critical,
        )

        # Phase 6: Validate final output
        final_h5ad = self.output_dir / "h5ad" / f"{gse_id}_final.h5ad"
        qc_result = self.validator.validate(final_h5ad) if final_h5ad.exists() else None

        # Phase 7: Build report
        report = self._build_report(
            gse_id=gse_id,
            manifest=manifest,
            route_plan=route_plan,
            step_results=step_results,
            qc_result=qc_result,
            final_h5ad=final_h5ad,
            started_at=started_at,
        )

        # Write report
        report_path = self.output_dir / f"{gse_id}_report.json"
        report_path.write_text(report.model_dump_json(indent=2))

        self._logger.info("Pipeline complete", gse_id=gse_id, report=str(report_path))
        return report

    # ------------------------------------------------------------------
    # Phase implementations
    # ------------------------------------------------------------------

    def _build_manifest(
        self,
        gse_id: str,
        sidecars: dict,
        study_context: StudyContext | None,
    ) -> DatasetManifest:
        """Build the dataset manifest from file sidecars and study context."""
        file_groups = self.classifier.group_files(sidecars)
        strategy = self.classifier.select_strategy(file_groups, prefer_processed=True)

        organism: list[str] = []
        if study_context:
            organism = study_context.organism

        unclassified = [
            fid for fid, sc in sidecars.items()
            if sc.detected_type.value == "unknown"
        ]

        manifest = DatasetManifest(
            dataset_id=gse_id,
            title=study_context.title if study_context else None,
            organism=organism,
            platform=study_context.platform if study_context else [],
            technology=study_context.technology if study_context else "scRNA-seq",
            tissue_or_organ=study_context.tissue_or_organ if study_context else [],
            disease=study_context.disease if study_context else [],
            condition_groups=study_context.condition_groups if study_context else [],
            publication=study_context.publication if study_context else {},
            metadata_confidence=study_context.confidence if study_context else {},
            files=sidecars,
            file_groups=file_groups,
            selected_strategy=strategy,
            unclassified_files=unclassified,
            manual_review_needed=len(unclassified) > 0 or strategy == Strategy.UNSOLVABLE,
            review_reasons=(
                [f"{len(unclassified)} unclassified files"] if unclassified else []
            ),
            source_accessions=self._collect_accessions(gse_id, study_context),
        )
        return manifest

    def _cross_validate_species(
        self,
        h5ad_path: Path,
        expected_species: str,
    ) -> list[str]:
        """Cross-validate species using gene symbol casing patterns.

        Reads gene names from the output h5ad and compares the gene-inferred
        species against the expected species from GEO metadata. Returns a
        list of warning strings if there is a mismatch.

        This catches GEO metadata errors like GSE149483 (plant scRNA-seq
        mislabeled as human/rat in GEO).

        Args:
            h5ad_path: Path to the output .h5ad file.
            expected_species: Expected organism from GEO metadata.

        Returns:
            List of warning strings, or empty list if no mismatch detected.
        """
        try:
            from ..probes.metadata import MetadataProbe
            import anndata

            adata = anndata.read_h5ad(h5ad_path)
            gene_species = MetadataProbe.detect_species_from_gene_symbols(
                list(adata.var_names)
            )

            if gene_species in ("Mixed_or_Unknown", "Unknown"):
                return []

            expected_lower = expected_species.lower()
            gene_lower = gene_species.lower()

            if "human" in gene_lower and "mouse" in expected_lower:
                return [
                    f"Species cross-validation mismatch: gene symbols look Human "
                    f"(UPPERCASE) but GEO metadata says '{expected_species}'. "
                    f"Verify organism assignment."
                ]
            if "mouse" in gene_lower and "human" in expected_lower:
                return [
                    f"Species cross-validation mismatch: gene symbols look Mouse "
                    f"(Titlecase) but GEO metadata says '{expected_species}'. "
                    f"Verify organism assignment."
                ]
        except Exception:
            pass  # Cross-validation is advisory; never blocks
        return []

    def _execute_plan(
        self,
        plan: RoutePlan,
        manifest: DatasetManifest,
        data_path: Path,
        skip_metadata_if_no_h5ad: bool = True,
        metadata_is_critical: bool = False,
    ) -> list[StepResult]:
        """Execute each step in the route plan.

        Args:
            plan: The route plan with ordered steps.
            manifest: The dataset manifest.
            data_path: Root data directory.
            skip_metadata_if_no_h5ad: If True, skip metadata merge when no prior
                h5ad was produced (avoids cascading "No h5ad found" errors).
            metadata_is_critical: If True, metadata failures propagate as FAILED.
                If False, metadata failures become SKIPPED (non-blocking).
        """
        results: list[StepResult] = []
        method_dispatcher = self._build_dispatcher(data_path, manifest)
        expected_species = ", ".join(manifest.organism) if manifest.organism else "unknown"

        # Track h5ad output paths from successful data conversion steps.
        # Metadata handlers receive this list explicitly rather than globbing
        # the output directory, ensuring they operate on the correct file.
        # When empty, metadata merge steps are skipped gracefully rather than
        # producing cascading "No h5ad found" errors.
        produced_h5ad_paths: list[Path] = []

        # Metadata step methods (merge/annotation — not data conversion)
        _METADATA_METHODS = {RouteMethod.MERGE_CELL_METADATA, RouteMethod.MERGE_GENE_ANNOTATION}

        for step in plan.steps:
            started = datetime.now(timezone.utc)
            self._logger.info("Executing step", step_id=step.step_id, method=step.method.value)

            # Agent-controlled: skip metadata steps when no h5ad exists
            if skip_metadata_if_no_h5ad and not produced_h5ad_paths and step.method in _METADATA_METHODS:
                self._logger.info(
                    "Skipping metadata step — no h5ad from prior steps",
                    step_id=step.step_id,
                    method=step.method.value,
                )
                result = StepResult(
                    step_id=step.step_id,
                    method=step.method.value,
                    status=StepStatus.SKIPPED,
                    started_at=started,
                    finished_at=datetime.now(timezone.utc),
                    error_message="Skipped: no h5ad produced by prior data steps",
                )
                results.append(result)
                continue

            try:
                handler = method_dispatcher.get(step.method)
                if handler is None:
                    result = StepResult(
                        step_id=step.step_id,
                        method=step.method.value,
                        status=StepStatus.SKIPPED,
                        started_at=started,
                        finished_at=datetime.now(timezone.utc),
                        error_message=f"No handler for method: {step.method.value}",
                    )
                else:
                    # Metadata handlers receive explicit h5ad paths from
                    # prior data steps instead of globbing the output dir.
                    if step.method in _METADATA_METHODS:
                        conversion = handler(step, manifest, produced_h5ad_paths)
                    else:
                        conversion = handler(step, manifest)
                    # If primary conversion failed, try fallback routes
                    if not conversion.success:
                        fallbacks = self.FALLBACK_ROUTES.get(step.method, [])
                        for fallback_method in fallbacks:
                            fallback_handler = method_dispatcher.get(fallback_method)
                            if fallback_handler is None:
                                continue
                            self._logger.info(
                                "Trying fallback route",
                                step_id=step.step_id,
                                original=step.method.value,
                                fallback=fallback_method.value,
                            )
                            try:
                                conversion = fallback_handler(step, manifest)
                                if conversion.success:
                                    conversion.warnings.append(
                                        f"Fallback from {step.method.value} to "
                                        f"{fallback_method.value} used"
                                    )
                                    break
                            except Exception:
                                continue
                    # Cross-validate species after successful data conversion steps
                    if conversion.success and conversion.output_path and expected_species != "unknown":
                        species_warnings = self._cross_validate_species(
                            conversion.output_path, expected_species
                        )
                        if species_warnings:
                            conversion.warnings.extend(species_warnings)

                    # Agent-controlled: metadata failures → SKIPPED when non-critical
                    if (
                        not conversion.success
                        and step.method in _METADATA_METHODS
                        and not metadata_is_critical
                    ):
                        result = StepResult(
                            step_id=step.step_id,
                            method=step.method.value,
                            status=StepStatus.SKIPPED,
                            started_at=started,
                            finished_at=datetime.now(timezone.utc),
                            output_artifacts=(
                                [str(conversion.output_path)] if conversion.output_path else []
                            ),
                            error_message=conversion.error,
                            warnings=conversion.warnings,
                        )
                    else:
                        result = self._build_step_result(step, conversion, started)

                    # Track h5ad output paths from successful data steps
                    if conversion.success and conversion.output_path:
                        produced_h5ad_paths.append(conversion.output_path)
            except Exception as exc:
                self._logger.error("Step failed", step_id=step.step_id, error=str(exc))
                # Agent-controlled: metadata exceptions → SKIPPED when non-critical
                if step.method in _METADATA_METHODS and not metadata_is_critical:
                    result = StepResult(
                        step_id=step.step_id,
                        method=step.method.value,
                        status=StepStatus.SKIPPED,
                        started_at=started,
                        finished_at=datetime.now(timezone.utc),
                        error_message=str(exc),
                    )
                else:
                    result = StepResult(
                        step_id=step.step_id,
                        method=step.method.value,
                        status=StepStatus.FAILED,
                        started_at=started,
                        finished_at=datetime.now(timezone.utc),
                        error_message=str(exc),
                    )
            results.append(result)

        return results

    def _build_report(
        self,
        gse_id: str,
        manifest: DatasetManifest,
        route_plan: RoutePlan,
        step_results: list[StepResult],
        qc_result: QCResult | None,
        final_h5ad: Path | None,
        started_at: datetime,
    ) -> ProcessingReport:
        """Build the final processing report."""
        finished_at = datetime.now(timezone.utc)
        duration = (finished_at - started_at).total_seconds()

        # Overall status — SKIPPED steps don't affect the outcome.
        # A pipeline with only SKIPPED metadata steps and successful data
        # conversion should not be marked FAILED.
        non_skipped = [r for r in step_results if r.status != StepStatus.SKIPPED]
        all_complete = all(r.status == StepStatus.COMPLETED for r in non_skipped) if non_skipped else False
        any_failed = any(r.status == StepStatus.FAILED for r in non_skipped)
        if any_failed:
            overall = StepStatus.FAILED
        elif all_complete:
            overall = StepStatus.COMPLETED
        elif not non_skipped:
            overall = StepStatus.MANUAL_REVIEW
        else:
            overall = StepStatus.MANUAL_REVIEW

        # QC metrics
        qc_metrics: list[QCMetric] = []
        if qc_result:
            for name, value in qc_result.metrics.items():
                qc_metrics.append(QCMetric(name=name, value=value))
            for failure in qc_result.failures:
                qc_metrics.append(QCMetric(name="failure", value=failure, passed=False))

        # Output summary
        output_summary = None
        if final_h5ad and final_h5ad.exists():
            output_summary = AnnDataSummary(
                n_obs=0, n_vars=0,  # Will be populated by validator
                file_path=str(final_h5ad),
                file_size_mb=final_h5ad.stat().st_size / 1e6,
            )

        return ProcessingReport(
            dataset_id=gse_id,
            report_id=f"{gse_id}_{started_at.strftime('%Y%m%d_%H%M%S')}",
            status=overall,
            strategy=manifest.selected_strategy.value,
            reference_genome=route_plan.reference_genome.value if route_plan.reference_genome else None,
            steps=step_results,
            qc_metrics=qc_metrics,
            output=output_summary,
            provenance={
                "pipeline_version": "0.1.0",
                "manifest_id": manifest.dataset_id,
                "route_plan_id": route_plan.dataset_id,
            },
            total_duration_seconds=duration,
            warnings=route_plan.warnings,
            source_manifest_id=manifest.dataset_id,
            source_route_plan_id=route_plan.dataset_id,
        )

    def _build_dry_run_report(
        self,
        gse_id: str,
        manifest: DatasetManifest,
        route_plan: RoutePlan,
        started_at: datetime,
    ) -> ProcessingReport:
        """Build a report for a dry run (plan only, no execution)."""
        return ProcessingReport(
            dataset_id=gse_id,
            report_id=f"{gse_id}_{started_at.strftime('%Y%m%d_%H%M%S')}_dryrun",
            status=StepStatus.PENDING,
            strategy=manifest.selected_strategy.value,
            steps=[
                StepResult(
                    step_id=s.step_id,
                    method=s.method.value,
                    status=StepStatus.PENDING,
                )
                for s in route_plan.steps
            ],
            provenance={"dry_run": True, "step_count": len(route_plan.steps)},
            warnings=route_plan.warnings,
        )

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------

    def _build_dispatcher(self, data_path: Path, manifest: DatasetManifest) -> dict:
        """Build a dispatch table: RouteMethod -> handler function.

        Each handler receives (step, manifest) — the manifest is needed for
        file group lookups via _get_group_for_step.
        """
        return {
            RouteMethod.CELL_RANGER_COUNT: lambda step, m: self._handle_fastq_skipped(step, m, data_path, "cellranger_count"),
            RouteMethod.STAR_SOLO: lambda step, m: self._handle_fastq_skipped(step, m, data_path, "starsolo"),
            RouteMethod.KALLISTO_BUSTOOLS: lambda step, m: self._handle_fastq_skipped(step, m, data_path, "kallisto_bustools"),
            RouteMethod.ALEVIN_FRY: lambda step, m: self._handle_fastq_skipped(step, m, data_path, "alevin_fry"),
            RouteMethod.READ_10X_MTX: lambda step, m: self._handle_10x_mtx_or_archive(step, m, data_path),
            RouteMethod.READ_10X_H5: lambda step, m: self._handle_10x_h5(step, m, data_path),
            RouteMethod.READ_H5AD: lambda step, m: self._handle_h5ad(step, m, data_path),
            RouteMethod.READ_LOOM: lambda step, m: self._handle_loom(step, m, data_path),
            RouteMethod.READ_CSV_MATRIX: lambda step, m: self._handle_csv(step, m, data_path),
            RouteMethod.READ_RAW_UMI_MATRIX: lambda step, m: self._handle_raw_umi(step, m, data_path),
            RouteMethod.READ_NORMALIZED_MATRIX: lambda step, m: self._handle_normalized(step, m, data_path),
            RouteMethod.READ_RDS_SEURAT: lambda step, m: self._handle_rds(step, m, data_path),
            RouteMethod.READ_RDS_SCE: lambda step, m: self._handle_rds(step, m, data_path),
            RouteMethod.MERGE_CELL_METADATA: lambda step, m, h5ads: self._handle_merge_cell_meta(step, m, h5ads, data_path),
            RouteMethod.MERGE_GENE_ANNOTATION: lambda step, m, h5ads: self._handle_merge_gene_anno(step, m, h5ads, data_path),
        }

    def _handle_10x_mtx(self, step, manifest: DatasetManifest, data_path: Path) -> ConversionResult:
        group = self._get_group_for_step(step, manifest)
        mtx = self._resolve_path(group, "matrix", manifest, data_path)
        barcodes = self._resolve_path(group, "barcodes", manifest, data_path)
        features = self._resolve_path(group, "features", manifest, data_path)
        return self.converters.convert_10x_mtx_to_h5ad(mtx, barcodes, features, sample_id=group.sample_id)

    def _handle_10x_mtx_or_archive(self, step, manifest: DatasetManifest, data_path: Path) -> ConversionResult:
        """Handle READ_10X_MTX — may come from a regular MEX triplet or an ARCHIVE_10X.

        For ARCHIVE_10X: extract the archive, probe for MEX triplets, then convert.
        For regular MEX: delegate to the standard handler.
        """
        group = self._get_group_for_step(step, manifest)
        is_archive = group.group_type == FileGroupType.ARCHIVE_10X

        if is_archive:
            # Extract the 10x archive and re-probe its contents
            archive_path = self._resolve_path(group, "file", manifest, data_path)
            extract_dir = data_path / f"_extracted_{archive_path.stem}"
            self._logger.info(
                "Extracting 10x archive",
                archive=str(archive_path),
                dest=str(extract_dir),
            )
            try:
                self.file_probe._extract_archive(archive_path, extract_dir)
            except Exception as exc:
                return ConversionResult(
                    success=False,
                    error=f"Failed to extract 10x archive: {exc}",
                )

            # Probe extracted files for MEX triplets
            extracted_sidecars = self.file_probe.probe_directory(extract_dir)
            mex_triplets = FileClassifier()._group_mex_triplets(
                {sc.file_id: sc for sc in extracted_sidecars}, set()
            )
            if not mex_triplets:
                return ConversionResult(
                    success=False,
                    error=(
                        f"No MEX triplets found in extracted 10x archive. "
                        f"Extracted {len(extracted_sidecars)} files."
                    ),
                )

            # Process the first MEX triplet found
            triplet = mex_triplets[0]
            mtx = data_path / manifest.files[triplet.files["matrix"]].path if triplet.files["matrix"] in manifest.files else None
            if mtx is None:
                # The extracted sidecars aren't in the manifest yet — resolve directly
                matrix_sc = next(
                    (sc for sc in extracted_sidecars if sc.detected_type == FileDetectedType.MEX_MATRIX_10X),
                    None,
                )
                barcodes_sc = next(
                    (sc for sc in extracted_sidecars if sc.detected_type == FileDetectedType.MEX_BARCODES_10X),
                    None,
                )
                features_sc = next(
                    (sc for sc in extracted_sidecars if sc.detected_type == FileDetectedType.MEX_FEATURES_10X),
                    None,
                )
                if not (matrix_sc and barcodes_sc and features_sc):
                    return ConversionResult(
                        success=False,
                        error="Could not resolve MEX triplet components from extracted archive",
                    )
                # Paths from probe_directory(extract_dir) are relative to extract_dir.
                # Try extract_dir-relative first; fall back to data_path-relative
                # in case an earlier extraction populated the data_path tree.
                mtx = extract_dir / matrix_sc.path
                barcodes = extract_dir / barcodes_sc.path
                features = extract_dir / features_sc.path
                if not mtx.exists():
                    mtx_fb = data_path / matrix_sc.path
                    if mtx_fb.exists():
                        mtx = mtx_fb
                        barcodes = data_path / barcodes_sc.path
                        features = data_path / features_sc.path
            else:
                barcodes = data_path / manifest.files[triplet.files["barcodes"]].path
                features = data_path / manifest.files[triplet.files["features"]].path

            return self.converters.convert_10x_mtx_to_h5ad(
                mtx, barcodes, features, sample_id=group.sample_id
            )

        # Regular MEX triplet — delegate to standard handler
        return self._handle_10x_mtx(step, manifest, data_path)

    def _handle_10x_h5(self, step, manifest: DatasetManifest, data_path: Path) -> ConversionResult:
        group = self._get_group_for_step(step, manifest)
        h5_path = self._resolve_path(group, "h5", manifest, data_path)
        return self.converters.convert_10x_h5_to_h5ad(h5_path, sample_id=group.sample_id)

    def _handle_h5ad(self, step, manifest: DatasetManifest, data_path: Path) -> ConversionResult:
        group = self._get_group_for_step(step, manifest)
        h5ad_path = self._resolve_path(group, "data", manifest, data_path)
        return self.converters.validate_and_standardize_h5ad(h5ad_path, sample_id=group.sample_id)

    def _handle_loom(self, step, manifest: DatasetManifest, data_path: Path) -> ConversionResult:
        group = self._get_group_for_step(step, manifest)
        loom_path = self._resolve_path(group, "data", manifest, data_path)
        return self.converters.convert_loom_to_h5ad(loom_path, sample_id=group.sample_id)

    def _handle_csv(self, step, manifest: DatasetManifest, data_path: Path) -> ConversionResult:
        group = self._get_group_for_step(step, manifest)
        file_path = self._resolve_path(group, "data", manifest, data_path)

        # Check if this is actually a MATLAB .mat file (mapped via MATRIX_MAT → CSV_MATRIX)
        file_ext = file_path.suffix.lower()
        if file_ext in (".mat",) or (file_ext == ".gz" and file_path.stem.lower().endswith(".mat")):
            return self.converters.convert_mat_to_h5ad(file_path, sample_id=group.sample_id)

        return self.converters.convert_csv_matrix_to_h5ad(file_path, sample_id=group.sample_id)

    def _handle_raw_umi(self, step, manifest: DatasetManifest, data_path: Path) -> ConversionResult:
        group = self._get_group_for_step(step, manifest)
        matrix_path = self._resolve_path(group, "data", manifest, data_path)
        return self.converters.convert_raw_umi_matrix_to_h5ad(matrix_path, sample_id=group.sample_id)

    def _handle_normalized(self, step, manifest: DatasetManifest, data_path: Path) -> ConversionResult:
        group = self._get_group_for_step(step, manifest)
        matrix_path = self._resolve_path(group, "data", manifest, data_path)
        return self.converters.convert_normalized_matrix_to_h5ad(matrix_path, sample_id=group.sample_id)

    def _handle_rds(self, step, manifest: DatasetManifest, data_path: Path) -> ConversionResult:
        group = self._get_group_for_step(step, manifest)
        rds_path = self._resolve_path(group, "data", manifest, data_path)
        return self.converters.convert_rds_to_h5ad(rds_path, sample_id=group.sample_id)

    def _handle_merge_cell_meta(
        self, step, manifest: DatasetManifest,
        h5ad_paths: list[Path], data_path: Path,
    ) -> ConversionResult:
        # Use the most recent h5ad from prior data steps (explicit path tracking).
        # The list is ordered by execution order, so [-1] is the most recent.
        if not h5ad_paths:
            return ConversionResult(
                success=False, error="No h5ad found to merge metadata into",
            )
        group = self._get_group_for_step(step, manifest)
        meta_path = self._resolve_path(group, "file", manifest, data_path)
        result = self.metadata_tools.merge_cell_metadata(h5ad_paths[-1], meta_path)
        return ConversionResult(
            success=result.success,
            output_path=result.output_path,
            error=result.error,
            provenance=result.provenance,
        )

    def _handle_merge_gene_anno(
        self, step, manifest: DatasetManifest,
        h5ad_paths: list[Path], data_path: Path,
    ) -> ConversionResult:
        # Use the most recent h5ad from prior data steps (explicit path tracking).
        if not h5ad_paths:
            return ConversionResult(
                success=False, error="No h5ad found to merge annotation into",
            )
        group = self._get_group_for_step(step, manifest)
        anno_path = self._resolve_path(group, "file", manifest, data_path)
        result = self.metadata_tools.merge_gene_annotation(h5ad_paths[-1], anno_path)
        return ConversionResult(
            success=result.success,
            output_path=result.output_path,
            error=result.error,
            provenance=result.provenance,
        )

    def _handle_fastq_skipped(
        self, step, manifest: DatasetManifest, data_path: Path, method_name: str
    ) -> ConversionResult:
        """Handle FASTQ/SRA routes by returning a clear skip result.

        Per the cardinal rule: raw FASTQ/SRA processing is disabled by default.
        If a dataset reaches this handler, it means no processed data was found
        and the user must explicitly opt in to raw pipeline processing.
        """
        return ConversionResult(
            success=False,
            error=(
                f"Raw FASTQ/SRA processing ({method_name}) is disabled by default. "
                "Per the cardinal rule: prefer processed counts over raw FASTQ. "
                "To enable: set skip_raw_fastq=False and install the raw-pipeline "
                "environment (see scripts/setup_raw_pipeline.sh). "
                "This requires WSL2 on Windows with Nextflow, Cell Ranger/STARsolo, "
                "and SRA Toolkit installed."
            ),
            warnings=[
                "FASTQ/SRA processing skipped — processed data should be used instead",
                "If no processed data exists in this dataset, manual review is required",
            ],
        )

    @staticmethod
    def _get_group_for_step(step, manifest: DatasetManifest) -> FileGroup:
        """Retrieve the FileGroup for a step from the manifest.

        Looks up the step's input_group_id in the manifest's file_groups list.
        """
        for group in manifest.file_groups:
            if group.group_id == step.input_group_id:
                return group
        raise KeyError(
            f"FileGroup '{step.input_group_id}' not found in manifest for step '{step.step_id}'. "
            f"Available groups: {[g.group_id for g in manifest.file_groups]}"
        )

    @staticmethod
    def _resolve_path(
        group: FileGroup, role: str, manifest: DatasetManifest, data_path: Path
    ) -> Path:
        """Resolve a file role in a group to an absolute filesystem path.

        FileGroup.files maps role -> file_id (sha256:...). This looks up the
        file_id in the manifest to get the actual relative path, then joins
        it with data_path.
        """
        file_id = group.files.get(role, "")
        if not file_id:
            raise KeyError(
                f"Role '{role}' not found in group '{group.group_id}'. "
                f"Available roles: {list(group.files.keys())}"
            )
        sidecar = manifest.files.get(file_id)
        if sidecar is None:
            raise KeyError(
                f"File ID '{file_id}' (role '{role}' in group '{group.group_id}') "
                f"not found in manifest.files"
            )
        return data_path / sidecar.path

    @staticmethod
    def _build_step_result(
        step, conversion: ConversionResult, started_at: datetime
    ) -> StepResult:
        return StepResult(
            step_id=step.step_id,
            method=step.method.value,
            status=StepStatus.COMPLETED if conversion.success else StepStatus.FAILED,
            started_at=started_at,
            finished_at=datetime.now(timezone.utc),
            output_artifacts=[str(conversion.output_path)] if conversion.output_path else [],
            error_message=conversion.error,
            warnings=conversion.warnings,
        )

    @staticmethod
    def _collect_accessions(gse_id: str, study_context: StudyContext | None) -> list[str]:
        accessions = [gse_id]
        if study_context:
            for sample in study_context.samples:
                accessions.append(sample.gsm)
                accessions.extend(sample.sra_runs)
        return accessions
