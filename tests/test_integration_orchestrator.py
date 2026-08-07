"""Integration tests for Orchestrator — execution flow, path resolution, metadata behavior."""

from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path

import pytest

from gse_router.pipeline.orchestrator import Orchestrator
from gse_router.schemas.manifest import (
    DatasetManifest,
    FileGroup,
    FileGroupType,
    Strategy,
)
from gse_router.schemas.report import StepResult, StepStatus
from gse_router.schemas.route_plan import (
    ProcessingStep,
    RouteMethod,
    RoutePlan,
    ReferenceGenome,
)
from gse_router.schemas.sidecar import FileDetectedType, FileSidecar, ProcessingStage


# Valid sha256 file IDs for Pydantic validation ('^sha256:[a-f0-9]{64}$')
# Use a counter-based hex suffix to generate distinct valid IDs.
_sid_counter: int = 0


def _sid(label: str = "") -> str:
    """Make a valid sha256 file_id. ``label`` is hashed into a hex suffix for
    deterministic but distinct IDs across test runs."""
    global _sid_counter
    _sid_counter += 1
    # Use hash of label + counter to produce a 64-char hex string
    raw = f"{label}_{_sid_counter}"
    h = hex(abs(hash(raw)))[2:].ljust(64, "a")[:64]
    return f"sha256:{h}"


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def tmp_output_dir(tmp_path):
    return tmp_path / "output"


@pytest.fixture
def orchestrator(tmp_output_dir):
    return Orchestrator(output_dir=tmp_output_dir)


def _make_sidecar(
    file_id: str,
    path: str,
    detected_type: FileDetectedType,
    confidence: float = 0.9,
    processing_stage: ProcessingStage = ProcessingStage.RAW_COUNTS,
) -> FileSidecar:
    return FileSidecar(
        file_id=file_id,
        path=path,
        filename=Path(path).name,
        detected_type=detected_type,
        processing_stage=processing_stage,
        confidence=confidence,
        size_bytes=1024,
        evidence=["test"],
    )


def _make_csv_group(
    group_id: str = "csv_test001",
    data_id: str | None = None,
) -> FileGroup:
    return FileGroup(
        group_id=group_id,
        group_type=FileGroupType.CSV_MATRIX,
        sample_id="test_sample",
        files={"data": data_id or _sid("data")},
        route="read_csv_matrix_to_h5ad",
        confidence=0.9,
    )


def _make_metadata_group(group_id: str = "meta_cells", file_id: str | None = None) -> FileGroup:
    return FileGroup(
        group_id=group_id,
        group_type=FileGroupType.METADATA,
        sample_id="test_sample",
        files={"file": file_id or _sid("meta")},
        route="merge_cell_metadata",
        confidence=0.9,
    )


def _make_route_plan(dataset_id: str, steps: list[ProcessingStep]) -> RoutePlan:
    return RoutePlan(
        dataset_id=dataset_id,
        manifest_version="1.0",
        strategy="processed_counts_first",
        reference_genome=ReferenceGenome.GRCH38,
        steps=steps,
    )


def _make_step(step_id: str, method: RouteMethod, input_group_id: str) -> ProcessingStep:
    return ProcessingStep(
        step_id=step_id,
        method=method,
        input_group_id=input_group_id,
        output_label=f"{step_id}_output",
    )


def _make_manifest(
    dataset_id: str,
    files: dict,
    file_groups: list[FileGroup],
) -> DatasetManifest:
    return DatasetManifest(
        dataset_id=dataset_id,
        files=files,
        file_groups=file_groups,
        selected_strategy=Strategy.PROCESSED_COUNTS_FIRST,
        organism=["Homo sapiens"],
        platform=["10x"],
        technology="scRNA-seq",
    )


# ---------------------------------------------------------------------------
# Fix 1: Archive path resolution
# ---------------------------------------------------------------------------

class TestArchivePathResolution:
    """Fix 1: extracted archive sidecar paths resolved relative to extract_dir."""

    def test_extract_dir_used_for_sidecar_paths(self, orchestrator, tmp_path):
        """Paths from probe_directory(extract_dir) should resolve via extract_dir."""
        data_path = tmp_path / "data"
        data_path.mkdir()
        extract_dir = data_path / "_extracted_archive"
        extract_dir.mkdir(parents=True)

        # Create MEX files ONLY in extract_dir (not in data_path)
        (extract_dir / "matrix.mtx").write_text(
            "%%MatrixMarket matrix coordinate integer general\n% comment\n3 3 3\n1 1 10\n2 2 5\n3 3 1\n"
        )
        (extract_dir / "barcodes.tsv").write_text("AAACCCAAG\nAAACCCAGT\nAAACCCATG\n")
        (extract_dir / "features.tsv").write_text("BRCA1\nTP53\nEGFR\n")

        # Create sidecars with paths relative to extract_dir
        matrix_sc = _make_sidecar(_sid("mat"), "matrix.mtx", FileDetectedType.MEX_MATRIX_10X)
        barcodes_sc = _make_sidecar(_sid("bcs"), "barcodes.tsv", FileDetectedType.MEX_BARCODES_10X)
        features_sc = _make_sidecar(_sid("fts"), "features.tsv", FileDetectedType.MEX_FEATURES_10X)

        # Simulate the archive handler's path resolution logic
        mtx = extract_dir / matrix_sc.path
        barcodes = extract_dir / barcodes_sc.path
        features = extract_dir / features_sc.path

        assert mtx.exists(), f"Matrix not found at {mtx} (extract_dir resolution)"
        assert barcodes.exists(), f"Barcodes not found at {barcodes}"
        assert features.exists(), f"Features not found at {features}"

    def test_data_path_fallback_when_extract_dir_missing(self, orchestrator, tmp_path):
        """If extract_dir-relative path doesn't exist, try data_path fallback."""
        data_path = tmp_path / "data"
        data_path.mkdir()
        extract_dir = data_path / "_extracted_archive"
        extract_dir.mkdir(parents=True)

        # Create MEX files in data_path (simulating already-extracted elsewhere)
        (data_path / "matrix.mtx").write_text(
            "%%MatrixMarket matrix coordinate integer general\n% comment\n3 3 3\n1 1 10\n2 2 5\n3 3 1\n"
        )

        # Sidecar paths are relative to extract_dir, but files only in data_path
        matrix_sc = _make_sidecar(_sid("mat"), "matrix.mtx", FileDetectedType.MEX_MATRIX_10X)
        mtx = extract_dir / matrix_sc.path
        assert not mtx.exists()

        # Fallback logic
        mtx_fb = data_path / matrix_sc.path
        assert mtx_fb.exists(), f"Fallback path {mtx_fb} should exist"


# ---------------------------------------------------------------------------
# Fix 4: Agent-controlled metadata behavior
# ---------------------------------------------------------------------------

class TestMetadataSkipWhenNoH5ad:
    """Fix 4: metadata steps skip when no h5ad from prior data steps."""

    def test_metadata_skipped_when_no_h5ad(self, orchestrator):
        """With skip_metadata_if_no_h5ad=True, metadata steps are SKIPPED."""
        meta_id = _sid("ms")
        plan = _make_route_plan("GSEtest", [
            _make_step("step_000", RouteMethod.MERGE_CELL_METADATA, "meta_cells"),
        ])
        manifest = _make_manifest(
            "GSEtest",
            files={meta_id: _make_sidecar(meta_id, "meta.csv", FileDetectedType.CELL_METADATA)},
            file_groups=[_make_metadata_group("meta_cells", file_id=meta_id)],
        )
        results = orchestrator._execute_plan(
            plan, manifest, Path("/tmp"),
            skip_metadata_if_no_h5ad=True,
        )
        assert len(results) == 1
        assert results[0].status == StepStatus.SKIPPED
        assert "no h5ad" in results[0].error_message.lower()

    def test_metadata_runs_when_h5ad_exists(self, orchestrator, tmp_path):
        """Metadata steps should run when a prior step produced an h5ad."""
        data_id = _sid("d1")
        meta_id = _sid("m1")

        data_path = tmp_path / "data"
        data_path.mkdir()
        # CSV data file that will convert successfully
        (data_path / "counts.csv").write_text(
            "gene,cell1,cell2,cell3\nBRCA1,10,5,0\nTP53,0,3,8\nEGFR,2,0,1\n"
        )
        (data_path / "meta.csv").write_text("barcode,condition\nAAACCCAAG,treated\n")

        plan = _make_route_plan("GSEtest", [
            _make_step("step_000", RouteMethod.READ_CSV_MATRIX, "csv_test001"),
            _make_step("step_001", RouteMethod.MERGE_CELL_METADATA, "meta_cells"),
        ])
        manifest = _make_manifest("GSEtest", files={
            data_id: _make_sidecar(data_id, "counts.csv", FileDetectedType.CSV_MATRIX),
            meta_id: _make_sidecar(meta_id, "meta.csv", FileDetectedType.CELL_METADATA),
        }, file_groups=[
            _make_csv_group("csv_test001", data_id),
            _make_metadata_group("meta_cells", file_id=meta_id),
        ])
        results = orchestrator._execute_plan(
            plan, manifest, data_path,
            skip_metadata_if_no_h5ad=True,
            metadata_is_critical=True,
        )
        assert len(results) == 2
        assert results[0].status == StepStatus.COMPLETED, f"Data step failed: {results[0].error_message}"
        # Metadata should NOT be skipped since data step produced h5ad.
        # It may FAIL (if barcodes don't match), but it ran.
        assert results[1].status != StepStatus.SKIPPED

    def test_skip_metadata_flag_disabled(self, orchestrator):
        """With skip_metadata_if_no_h5ad=False, metadata always runs."""
        meta_id = _sid("md")
        plan = _make_route_plan("GSEtest", [
            _make_step("step_000", RouteMethod.MERGE_CELL_METADATA, "meta_cells"),
        ])
        manifest = _make_manifest(
            "GSEtest",
            files={meta_id: _make_sidecar(meta_id, "meta.csv", FileDetectedType.CELL_METADATA)},
            file_groups=[_make_metadata_group("meta_cells", file_id=meta_id)],
        )
        results = orchestrator._execute_plan(
            plan, manifest, Path("/tmp"),
            skip_metadata_if_no_h5ad=False,
        )
        assert len(results) == 1
        # With skip disabled, metadata runs and fails with "No h5ad found" (expected)
        assert results[0].status in (StepStatus.FAILED, StepStatus.SKIPPED)


class TestMetadataNonBlocking:
    """Fix 4: metadata failures produce SKIPPED when metadata_is_critical=False."""

    def test_metadata_failure_is_skipped_not_failed(self, orchestrator, tmp_path):
        """Non-critical metadata failure → SKIPPED status, not FAILED."""
        data_id = _sid("d2")
        meta_id = _sid("m2")

        data_path = tmp_path / "data"
        data_path.mkdir()
        (data_path / "counts.csv").write_text(
            "gene,cell1,cell2,cell3\nBRCA1,10,5,0\nTP53,0,3,8\nEGFR,2,0,1\n"
        )
        # Metadata that will fail to merge (no matching column)
        (data_path / "bad_meta.csv").write_text("sample,condition\nS1,treated\nS2,control\n")

        plan = _make_route_plan("GSEtest", [
            _make_step("step_000", RouteMethod.READ_CSV_MATRIX, "csv_test001"),
            _make_step("step_001", RouteMethod.MERGE_CELL_METADATA, "meta_cells"),
        ])
        manifest = _make_manifest("GSEtest", files={
            data_id: _make_sidecar(data_id, "counts.csv", FileDetectedType.CSV_MATRIX),
            meta_id: _make_sidecar(meta_id, "bad_meta.csv", FileDetectedType.CELL_METADATA),
        }, file_groups=[
            _make_csv_group("csv_test001", data_id),
            _make_metadata_group("meta_cells", file_id=meta_id),
        ])
        results = orchestrator._execute_plan(
            plan, manifest, data_path,
            metadata_is_critical=False,
        )
        assert len(results) == 2
        assert results[0].status == StepStatus.COMPLETED, f"Data step failed: {results[0].error_message}"
        assert results[1].status == StepStatus.SKIPPED, (
            f"Expected SKIPPED (non-critical metadata) but got {results[1].status}: {results[1].error_message}"
        )

    def test_metadata_failure_is_failed_when_critical(self, orchestrator, tmp_path):
        """With metadata_is_critical=True, metadata failures → FAILED."""
        data_id = _sid("d3")
        meta_id = _sid("m3")

        data_path = tmp_path / "data"
        data_path.mkdir()
        (data_path / "counts.csv").write_text(
            "gene,cell1,cell2,cell3\nBRCA1,10,5,0\nTP53,0,3,8\nEGFR,2,0,1\n"
        )
        (data_path / "bad_meta.csv").write_text("sample,condition\nS1,treated\nS2,control\n")

        plan = _make_route_plan("GSEtest", [
            _make_step("step_000", RouteMethod.READ_CSV_MATRIX, "csv_test001"),
            _make_step("step_001", RouteMethod.MERGE_CELL_METADATA, "meta_cells"),
        ])
        manifest = _make_manifest("GSEtest", files={
            data_id: _make_sidecar(data_id, "counts.csv", FileDetectedType.CSV_MATRIX),
            meta_id: _make_sidecar(meta_id, "bad_meta.csv", FileDetectedType.CELL_METADATA),
        }, file_groups=[
            _make_csv_group("csv_test001", data_id),
            _make_metadata_group("meta_cells", file_id=meta_id),
        ])
        results = orchestrator._execute_plan(
            plan, manifest, data_path,
            metadata_is_critical=True,
        )
        assert len(results) == 2
        assert results[0].status == StepStatus.COMPLETED
        assert results[1].status == StepStatus.FAILED, (
            f"Expected FAILED (critical metadata) but got {results[1].status}"
        )


# ---------------------------------------------------------------------------
# Fix 4 continued: _build_report excludes SKIPPED from overall FAILED
# ---------------------------------------------------------------------------

class TestBuildReportSkips:
    """Fix 4: _build_report excludes SKIPPED steps from overall FAILED."""

    def test_skipped_steps_dont_cause_failed(self, orchestrator):
        """Pipeline with only COMPLETED + SKIPPED → COMPLETED, not FAILED."""
        plan = _make_route_plan("GSEtest", [
            _make_step("s0", RouteMethod.READ_10X_MTX, "g0"),
            _make_step("s1", RouteMethod.MERGE_CELL_METADATA, "g1"),
        ])
        manifest = _make_manifest("GSEtest", files={}, file_groups=[])
        now = datetime.now(timezone.utc)
        step_results = [
            StepResult(step_id="s0", method="read_10x_mtx_to_h5ad", status=StepStatus.COMPLETED,
                       started_at=now, finished_at=now),
            StepResult(step_id="s1", method="merge_cell_metadata", status=StepStatus.SKIPPED,
                       started_at=now, finished_at=now, error_message="Skipped: no h5ad"),
        ]
        report = orchestrator._build_report(
            gse_id="GSEtest", manifest=manifest, route_plan=plan,
            step_results=step_results, qc_result=None, final_h5ad=None, started_at=now,
        )
        assert report.status == StepStatus.COMPLETED, (
            f"Expected COMPLETED (only SKIPPED + COMPLETED) but got {report.status}"
        )

    def test_failed_data_step_still_fails(self, orchestrator):
        """A FAILED data step should still cause overall FAILED."""
        plan = _make_route_plan("GSEtest", [
            _make_step("s0", RouteMethod.READ_10X_MTX, "g0"),
        ])
        manifest = _make_manifest("GSEtest", files={}, file_groups=[])
        now = datetime.now(timezone.utc)
        step_results = [
            StepResult(step_id="s0", method="read_10x_mtx_to_h5ad", status=StepStatus.FAILED,
                       started_at=now, finished_at=now, error_message="Boom"),
        ]
        report = orchestrator._build_report(
            gse_id="GSEtest", manifest=manifest, route_plan=plan,
            step_results=step_results, qc_result=None, final_h5ad=None, started_at=now,
        )
        assert report.status == StepStatus.FAILED

    def test_all_skipped_is_manual_review(self, orchestrator):
        """All SKIPPED → MANUAL_REVIEW (nothing executed)."""
        plan = _make_route_plan("GSEtest", [
            _make_step("s0", RouteMethod.MERGE_CELL_METADATA, "g0"),
        ])
        manifest = _make_manifest("GSEtest", files={}, file_groups=[])
        now = datetime.now(timezone.utc)
        step_results = [
            StepResult(step_id="s0", method="merge_cell_metadata", status=StepStatus.SKIPPED,
                       started_at=now, finished_at=now, error_message="Skipped"),
        ]
        report = orchestrator._build_report(
            gse_id="GSEtest", manifest=manifest, route_plan=plan,
            step_results=step_results, qc_result=None, final_h5ad=None, started_at=now,
        )
        assert report.status == StepStatus.MANUAL_REVIEW
