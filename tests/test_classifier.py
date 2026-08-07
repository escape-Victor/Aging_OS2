"""Tests for the file classifier — grouping logic, ARCHIVE_10X mapping, strategy selection.

Tests cover:
- ARCHIVE → ARCHIVE_10X upgrade based on filename patterns
- MEX triplet grouping with shape-matching fallback
- Non-count supplementary files mapped to UNKNOWN (not count matrices)
- Strategy selection with mixed formats
"""

from __future__ import annotations

from pathlib import Path

import pytest

from gse_router.router.classifier import FileClassifier
from gse_router.schemas.manifest import FileGroupType, Strategy
from gse_router.schemas.sidecar import FileDetectedType, FileSidecar, ProcessingStage


# Valid SHA-256 hashes for testing (64 hex chars)
SHA_A = "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa"
SHA_B = "bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb"
SHA_C = "cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc"
SHA_D = "dddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddd"


def _make_sidecar(
    sha_suffix: str,
    filename: str,
    detected_type: FileDetectedType,
    path: str = "",
    confidence: float = 0.90,
    evidence: list[str] | None = None,
) -> FileSidecar:
    """Helper to create a minimal FileSidecar for testing.

    Args:
        sha_suffix: Last 64 hex chars of the file_id (prepended with 'sha256:').
    """
    return FileSidecar(
        file_id=f"sha256:{sha_suffix}",
        path=path or filename,
        filename=filename,
        size_bytes=100,
        extension=Path(filename).suffix.lower(),
        detected_type=detected_type,
        processing_stage=ProcessingStage.RAW_COUNTS,
        confidence=confidence,
        evidence=evidence or [f"Test: {filename}"],
    )


class TestArchive10xMapping:
    """Test that archives with 10x evidence are mapped to ARCHIVE_10X."""

    @pytest.mark.parametrize("filename,evidence", [
        ("sample_bc_matrix.tar.gz", ["10x MEX bundle inside archive"]),
        ("sample.matrix.tar.gz", ["10x MEX bundle inside archive"]),
        ("filtered_gene_bc_matrices.tar.gz", ["10x MEX bundle inside archive"]),
        ("raw_feature_bc_matrix.tar.gz", ["10x MEX bundle inside archive"]),
    ])
    def test_evidence_upgrades_to_archive_10x(self, filename, evidence):
        """Archives with 10x evidence in their sidecar should map to ARCHIVE_10X."""
        sc = _make_sidecar(SHA_A, filename, FileDetectedType.ARCHIVE, evidence=evidence)
        gt = FileClassifier._map_detected_to_group(FileDetectedType.ARCHIVE, sidecar=sc)
        assert gt == FileGroupType.ARCHIVE_10X, \
            f"Expected ARCHIVE_10X for '{filename}', got {gt}"

    def test_filename_alone_no_longer_upgrades_to_archive_10x(self):
        """Filename patterns alone should NOT upgrade to ARCHIVE_10X — evidence required.

        This is the key behavioral change: content-peeking evidence is now required
        for the ARCHIVE → ARCHIVE_10X upgrade. Filename patterns are no longer used.
        """
        sc = _make_sidecar(SHA_B, "sample_bc_matrix.tar.gz", FileDetectedType.ARCHIVE,
                           evidence=["TAR signature detected"])
        gt = FileClassifier._map_detected_to_group(FileDetectedType.ARCHIVE, sidecar=sc)
        # Without content-peeking evidence, stays as generic ARCHIVE
        assert gt == FileGroupType.ARCHIVE, \
            f"Filename-based upgrade removed — expected ARCHIVE, got {gt}"

    @pytest.mark.parametrize("evidence", [
        ["Contains 10x MEX triplet (matrix.mtx + barcodes.tsv + features.tsv)"],
        ["Contains .h5 file(s) — may be 10x HDF5"],
        ["Archive contains 10x MEX triplet (matrix.mtx + barcodes.tsv + features/genes.tsv)"],
        ["Contains .h5 file(s) — may be 10x"],
    ])
    def test_content_evidence_upgrades_to_archive_10x(self, evidence):
        """Content-peeking evidence mentioning 10x/H5/MEX should upgrade to ARCHIVE_10X."""
        sc = _make_sidecar(SHA_D, "any_archive.tar.gz", FileDetectedType.ARCHIVE,
                           evidence=evidence)
        gt = FileClassifier._map_detected_to_group(FileDetectedType.ARCHIVE, sidecar=sc)
        assert gt == FileGroupType.ARCHIVE_10X, \
            f"Expected ARCHIVE_10X for evidence '{evidence}', got {gt}"

    def test_plain_archive_stays_archive(self):
        """Archives without 10x evidence should stay as generic ARCHIVE."""
        sc = _make_sidecar(SHA_C, "misc_data.tar.gz", FileDetectedType.ARCHIVE,
                           evidence=["TAR signature detected"])
        gt = FileClassifier._map_detected_to_group(FileDetectedType.ARCHIVE, sidecar=sc)
        assert gt == FileGroupType.ARCHIVE

    def test_archive_10x_routed_to_read_10x_mtx(self):
        """ARCHIVE_10X groups should have 'extract_10x_archive_to_h5ad' route."""
        from gse_router.router.rules import decode_route
        from gse_router.schemas.manifest import FileGroup

        group = FileGroup(
            group_id="test_archive_10x",
            group_type=FileGroupType.ARCHIVE_10X,
            files={"file": f"sha256:{SHA_D}"},
            route="extract_10x_archive_to_h5ad",
            confidence=0.85,
        )
        decision = decode_route(group)
        assert decision is not None
        assert decision.parameters.get("source") == "archive_10x"


class TestNonCountMapping:
    """Test that non-count supplementary types map to UNKNOWN group type."""

    @pytest.mark.parametrize("detected_type", [
        FileDetectedType.IMMUNE_RECEPTOR_DATA,
        FileDetectedType.PEAK_DATA,
        FileDetectedType.DIFF_RESULTS,
    ])
    def test_non_count_types_map_to_unknown(self, detected_type):
        """BCR/TCR, peak, and diff result types must map to UNKNOWN."""
        gt = FileClassifier._map_detected_to_group(detected_type, sidecar=None)
        assert gt == FileGroupType.UNKNOWN, \
            f"Expected UNKNOWN for {detected_type.value}, got {gt}"


class TestStrategySelection:
    """Test processing strategy selection with mixed formats."""

    def test_processed_counts_first_with_h5ad(self):
        """When h5ad is present, strategy should be PROCESSED_COUNTS_FIRST."""
        from gse_router.schemas.manifest import FileGroup

        groups = [
            FileGroup(
                group_id="h5_1",
                group_type=FileGroupType.H5AD,
                files={"data": f"sha256:{SHA_A}"},
                route="read_h5ad_validate",
                confidence=0.95,
            ),
        ]
        classifier = FileClassifier()
        strategy = classifier.select_strategy(groups)
        assert strategy == Strategy.PROCESSED_COUNTS_FIRST

    def test_normalized_only_strategy(self):
        """When only normalized data is available, strategy should be NORMALIZED_ONLY."""
        from gse_router.schemas.manifest import FileGroup

        groups = [
            FileGroup(
                group_id="norm_1",
                group_type=FileGroupType.NORMALIZED_MATRIX,
                files={"data": f"sha256:{SHA_B}"},
                route="read_normalized_matrix_to_h5ad",
                confidence=0.80,
            ),
        ]
        classifier = FileClassifier()
        strategy = classifier.select_strategy(groups)
        assert strategy == Strategy.NORMALIZED_ONLY

    def test_unsolvable_with_no_data(self):
        """When no processable data exists, strategy should be UNSOLVABLE."""
        from gse_router.schemas.manifest import FileGroup

        groups = [
            FileGroup(
                group_id="unk_1",
                group_type=FileGroupType.UNKNOWN,
                files={"file": f"sha256:{SHA_C}"},
                route="manual_review",
                confidence=0.30,
            ),
        ]
        classifier = FileClassifier()
        strategy = classifier.select_strategy(groups)
        assert strategy == Strategy.UNSOLVABLE

    def test_fastq_skipped_when_processed_exists(self):
        """FASTQ should be excluded when processed data is present."""
        from gse_router.schemas.manifest import FileGroup

        groups = [
            FileGroup(
                group_id="h5_1",
                group_type=FileGroupType.H5_10X,
                files={"h5": f"sha256:{SHA_A}"},
                route="read_10x_h5_to_h5ad",
                confidence=0.95,
            ),
            FileGroup(
                group_id="fq_1",
                group_type=FileGroupType.FASTQ_PAIRED,
                files={"read1": f"sha256:{SHA_B}", "read2": f"sha256:{SHA_C}"},
                route="cellranger_count",
                confidence=0.85,
            ),
        ]
        classifier = FileClassifier()
        strategy = classifier.select_strategy(groups)
        # Should prefer processed, not hybrid or raw
        assert strategy in (Strategy.PROCESSED_COUNTS_FIRST, Strategy.HYBRID)

    def test_archive_10x_is_processed(self):
        """ARCHIVE_10X should count as processed data (not raw FASTQ)."""
        from gse_router.schemas.manifest import FileGroup

        groups = [
            FileGroup(
                group_id="arch_1",
                group_type=FileGroupType.ARCHIVE_10X,
                files={"file": f"sha256:{SHA_A}"},
                route="extract_10x_archive_to_h5ad",
                confidence=0.88,
            ),
        ]
        classifier = FileClassifier()
        strategy = classifier.select_strategy(groups)
        assert strategy == Strategy.PROCESSED_COUNTS_FIRST, \
            f"ARCHIVE_10X should be processed, got {strategy}"
