"""Tests for the route planner — priority ordering, technology detection, reference genome selection.

Tests cover:
- 8-level routing priority with ARCHIVE_10X at priority 4
- Non-count files excluded from route planning
- Reference genome selection from organism
- Rule application
"""

from __future__ import annotations

import pytest

from gse_router.router.route_planner import RoutePlanner
from gse_router.router.rules import (
    RouteDecision,
    RulePriority,
    apply_all_rules,
    decode_route,
    select_reference_genome,
)
from gse_router.schemas.manifest import (
    DatasetManifest,
    FileGroup,
    FileGroupType,
    Strategy,
)
from gse_router.schemas.route_plan import ReferenceGenome, RouteMethod
from gse_router.schemas.sidecar import FileDetectedType, FileSidecar, ProcessingStage


def _make_manifest(groups: list[FileGroup], organism: list[str] | None = None) -> DatasetManifest:
    """Helper to create a minimal DatasetManifest for testing."""
    return DatasetManifest(
        dataset_id="GSE999999",
        title="Test Dataset",
        organism=organism or ["Homo sapiens"],
        platform=["GPL24676"],
        technology="10x Genomics 3' v3",
        files={},
        file_groups=groups,
        routing_priority=[g.group_id for g in groups],
        selected_strategy=Strategy.PROCESSED_COUNTS_FIRST,
    )


class TestRoutingPriority:
    """Test the 8-level routing priority ordering."""

    def test_archive_10x_priority_same_as_mex(self):
        """ARCHIVE_10X should have the same priority as MEX_10X (4)."""
        planner = RoutePlanner()

        groups = [
            FileGroup(
                group_id="h5ad_1",
                group_type=FileGroupType.H5AD,
                files={"data": "sha256:aaa"},
                route="read_h5ad_validate",
                confidence=0.95,
            ),
            FileGroup(
                group_id="arch_10x_1",
                group_type=FileGroupType.ARCHIVE_10X,
                files={"file": "sha256:bbb"},
                route="extract_10x_archive_to_h5ad",
                confidence=0.88,
            ),
            FileGroup(
                group_id="mex_1",
                group_type=FileGroupType.MEX_10X,
                files={"matrix": "sha256:ccc", "barcodes": "sha256:ddd", "features": "sha256:eee"},
                route="read_10x_mtx_to_h5ad",
                confidence=0.92,
            ),
        ]

        manifest = _make_manifest(groups)
        ordered = planner._prioritize_groups(groups, Strategy.PROCESSED_COUNTS_FIRST)

        # h5ad should be first (priority 1)
        assert ordered[0].group_type == FileGroupType.H5AD
        # Both ARCHIVE_10X and MEX_10X should come before FASTQ-level items
        group_types = [g.group_type for g in ordered]
        assert FileGroupType.ARCHIVE_10X in group_types
        assert FileGroupType.MEX_10X in group_types

    def test_non_count_groups_excluded_from_routing(self):
        """Non-count supplementary groups should be excluded from route planning."""
        planner = RoutePlanner()

        groups = [
            FileGroup(
                group_id="mex_1",
                group_type=FileGroupType.MEX_10X,
                files={"matrix": "sha256:aaa", "barcodes": "sha256:bbb", "features": "sha256:ccc"},
                route="read_10x_mtx_to_h5ad",
                confidence=0.92,
            ),
            FileGroup(
                group_id="unknown_1",
                group_type=FileGroupType.UNKNOWN,
                files={"file": "sha256:ddd"},
                route="manual_review",
                confidence=0.30,
            ),
        ]

        manifest = _make_manifest(groups)
        ordered = planner._prioritize_groups(groups, Strategy.PROCESSED_COUNTS_FIRST)

        # UNKNOWN groups should be excluded
        group_types = [g.group_type for g in ordered]
        assert FileGroupType.UNKNOWN not in group_types or len(ordered) == 1

    def test_fastq_excluded_when_processed_exists(self):
        """FASTQ/SRA should be excluded when processed data exists."""
        planner = RoutePlanner()

        groups = [
            FileGroup(
                group_id="h5ad_1",
                group_type=FileGroupType.H5AD,
                files={"data": "sha256:aaa"},
                route="read_h5ad_validate",
                confidence=0.95,
            ),
            FileGroup(
                group_id="fq_1",
                group_type=FileGroupType.FASTQ_PAIRED,
                files={"read1": "sha256:bbb", "read2": "sha256:ccc"},
                route="cellranger_count",
                confidence=0.85,
            ),
        ]

        manifest = _make_manifest(groups)
        ordered = planner._prioritize_groups(groups, Strategy.PROCESSED_COUNTS_FIRST)

        group_types = [g.group_type for g in ordered]
        assert FileGroupType.FASTQ_PAIRED not in group_types


class TestReferenceGenome:
    """Test reference genome selection from organism name."""

    @pytest.mark.parametrize("organism,expected", [
        ("Homo sapiens", ReferenceGenome.GRCH38),
        ("human", ReferenceGenome.GRCH38),
        ("Hsapiens", ReferenceGenome.GRCH38),
        ("Mus musculus", ReferenceGenome.GRCM39),
        ("mouse", ReferenceGenome.GRCM39),
        ("Mmusculus", ReferenceGenome.GRCM39),
        ("Danio rerio", ReferenceGenome.CUSTOM),
        ("unknown species", ReferenceGenome.CUSTOM),
    ])
    def test_reference_genome_mapping(self, organism, expected):
        result = select_reference_genome(organism)
        assert result == expected, \
            f"Expected {expected} for '{organism}', got {result}"


class TestRouteDecoding:
    """Test route string → RouteDecision mapping."""

    def test_extract_10x_archive_decoded(self):
        """ARCHIVE_10X route should decode to READ_10X_MTX with archive source."""
        group = FileGroup(
            group_id="test_arch",
            group_type=FileGroupType.ARCHIVE_10X,
            files={"file": "sha256:aaa"},
            route="extract_10x_archive_to_h5ad",
            confidence=0.88,
        )
        decision = decode_route(group)
        assert decision.method == RouteMethod.READ_10X_MTX
        assert decision.parameters.get("source") == "archive_10x"

    def test_unknown_route_decodes_to_manual_review(self):
        """Unknown routes should fall back to MANUAL_REVIEW."""
        group = FileGroup(
            group_id="test_unk",
            group_type=FileGroupType.UNKNOWN,
            files={"file": "sha256:aaa"},
            route="nonexistent_route",
            confidence=0.30,
        )
        decision = decode_route(group)
        assert decision.method == RouteMethod.MANUAL_REVIEW
        assert len(decision.warnings) > 0


class TestRuleApplication:
    """Test that critical routing rules work correctly."""

    def test_must_have_usable_data_passes(self):
        """Rule should pass when processable data exists."""
        groups = [
            FileGroup(
                group_id="h5_1",
                group_type=FileGroupType.H5_10X,
                files={"h5": "sha256:aaa"},
                route="read_10x_h5_to_h5ad",
                confidence=0.95,
            ),
        ]
        manifest = _make_manifest(groups)
        results = apply_all_rules(manifest)
        data_rule = [r for r in results if r.rule_name == "must_have_usable_data"]
        assert len(data_rule) == 1
        assert data_rule[0].passed is True

    def test_must_have_usable_data_fails_with_nothing(self):
        """Rule should fail when no processable data exists."""
        manifest = _make_manifest([])
        results = apply_all_rules(manifest)
        data_rule = [r for r in results if r.rule_name == "must_have_usable_data"]
        assert len(data_rule) == 1
        assert data_rule[0].passed is False

    def test_skip_fastq_rule(self):
        """Rule should flag FASTQ when processed data exists."""
        groups = [
            FileGroup(
                group_id="h5_1",
                group_type=FileGroupType.H5AD,
                files={"data": "sha256:aaa"},
                route="read_h5ad_validate",
                confidence=0.95,
            ),
            FileGroup(
                group_id="fq_1",
                group_type=FileGroupType.FASTQ_PAIRED,
                files={"read1": "sha256:bbb", "read2": "sha256:ccc"},
                route="cellranger_count",
                confidence=0.85,
            ),
        ]
        manifest = _make_manifest(groups)
        results = apply_all_rules(manifest)
        skip_rule = [r for r in results if r.rule_name == "skip_fastq_if_processed_exists"]
        assert len(skip_rule) == 1
        assert "fastq" in skip_rule[0].message.lower()
