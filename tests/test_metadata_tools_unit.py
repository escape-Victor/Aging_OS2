"""Unit tests for MetadataTools — join column detection and error handling."""

from __future__ import annotations

import tempfile
from pathlib import Path

import anndata
import numpy as np
import pandas as pd
import pytest

from gse_router.tools.metadata_tools import (
    MetadataTools,
    MetadataMergeResult,
    _auto_detect_join_column,
)


class TestAutoDetectJoinColumn:
    """Fix 3: _auto_detect_join_column returns None (not 'barcode') when no match."""

    def test_returns_barcode_when_present(self):
        """Should detect 'barcode' column when it exists."""
        meta_df = pd.DataFrame({
            "barcode": ["AAACCCAAG", "AAACCCAGT"],
            "condition": ["treated", "control"],
        })
        result = _auto_detect_join_column(meta_df, pd.Index(["AAACCCAAG", "AAACCCAGT"]))
        assert result == "barcode"

    def test_returns_cell_when_barcode_missing(self):
        """Should detect 'cell' column when 'barcode' not present."""
        meta_df = pd.DataFrame({
            "cell": ["AAACCCAAG", "AAACCCAGT"],
            "condition": ["treated", "control"],
        })
        result = _auto_detect_join_column(meta_df, pd.Index(["AAACCCAAG", "AAACCCAGT"]))
        assert result == "cell"

    def test_returns_none_when_no_match(self):
        """Should return None when no column matches — NOT hard-coded 'barcode'."""
        meta_df = pd.DataFrame({
            "sample": ["S1", "S2"],
            "condition": ["treated", "control"],
            "batch": ["1", "2"],
        })
        obs_index = pd.Index(["AAACCCAAG", "AAACCCAGT"])
        result = _auto_detect_join_column(meta_df, obs_index)
        assert result is None, (
            f"Expected None (no match found) but got '{result}'. "
            "Hard-coded fallback to 'barcode' was removed."
        )

    def test_returns_none_for_empty_metadata(self):
        """Should return None for empty metadata DataFrame."""
        meta_df = pd.DataFrame()
        result = _auto_detect_join_column(meta_df, pd.Index(["AAACCCAAG"]))
        assert result is None

    def test_overlap_detection_finds_best_column(self):
        """Overlap-based fallback should find column with highest index overlap."""
        meta_df = pd.DataFrame({
            "id": ["AAACCCAAG", "BBBBCCDD"],
            "extra": ["X", "Y"],
        })
        obs_index = pd.Index(["AAACCCAAG", "AAACCCAGT"])
        result = _auto_detect_join_column(meta_df, obs_index)
        assert result == "id", f"Expected 'id' (highest overlap) but got '{result}'"

    def test_case_insensitive_candidates(self):
        """Candidate column names are checked case-sensitively first, then overlap."""
        meta_df = pd.DataFrame({
            "Barcode": ["AAACCCAAG", "AAACCCAGT"],
            "Condition": ["treated", "control"],
        })
        result = _auto_detect_join_column(meta_df, pd.Index(["AAACCCAAG", "AAACCCAGT"]))
        assert result == "Barcode"

    def test_cell_name_variants(self):
        """CellName/cellname/cell_id variants should all be detected."""
        for col in ["CellName", "cellname", "cell_id", "CellID", "cellid"]:
            meta_df = pd.DataFrame({
                col: ["AAACCCAAG", "AAACCCAGT"],
                "extra": [1, 2],
            })
            result = _auto_detect_join_column(meta_df, pd.Index(["AAACCCAAG", "AAACCCAGT"]))
            assert result == col, f"Should detect '{col}'"


class TestMergeCellMetadata:
    """Fix 3: merge_cell_metadata produces clear errors when no join column found."""

    @pytest.fixture
    def tools(self):
        return MetadataTools()

    @pytest.fixture
    def sample_h5ad(self, tmp_path):
        """Create a minimal h5ad file for testing."""
        adata = anndata.AnnData(
            X=np.eye(3, dtype=np.float32),
            obs=pd.DataFrame(index=["AAACCCAAG", "AAACCCAGT", "AAACCCATG"]),
            var=pd.DataFrame(index=["BRCA1", "TP53", "EGFR"]),
        )
        path = tmp_path / "test.h5ad"
        adata.write_h5ad(path)
        return path

    def test_merge_with_auto_detected_column(self, tools, sample_h5ad, tmp_path):
        """Merge should succeed when barcode column is auto-detected."""
        meta_path = tmp_path / "meta.csv"
        meta_path.write_text("barcode,condition\nAAACCCAAG,treated\nAAACCCAGT,control\n")

        result = tools.merge_cell_metadata(sample_h5ad, meta_path)
        assert result.success, f"Expected success but got: {result.error}"
        assert result.n_cells_matched > 0

    def test_merge_fails_clearly_when_no_column_found(self, tools, sample_h5ad, tmp_path):
        """Clear error with column list when no join column can be found."""
        meta_path = tmp_path / "meta_no_barcode.csv"
        meta_path.write_text("sample,condition,batch\nS1,treated,1\nS2,control,2\n")

        result = tools.merge_cell_metadata(sample_h5ad, meta_path)
        assert not result.success
        assert "No join column found" in result.error
        assert "sample" in result.error or "Available columns" in result.error

    def test_merge_with_explicit_join_column(self, tools, sample_h5ad, tmp_path):
        """Explicit join_column should be used when provided."""
        meta_path = tmp_path / "meta_custom.csv"
        meta_path.write_text("cell_id,condition\nAAACCCAAG,treated\nAAACCCAGT,control\n")

        result = tools.merge_cell_metadata(sample_h5ad, meta_path, join_column="cell_id")
        assert result.success, f"Expected success but got: {result.error}"

    def test_merge_exception_includes_column_context(self, tools, sample_h5ad, tmp_path):
        """Exception handler should include column info for agent debugging."""
        meta_path = tmp_path / "bad_meta.csv"
        meta_path.write_text("x,y\n1,2\n")  # Won't match anything

        result = tools.merge_cell_metadata(sample_h5ad, meta_path, join_column="nonexistent")
        assert not result.success
        # The error should contain column info since auto-detection was attempted
        assert "x" in result.error.lower() or "column" in result.error.lower() or "nonexistent" in result.error.lower()

    def test_merge_preserves_existing_obs_columns(self, tools, sample_h5ad, tmp_path):
        """Existing obs columns should survive the merge."""
        # Add a pre-existing column
        import anndata as ad
        adata = ad.read_h5ad(sample_h5ad)
        adata.obs["pre_existing"] = ["a", "b", "c"]
        adata.write_h5ad(sample_h5ad)

        meta_path = tmp_path / "meta.csv"
        meta_path.write_text("barcode,condition\nAAACCCAAG,treated\n")

        result = tools.merge_cell_metadata(sample_h5ad, meta_path)
        assert result.success
        result_adata = ad.read_h5ad(result.output_path)
        assert "pre_existing" in result_adata.obs.columns
