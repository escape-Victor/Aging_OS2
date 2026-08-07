"""Unit tests for FormatConverters — CSV type handling and numeric coercion."""

from __future__ import annotations

import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from gse_router.tools.converters import FormatConverters, ConversionResult


class TestCsvMixedTypeHandling:
    """Fix 2: CSV mixed-type crashes when pd.read_csv produces object-dtype columns."""

    @pytest.fixture
    def converters(self, tmp_path):
        return FormatConverters(output_dir=tmp_path)

    def _write_csv(self, dir_path: Path, name: str, data: str) -> Path:
        path = dir_path / name
        path.write_text(data)
        return path

    # ------------------------------------------------------------------
    # coerce_numeric=True (default) — handles mixed types gracefully
    # ------------------------------------------------------------------

    def test_numeric_coercion_handles_string_footer(self, converters, tmp_path):
        """CSV with a non-numeric footer row should not crash."""
        csv_path = self._write_csv(tmp_path, "with_footer.csv", (
            "gene,cell1,cell2,cell3\n"
            "BRCA1,10,5,0\n"
            "TP53,0,3,8\n"
            "EGFR,2,0,1\n"
            "# comment row with text,skipped,by,coercion\n"
        ))
        result = converters.convert_csv_matrix_to_h5ad(csv_path, coerce_numeric=True)
        assert result.success, f"Expected success but got: {result.error}"
        assert result.output_path.exists()
        assert result.provenance["n_cells"] == 3
        assert result.provenance["n_genes"] == 3

    def test_numeric_coercion_handles_mixed_column(self, converters, tmp_path):
        """Column with both numbers and strings should be coerced, not crash."""
        csv_path = self._write_csv(tmp_path, "mixed_col.csv", (
            "gene,cell1,cell2,cell3\n"
            "BRCA1,10,5,0\n"
            "TP53,NA,3,8\n"
            "EGFR,2,invalid,1\n"
            "GAPDH,0,0,0\n"
        ))
        result = converters.convert_csv_matrix_to_h5ad(csv_path, coerce_numeric=True)
        assert result.success, f"Expected success but got: {result.error}"
        assert result.output_path.exists()

    def test_all_non_numeric_returns_clear_error(self, converters, tmp_path):
        """A CSV with only non-numeric data should produce a clear error."""
        csv_path = self._write_csv(tmp_path, "all_metadata.csv", (
            "gene,description,category,notes\n"
            "BRCA1,breast cancer,oncogene,important\n"
            "TP53,tumor suppressor,tumor,key\n"
        ))
        result = converters.convert_csv_matrix_to_h5ad(csv_path, coerce_numeric=True)
        assert not result.success
        assert "No numeric data found" in result.error
        assert "coerce_numeric" in result.error

    def test_coercion_preserves_integer_counts(self, converters, tmp_path):
        """Integer count values should survive coercion unchanged."""
        csv_path = self._write_csv(tmp_path, "counts.csv", (
            "gene,cell1,cell2,cell3,cell4\n"
            "BRCA1,10,5,0,2\n"
            "TP53,0,3,8,1\n"
            "EGFR,2,0,1,0\n"
        ))
        result = converters.convert_csv_matrix_to_h5ad(csv_path, coerce_numeric=True)
        assert result.success
        # Read back and verify counts
        import anndata
        adata = anndata.read_h5ad(result.output_path)
        assert "counts" in adata.layers
        assert adata.n_obs + adata.n_vars == 7  # 3 genes + 4 cells
        # BRCA1 should be in either obs_names or var_names
        assert "BRCA1" in list(adata.obs_names) + list(adata.var_names)

    # ------------------------------------------------------------------
    # coerce_numeric=False — strict mode
    # ------------------------------------------------------------------

    def test_strict_mode_preserves_original_behavior(self, converters, tmp_path):
        """With coerce_numeric=False, well-formed CSV should still work."""
        csv_path = self._write_csv(tmp_path, "clean.csv", (
            "gene,cell1,cell2,cell3\n"
            "BRCA1,10,5,0\n"
            "TP53,0,3,8\n"
            "EGFR,2,0,1\n"
        ))
        result = converters.convert_csv_matrix_to_h5ad(csv_path, coerce_numeric=False)
        assert result.success
        assert result.output_path.exists()

    def test_tsv_with_tab_delimiter(self, converters, tmp_path):
        """Tab-delimited CSV should also get coercion."""
        csv_path = self._write_csv(tmp_path, "counts.tsv", (
            "gene\tcell1\tcell2\tcell3\n"
            "BRCA1\t10\t5\t0\n"
            "TP53\t0\t3\t8\n"
            "comment_row\ttext\there\tnow\n"
        ))
        result = converters.convert_csv_matrix_to_h5ad(csv_path, coerce_numeric=True)
        assert result.success
        assert result.provenance["delimiter"] == "\t"

    # ------------------------------------------------------------------
    # Orientation detection robustness
    # ------------------------------------------------------------------

    def test_orientation_detection_survives_mixed_types(self, converters, tmp_path):
        """_index_looks_like_genes/_columns_look_like_cells wrapped in try/except."""
        csv_path = self._write_csv(tmp_path, "odd_index.csv", (
            ",cell1,cell2,cell3\n"
            "BRCA1,10,5,0\n"
            "TP53,0,3,8\n"
            "EGFR,2,0,1\n"
        ))
        result = converters.convert_csv_matrix_to_h5ad(csv_path, coerce_numeric=True)
        assert result.success

    def test_empty_dataframe_after_coercion(self, converters, tmp_path):
        """DataFrame with only NaN after coercion should error clearly."""
        csv_path = self._write_csv(tmp_path, "all_text.csv", (
            "gene,comment1,comment2\n"
            "note,hello,world\n"
            "footnote,foo,bar\n"
        ))
        result = converters.convert_csv_matrix_to_h5ad(csv_path, coerce_numeric=True)
        assert not result.success
        assert "No numeric data" in result.error

    # ------------------------------------------------------------------
    # Large-value and edge cases
    # ------------------------------------------------------------------

    def test_large_integers_preserved(self, converters, tmp_path):
        """Large UMI counts should not be corrupted by float conversion."""
        csv_path = self._write_csv(tmp_path, "large_counts.csv", (
            "gene,cell1,cell2\n"
            "BRCA1,999999,0\n"
            "TP53,0,888888\n"
            "EGFR,123456,0\n"
        ))
        result = converters.convert_csv_matrix_to_h5ad(csv_path, coerce_numeric=True)
        assert result.success
        import anndata
        adata = anndata.read_h5ad(result.output_path)
        assert adata.X.max() >= 999999

    def test_single_cell_csv(self, converters, tmp_path):
        """Single-cell CSV (one data column) should work."""
        csv_path = self._write_csv(tmp_path, "one_cell.csv", (
            "gene,cell1\n"
            "BRCA1,10\n"
            "TP53,5\n"
        ))
        result = converters.convert_csv_matrix_to_h5ad(csv_path, coerce_numeric=True)
        assert result.success
        # Orientation: index_looks_gene=True, columns_look_cells=False
        # (cell1 is not a valid barcode pattern). Shape heuristic 2 < 1 → no transpose.
        # Result stays as-is: genes=obs(index), cells=var(columns)
        assert result.provenance["n_cells"] == 2
        assert result.provenance["n_genes"] == 1

    def test_single_gene_csv(self, converters, tmp_path):
        """Single-gene CSV (one data row) should work."""
        csv_path = self._write_csv(tmp_path, "one_gene.csv", (
            "gene,cell1,cell2,cell3\n"
            "BRCA1,10,5,0\n"
        ))
        result = converters.convert_csv_matrix_to_h5ad(csv_path, coerce_numeric=True)
        assert result.success
        assert result.output_path.exists()
        assert result.provenance["n_cells"] + result.provenance["n_genes"] == 4
