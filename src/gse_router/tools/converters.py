"""Format converters — deterministic tools for converting between single-cell formats.

Each converter is a self-contained tool that:
1. Takes a file group as input
2. Produces an AnnData object and a provenance record as output
3. Is idempotent where possible
4. Reports failures clearly rather than silently producing bad output
"""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

import structlog

logger = structlog.get_logger(__name__)

# Module-level paths to R conversion scripts.
# Resolved relative to the project root (where pyproject.toml lives).
def _find_project_root() -> Path:
    """Find the project root directory by looking for pyproject.toml."""
    current = Path(__file__).resolve().parent
    for _ in range(8):  # search up to 8 levels
        if (current / "pyproject.toml").exists():
            return current
        if (current / "tools").is_dir() and (current / "tools" / "convert_seurat_rds.R").exists():
            return current
        current = current.parent
    # Fallback: assume standard layout (src/gse_router/tools/converters.py -> root)
    return Path(__file__).resolve().parent.parent.parent.parent


_PROJECT_ROOT: Path = _find_project_root()


def _ensure_features_tsv_3cols(
    mtx_dir: Path,
    features_path: Path,
    _temp_links: list[Path],
) -> None:
    """Ensure features.tsv(.gz) has 3 columns (scanpy >= 1.10 requirement).

    Older 10x datasets use ``genes.tsv`` with 2 columns:
        gene_id  gene_symbol

    Newer format uses ``features.tsv`` with 3 columns:
        gene_id  gene_symbol  feature_type

    When a 2-column file is detected, upgrade it in-place to 3 columns
    by appending "Gene Expression" as the feature type for every row.
    """
    import gzip

    for candidate, is_gz in [
        (features_path, features_path.suffix == ".gz"),
        (mtx_dir / "features.tsv", False),
        (mtx_dir / "genes.tsv.gz", True),
        (mtx_dir / "genes.tsv", False),
    ]:
        if not candidate.exists():
            continue
        try:
            opener = gzip.open if is_gz else open
            with opener(str(candidate), "rt") as fh:
                first_line = fh.readline().rstrip("\n")
            ncols = len(first_line.split("\t"))
            if ncols >= 3:
                return  # Already 3+ columns, nothing to do
            # Upgrade to 3 columns: read all, append "Gene Expression", write back
            with opener(str(candidate), "rt") as fh:
                lines = [line.rstrip("\n") + "\tGene Expression" for line in fh]
            with opener(str(candidate), "wt") as fh:
                for line in lines:
                    fh.write(line + "\n")
            _temp_links.append(candidate)  # Track for cleanup
            return
        except Exception:
            continue


def _make_temp_link(source: Path, dest: Path) -> None:
    """Create a temporary hardlink or copy so scanpy can find files with standard names.

    Many GEO datasets use non-standard filenames for MEX components
    (e.g., ``GSM123_barcodes.tsv.gz`` instead of ``barcodes.tsv.gz``).
    scanpy.read_10x_mtx looks for exact filenames in the matrix directory,
    so we create temporary links with standard names.
    """
    try:
        # Try hardlink first (same filesystem, no extra disk space)
        import os
        os.link(str(source), str(dest))
    except (OSError, NotImplementedError):
        # Fallback: copy the file
        shutil.copy2(source, dest)


R_SCRIPTS: dict[str, Path] = {
    "seurat": _PROJECT_ROOT / "tools" / "convert_seurat_rds.R",
    "sce": _PROJECT_ROOT / "tools" / "convert_sce_rds.R",
    "h5seurat": _PROJECT_ROOT / "tools" / "convert_h5seurat.R",
}


class ConversionResult:
    """Result of a format conversion operation."""

    def __init__(
        self,
        success: bool,
        output_path: Path | None = None,
        provenance: dict | None = None,
        error: str | None = None,
        warnings: list[str] | None = None,
    ) -> None:
        self.success = success
        self.output_path = output_path
        self.provenance = provenance or {}
        self.error = error
        self.warnings = warnings or []


def _read_first_line(path: Path) -> str:
    """Read the first line of a file, handling gzip transparently."""
    if path.suffix.lower() == ".gz":
        import gzip
        with gzip.open(path, "rt", encoding="utf-8", errors="replace") as fh:
            return fh.readline()
    else:
        with open(path, encoding="utf-8", errors="replace") as fh:
            return fh.readline()


def _index_looks_like_genes(df: "pd.DataFrame") -> bool:
    """Check if DataFrame index values look like gene symbols/IDs.

    Gene symbols are typically short alphanumeric strings like BRCA1, TP53,
    ENSG00000139618, RP11-34P13.7. Cell barcodes are typically long nucleotide
    sequences (ACGT only) or sample_barcode_# combinations.
    """
    import re
    # Check a sample of index entries (up to 100)
    sample = list(df.index[:100])
    if not sample:
        return False
    gene_score = 0
    barcode_score = 0
    for val in sample:
        val = str(val).strip().strip('"')
        if not val:
            continue
        # Barcode patterns: pure ACGT of 14+ chars, or Sample_Barcode_# format
        if re.match(r"^[ACGT]{14,}$", val, re.IGNORECASE):
            barcode_score += 1
        elif re.match(r"^.+\_[ACGT]{10,}.*\d$", val, re.IGNORECASE):
            barcode_score += 1
        # Gene patterns: short alphanumeric, often with dots/hyphens, starts with letter/digit
        elif re.match(r"^[A-Z0-9][A-Z0-9\-\.]{1,14}$", val, re.IGNORECASE):
            gene_score += 1
        elif re.match(r"^ENS[GT]\d{11}$", val, re.IGNORECASE):  # Ensembl IDs
            gene_score += 1
        elif len(val) <= 15 and not re.match(r"^[ACGT]+$", val, re.IGNORECASE):
            gene_score += 1
        else:
            barcode_score += 1
    return gene_score > barcode_score


def _columns_look_like_cells(df: "pd.DataFrame") -> bool:
    """Check if DataFrame column headers look like cell barcodes/sample IDs.

    Cell barcodes are typically 10x-style 16bp sequences, or Sample_Barcode_#
    compound IDs. Gene symbols almost never appear as column headers in
    single-cell count matrices.
    """
    import re
    sample = list(df.columns[:100])
    if not sample:
        return False
    barcode_score = 0
    gene_score = 0
    for val in sample:
        val = str(val).strip().strip('"')
        if not val:
            continue
        # Cell barcode patterns
        if re.match(r"^[ACGT]{14,}$", val, re.IGNORECASE):
            barcode_score += 1
        elif re.match(r"^.+\_[ACGT]{10,}.*", val, re.IGNORECASE):
            barcode_score += 1
        elif re.match(r"^[ACGT]+\-\d+$", val, re.IGNORECASE):
            barcode_score += 1
        # Gene patterns — should NOT appear as column headers
        elif re.match(r"^[A-Z][A-Z0-9\-\.]{1,14}$", val, re.IGNORECASE):
            gene_score += 1
        elif re.match(r"^ENS[GT]\d{11}$", val, re.IGNORECASE):
            gene_score += 1
        else:
            # Unknown — lean toward barcode (most column headers in scRNA-seq
            # matrices are cell identifiers, not gene symbols)
            barcode_score += 0.5
    return barcode_score > gene_score


def wash_gene_symbols(adata: "AnnData") -> list[str]:
    """Clean gene symbol index by stripping genome prefixes and resolving IDs.

    Ported from v1 ``AutomatedSemanticAnnotator._wash_gene_symbols()``.

    After conversion, gene symbols may contain undesirable artifacts:
    - Genome assembly prefixes: ``mm10_Foxp3`` → ``Foxp3``
    - Ensembl ID concatenations: ``ENSG00000141510_TP53`` → ``TP53``
    - Numeric or Ensembl-only indices where the real gene symbols are in
      a ``.var`` metadata column.

    This function cleans the index in-place and returns a list of
    informational messages describing what was changed.

    Steps:
    1. If the index is predominantly numeric or Ensembl IDs, scan ``.var``
       columns for one containing valid gene symbols and remap.
    2. Strip common genome assembly prefixes from each gene name.
    3. Split Ensembl concatenations on ``_`` and keep the trailing symbol.
    4. Call ``adata.var_names_make_unique()`` to handle collisions.

    Args:
        adata: AnnData object to clean in-place.

    Returns:
        List of info/warning messages about what was changed.
    """
    import re

    import pandas as pd

    messages: list[str] = []
    var_series = pd.Series(adata.var_names.astype(str))

    # --- Step 1: Detect numeric/Ensembl index and remap from var columns ---
    is_numeric = var_series.str.isnumeric().sum() / max(len(var_series), 1) > 0.5
    is_ensembl = (
        var_series.str.startswith("ENS").sum() / max(len(var_series), 1) > 0.5
    )

    if is_numeric or is_ensembl:
        best_col: str | None = None
        best_score = -1
        for col in adata.var.columns:
            col_series = adata.var[col].astype(str)
            valid_mask = (
                (~col_series.str.isnumeric())
                & (~col_series.str.startswith("ENS", na=False))
                & (col_series.str.contains("[a-zA-Z]", regex=True, na=False))
            )
            score = int(valid_mask.sum())
            if score > best_score:
                best_score = score
                best_col = str(col)

        if best_col is not None and best_score > len(adata.var) * 0.05:
            messages.append(
                f"Remapped gene symbols from var column '{best_col}' "
                f"({best_score}/{len(adata.var)} symbols found)"
            )
            # Preserve original names before overwriting
            adata.var["_original_name"] = adata.var_names.astype(str)
            adata.var_names = adata.var[best_col].astype(str)
            adata.var_names_make_unique()
        else:
            messages.append(
                "Warning: numeric/Ensembl index detected; no valid gene symbol column found in var"
            )

    # --- Step 2: Strip genome assembly prefixes ---
    clean_names: list[str] = []
    for gene in adata.var_names:
        g_str = str(gene)
        # Remove common genome assembly prefixes
        g_str = re.sub(
            r"^(GRCh38|GRCm38|GRCm39|hg19|hg38|mm9|mm10|mm39|human|mouse|danRer11|galGal6)[_\-]",
            "",
            g_str,
            flags=re.IGNORECASE,
        )
        # --- Step 3: Split Ensembl concatenations ---
        #   e.g., ENSG00000141510_TP53 → TP53
        if "_" in g_str and g_str.upper().startswith("ENS"):
            parts = g_str.split("_")
            # Keep the last part that doesn't look like an Ensembl ID
            for part in reversed(parts):
                if not re.match(r"^ENS[GT]\d{11}$", part, re.IGNORECASE):
                    g_str = part
                    break
            else:
                g_str = parts[-1]  # fallback: use last segment
        clean_names.append(g_str)

    old_names = list(adata.var_names)
    adata.var_names = clean_names
    adata.var_names_make_unique()

    # Count meaningful changes
    changed = sum(1 for o, n in zip(old_names, adata.var_names) if str(o) != str(n))
    if changed > 0:
        messages.append(f"Cleaned {changed}/{len(old_names)} gene symbols (prefixes, concatenations)")

    return messages


class FormatConverters:
    """Collection of single-cell format conversion tools.

    These are MCP-style executable tools — they do one thing, they do it
    deterministically, and they report provenance.
    """

    def __init__(self, output_dir: Path | None = None) -> None:
        self.output_dir = output_dir or Path.cwd() / "output"
        self._logger = logger.bind(component="FormatConverters")

    # ------------------------------------------------------------------
    # 10x MEX -> h5ad
    # ------------------------------------------------------------------

    def convert_10x_mtx_to_h5ad(
        self,
        mtx_path: Path,
        barcodes_path: Path,
        features_path: Path,
        sample_id: str | None = None,
        var_names: str = "gene_symbols",
    ) -> ConversionResult:
        """Convert a 10x MEX triplet (matrix.mtx + barcodes.tsv + features.tsv)
        to a standard AnnData h5ad file.

        Uses scanpy.read_10x_mtx under the hood.

        Args:
            mtx_path: Path to the MatrixMarket .mtx file.
            barcodes_path: Path to barcodes.tsv(.gz).
            features_path: Path to features/genes.tsv(.gz).
            sample_id: Optional sample identifier.
            var_names: How to name variables ('gene_symbols' or 'gene_ids').

        Returns:
            ConversionResult with path to output .h5ad.
        """
        # --- Pre-flight checks ---
        for file_path, desc in [
            (mtx_path, "matrix.mtx"),
            (barcodes_path, "barcodes.tsv"),
            (features_path, "features.tsv"),
        ]:
            ok, err = self._check_file_exists(file_path, desc)
            if not ok:
                return ConversionResult(success=False, error=err)

        self._logger.info(
            "Converting 10x MEX to h5ad",
            mtx=str(mtx_path),
            sample_id=sample_id,
        )
        try:
            import scanpy as sc

            mtx_dir = mtx_path.parent

            # scanpy.read_10x_mtx expects standard filenames in the matrix directory:
            #   matrix.mtx[.gz], barcodes.tsv[.gz], genes.tsv[.gz] / features.tsv[.gz]
            # Many GEO datasets use non-standard names (e.g. GSM123_barcodes.tsv.gz).
            # When the provided paths differ from standard names, create temporary
            # symlinks/copies so scanpy can find them.
            expected_matrix = mtx_dir / "matrix.mtx.gz"
            expected_barcodes = mtx_dir / "barcodes.tsv.gz"
            expected_features = mtx_dir / "features.tsv.gz"
            _temp_links: list[Path] = []

            # --- Matrix file: handle non-standard names ---
            if not expected_matrix.exists() and not (mtx_dir / "matrix.mtx").exists():
                if mtx_path.exists():
                    _make_temp_link(mtx_path, expected_matrix)
                    _temp_links.append(expected_matrix)

            # --- Barcodes file: handle non-standard names ---
            if not expected_barcodes.exists() and not (mtx_dir / "barcodes.tsv").exists():
                if barcodes_path.exists():
                    _make_temp_link(barcodes_path, expected_barcodes)
                    _temp_links.append(expected_barcodes)

            # --- Features/genes file: handle non-standard names ---
            if not expected_features.exists() and not (mtx_dir / "features.tsv").exists():
                if not (mtx_dir / "genes.tsv.gz").exists() and not (mtx_dir / "genes.tsv").exists():
                    if features_path.exists():
                        _make_temp_link(features_path, expected_features)
                        _temp_links.append(expected_features)

            # --- Old-format genes.tsv (2 columns) → features.tsv (3 columns) ---
            # Older 10x datasets use genes.tsv with only gene_id + gene_symbol.
            # scanpy >= 1.10 expects 3 columns (gene_id, gene_symbol, feature_type).
            # Detect and upgrade the 2-column format.
            _ensure_features_tsv_3cols(mtx_dir, expected_features, _temp_links)

            adata = sc.read_10x_mtx(
                str(mtx_dir),
                var_names=var_names,
                gex_only=True,
            )
            if sample_id:
                adata.obs["sample_id"] = sample_id

            # Preserve raw counts
            adata.layers["counts"] = adata.X.copy()

            wash_msgs = wash_gene_symbols(adata)

            output_path = self._output_path(mtx_path, "10x_mtx", sample_id)
            adata.write_h5ad(output_path)

            # Clean up temp links
            for link in _temp_links:
                try:
                    link.unlink()
                except Exception:
                    pass

            return ConversionResult(
                success=True,
                output_path=output_path,
                provenance={
                    "method": "scanpy.read_10x_mtx",
                    "input_matrix": str(mtx_path),
                    "input_barcodes": str(barcodes_path),
                    "input_features": str(features_path),
                    "var_names": var_names,
                    "n_cells": adata.n_obs,
                    "n_genes": adata.n_vars,
                },
                warnings=wash_msgs if wash_msgs else None,
            )
        except Exception as exc:
            # Clean up temp links on failure too
            for link in _temp_links:
                try:
                    link.unlink()
                except Exception:
                    pass
            self._logger.error("10x MEX conversion failed", error=str(exc))
            return ConversionResult(success=False, error=str(exc))

    # ------------------------------------------------------------------
    # 10x H5 -> h5ad
    # ------------------------------------------------------------------

    def convert_10x_h5_to_h5ad(
        self,
        h5_path: Path,
        sample_id: str | None = None,
        genome: str | None = None,
    ) -> ConversionResult:
        """Convert a 10x HDF5 file to standard AnnData h5ad.

        Uses scanpy.read_10x_h5 under the hood.

        Args:
            h5_path: Path to the .h5 file.
            sample_id: Optional sample identifier.
            genome: Optional genome string for multi-genome H5 files.

        Returns:
            ConversionResult with path to output .h5ad.
        """
        # --- Pre-flight checks ---
        ok, err = self._check_file_exists(h5_path, "10x HDF5 file")
        if not ok:
            return ConversionResult(success=False, error=err)

        # Validate HDF5 magic bytes (\x89HDF\r\n\x1a\n)
        with open(h5_path, "rb") as fh:
            magic = fh.read(8)
        if magic != b"\x89HDF\r\n\x1a\n":
            return ConversionResult(
                success=False,
                error=f"Not a valid HDF5 file (bad magic bytes): '{h5_path}'",
            )

        self._logger.info("Converting 10x H5 to h5ad", h5=str(h5_path), sample_id=sample_id)
        try:
            import scanpy as sc

            adata = sc.read_10x_h5(str(h5_path), genome=genome)
            if sample_id:
                adata.obs["sample_id"] = sample_id

            adata.layers["counts"] = adata.X.copy()

            wash_msgs = wash_gene_symbols(adata)

            output_path = self._output_path(h5_path, "10x_h5", sample_id)
            adata.write_h5ad(output_path)

            return ConversionResult(
                success=True,
                output_path=output_path,
                provenance={
                    "method": "scanpy.read_10x_h5",
                    "input": str(h5_path),
                    "genome": genome,
                    "n_cells": adata.n_obs,
                    "n_genes": adata.n_vars,
                },
                warnings=wash_msgs if wash_msgs else None,
            )
        except Exception as exc:
            self._logger.error("10x H5 conversion failed", error=str(exc))
            return ConversionResult(success=False, error=str(exc))

    # ------------------------------------------------------------------
    # h5ad validate & standardize
    # ------------------------------------------------------------------

    def validate_and_standardize_h5ad(
        self,
        h5ad_path: Path,
        sample_id: str | None = None,
    ) -> ConversionResult:
        """Validate an existing h5ad and standardize its schema.

        Checks for required elements (X, obs, var), ensures raw counts
        are preserved in layers['counts'], and normalizes metadata columns.

        Args:
            h5ad_path: Path to the .h5ad file.
            sample_id: Optional sample identifier to add to obs.

        Returns:
            ConversionResult with path to standardized output .h5ad.
        """
        # --- Pre-flight checks ---
        ok, err = self._check_file_exists(h5ad_path, "h5ad file")
        if not ok:
            return ConversionResult(success=False, error=err)

        # Validate HDF5 magic bytes (h5ad is HDF5-based)
        with open(h5ad_path, "rb") as fh:
            magic = fh.read(8)
        if magic != b"\x89HDF\r\n\x1a\n":
            return ConversionResult(
                success=False,
                error=f"Not a valid HDF5/h5ad file (bad magic bytes): '{h5ad_path}'",
            )

        self._logger.info("Validating and standardizing h5ad", path=str(h5ad_path))
        try:
            import anndata

            adata = anndata.read_h5ad(h5ad_path)

            # Basic validation
            if adata.n_obs == 0:
                return ConversionResult(success=False, error="Empty h5ad: n_obs == 0")
            if adata.n_vars == 0:
                return ConversionResult(success=False, error="Empty h5ad: n_vars == 0")

            # Preserve raw counts if not already in layers
            if "counts" not in adata.layers and adata.X is not None:
                adata.layers["counts"] = adata.X.copy()

            if sample_id:
                adata.obs["sample_id"] = sample_id

            wash_msgs = wash_gene_symbols(adata)

            output_path = self._output_path(h5ad_path, "standardized", sample_id)
            adata.write_h5ad(output_path)

            return ConversionResult(
                success=True,
                output_path=output_path,
                provenance={
                    "method": "validate_and_standardize",
                    "input": str(h5ad_path),
                    "n_cells": adata.n_obs,
                    "n_genes": adata.n_vars,
                    "layers": list(adata.layers.keys()),
                },
                warnings=wash_msgs if wash_msgs else None,
            )
        except Exception as exc:
            self._logger.error("h5ad validation failed", error=str(exc))
            return ConversionResult(success=False, error=str(exc))

    # ------------------------------------------------------------------
    # Loom -> h5ad
    # ------------------------------------------------------------------

    def convert_loom_to_h5ad(
        self,
        loom_path: Path,
        sample_id: str | None = None,
    ) -> ConversionResult:
        """Convert a Loom file to standard AnnData h5ad.

        Args:
            loom_path: Path to the .loom file.
            sample_id: Optional sample identifier.

        Returns:
            ConversionResult with path to output .h5ad.
        """
        # --- Pre-flight checks ---
        ok, err = self._check_file_exists(loom_path, "loom file")
        if not ok:
            return ConversionResult(success=False, error=err)

        # Validate file extension (.loom or .loom.gz)
        loom_name = loom_path.name
        if not (loom_name.endswith(".loom") or loom_name.endswith(".loom.gz")):
            return ConversionResult(
                success=False,
                error=f"Not a recognised loom file extension (expected .loom or .loom.gz): '{loom_path}'",
            )

        self._logger.info("Converting Loom to h5ad", path=str(loom_path))
        try:
            import scanpy as sc

            adata = sc.read_loom(str(loom_path))
            if sample_id:
                adata.obs["sample_id"] = sample_id

            adata.layers["counts"] = adata.X.copy()

            wash_msgs = wash_gene_symbols(adata)

            output_path = self._output_path(loom_path, "loom", sample_id)
            adata.write_h5ad(output_path)

            return ConversionResult(
                success=True,
                output_path=output_path,
                provenance={
                    "method": "scanpy.read_loom",
                    "input": str(loom_path),
                    "n_cells": adata.n_obs,
                    "n_genes": adata.n_vars,
                },
                warnings=wash_msgs if wash_msgs else None,
            )
        except Exception as exc:
            self._logger.error("Loom conversion failed", error=str(exc))
            return ConversionResult(success=False, error=str(exc))

    # ------------------------------------------------------------------
    # Raw UMI matrix -> h5ad (priority 5)
    # ------------------------------------------------------------------

    def convert_raw_umi_matrix_to_h5ad(
        self,
        matrix_path: Path,
        sample_id: str | None = None,
        delimiter: str | None = None,
        transpose: bool = False,
        file_format: str = "txt",
    ) -> ConversionResult:
        """Convert a raw UMI count matrix (TXT/CSV/RDS) to standard AnnData h5ad.

        This is the priority-5 route: explicitly raw counts in a delimited format.
        Differs from generic CSV conversion because we KNOW these are raw counts
        (e.g., from GEO examples like GSE131907 raw_UMI_matrix.txt.gz).

        Args:
            matrix_path: Path to the count matrix file.
            sample_id: Optional sample identifier.
            delimiter: Column delimiter (auto-detected if None).
            transpose: If True, transpose the matrix.
            file_format: 'txt', 'csv', 'tsv', or 'rds'.

        Returns:
            ConversionResult with path to output .h5ad.
        """
        self._logger.info("Converting raw UMI matrix to h5ad", path=str(matrix_path), format=file_format)
        try:
            import anndata
            import numpy as np
            import pandas as pd

            if file_format == "rds":
                # Raw UMI matrix in RDS format — needs R bridge
                return self.convert_rds_to_h5ad(matrix_path, sample_id=sample_id, object_type="matrix")

            if delimiter is None:
                first_line = _read_first_line(matrix_path)
                delimiter = "\t" if "\t" in first_line else ","

            df = pd.read_csv(matrix_path, sep=delimiter, index_col=0)

            # Heuristic: orient genes as rows, cells as columns
            if transpose or (df.shape[0] < df.shape[1]):
                df = df.T

            adata = anndata.AnnData(
                X=df.values.astype(np.float32),
                obs=pd.DataFrame(index=df.index),
                var=pd.DataFrame(index=df.columns),
            )
            if sample_id:
                adata.obs["sample_id"] = sample_id

            # These ARE raw counts — store in X and layers['counts']
            adata.layers["counts"] = adata.X.copy()

            wash_msgs = wash_gene_symbols(adata)

            output_path = self._output_path(matrix_path, "raw_umi", sample_id)
            adata.write_h5ad(output_path)

            return ConversionResult(
                success=True,
                output_path=output_path,
                provenance={
                    "method": "raw_umi_matrix_to_h5ad",
                    "input": str(matrix_path),
                    "delimiter": delimiter,
                    "file_format": file_format,
                    "is_raw_counts": True,
                    "n_cells": adata.n_obs,
                    "n_genes": adata.n_vars,
                },
                warnings=wash_msgs if wash_msgs else None,
            )
        except Exception as exc:
            self._logger.error("Raw UMI matrix conversion failed", error=str(exc))
            return ConversionResult(success=False, error=str(exc))

    # ------------------------------------------------------------------
    # Normalized matrix -> h5ad (priority 7 — limited use)
    # ------------------------------------------------------------------

    def convert_normalized_matrix_to_h5ad(
        self,
        matrix_path: Path,
        sample_id: str | None = None,
        delimiter: str | None = None,
        normalization_method: str = "unknown",
    ) -> ConversionResult:
        """Convert a normalized expression matrix to h5ad — LIMITED USE.

        This handles cases like GSE131907 normalized_log2TPM_matrix.txt.gz.
        The data CAN be loaded but MUST be marked as not-raw-counts.
        Raw counts are NOT preserved in layers['counts'].

        Args:
            matrix_path: Path to the normalized matrix file.
            sample_id: Optional sample identifier.
            delimiter: Column delimiter (auto-detected if None).
            normalization_method: e.g. 'log2TPM', 'FPKM', 'log-normalized'.

        Returns:
            ConversionResult with path to output .h5ad (marked limited-use).
        """
        self._logger.info(
            "Converting normalized matrix to h5ad (limited-use)",
            path=str(matrix_path),
            method=normalization_method,
        )
        try:
            import anndata
            import numpy as np
            import pandas as pd

            if delimiter is None:
                first_line = _read_first_line(matrix_path)
                delimiter = "\t" if "\t" in first_line else ","

            df = pd.read_csv(matrix_path, sep=delimiter, index_col=0)

            if df.shape[0] < df.shape[1]:
                df = df.T

            adata = anndata.AnnData(
                X=df.values.astype(np.float32),
                obs=pd.DataFrame(index=df.index),
                var=pd.DataFrame(index=df.columns),
            )
            if sample_id:
                adata.obs["sample_id"] = sample_id

            # CRITICAL: do NOT copy to layers['counts'] — these are NOT raw counts
            adata.layers["normalized"] = adata.X.copy()
            adata.uns["normalization_method"] = normalization_method
            adata.uns["limited_use"] = True
            adata.uns["limited_use_reason"] = (
                f"Normalized data ({normalization_method}) — raw counts not available. "
                f"Do not use for methods requiring raw counts (e.g., differential expression with DESeq2)."
            )

            wash_msgs = wash_gene_symbols(adata)
            norm_warnings = [
                f"Input is normalized ({normalization_method}), not raw counts. "
                "Output marked as limited-use."
            ]
            if wash_msgs:
                norm_warnings.extend(wash_msgs)

            output_path = self._output_path(matrix_path, "normalized_limited", sample_id)
            adata.write_h5ad(output_path)

            return ConversionResult(
                success=True,
                output_path=output_path,
                provenance={
                    "method": "normalized_matrix_to_h5ad",
                    "input": str(matrix_path),
                    "normalization_method": normalization_method,
                    "is_raw_counts": False,
                    "limited_use": True,
                    "n_cells": adata.n_obs,
                    "n_genes": adata.n_vars,
                },
                warnings=norm_warnings if norm_warnings else None,
            )
        except Exception as exc:
            self._logger.error("Normalized matrix conversion failed", error=str(exc))
            return ConversionResult(success=False, error=str(exc))

    # ------------------------------------------------------------------
    # CSV/TSV matrix -> h5ad (generic — ambiguous)
    # ------------------------------------------------------------------

    def convert_csv_matrix_to_h5ad(
        self,
        csv_path: Path,
        sample_id: str | None = None,
        delimiter: str | None = None,
        transpose: bool = False,
        is_raw_counts: bool | None = None,
        coerce_numeric: bool = True,
    ) -> ConversionResult:
        """Convert a CSV/TSV count matrix to standard AnnData h5ad.

        Attempts to auto-detect gene vs cell orientation by column count.
        If is_raw_counts is unknown, the converter treats data as counts
        but adds a warning that the provenance is uncertain.

        Args:
            csv_path: Path to the CSV/TSV file.
            sample_id: Optional sample identifier.
            delimiter: Column delimiter (auto-detected if None).
            transpose: If True, transpose the matrix.
            is_raw_counts: Whether these are known raw counts. None = unknown.
            coerce_numeric: If True, convert all data columns to numeric with
                ``pd.to_numeric(errors='coerce')`` after CSV read. This prevents
                mixed-type crashes when the CSV contains non-numeric metadata
                rows (e.g., gene annotation footers). Set to False for strict mode
                where non-numeric data should cause a hard error.

        Returns:
            ConversionResult with path to output .h5ad.
        """
        self._logger.info("Converting CSV matrix to h5ad", path=str(csv_path))
        try:
            import anndata
            import numpy as np
            import pandas as pd

            if delimiter is None:
                first_line = _read_first_line(csv_path)
                delimiter = "\t" if "\t" in first_line else ","

            df = pd.read_csv(csv_path, sep=delimiter, index_col=0)

            # Coerce data columns to numeric to handle mixed-type CSVs.
            # Real-world GEO data often has non-numeric footer rows (e.g., gene
            # annotation comments, metadata embedded in count matrices). Without
            # coercion, numpy ops on object-dtype columns crash with "'str' vs 'int'"
            # or "Can't implicitly convert non-string objects to strings".
            if coerce_numeric:
                for col in df.columns:
                    df[col] = pd.to_numeric(df[col], errors="coerce")
                df = df.dropna(how="all")
                if df.empty or df.shape[1] == 0:
                    return ConversionResult(
                        success=False,
                        error=(
                            f"No numeric data found in '{csv_path}'. "
                            "All columns were non-numeric after coercion. "
                            "The file may contain only metadata or annotation rows. "
                            "Try setting coerce_numeric=False for strict parsing."
                        ),
                    )

            # Orientation detection: genes should be var (columns), cells should be obs (rows).
            # In standard CSV format, genes are rows and cells are columns — so
            # after read_csv(index_col=0), genes are df.index and cells are df.columns.
            # We need to transpose if genes are in rows (so they become var/columns).
            #
            # Heuristic: check if the index looks like gene symbols and the columns
            # look like cell barcodes. If both agree, transpose is needed.
            # Wrapped in try/except — mixed-type CSVs can cause these helpers to
            # encounter unexpected values during regex matching.
            try:
                index_looks_gene = _index_looks_like_genes(df)
            except Exception:
                index_looks_gene = False
            try:
                columns_look_cells = _columns_look_like_cells(df)
            except Exception:
                columns_look_cells = False

            if transpose:
                df = df.T
            elif index_looks_gene and columns_look_cells:
                # Genes as rows, cells as columns → transpose to standard orientation
                df = df.T
            elif not index_looks_gene and not columns_look_cells:
                # Ambiguous — fall back to shape heuristic
                if df.shape[0] < df.shape[1]:
                    df = df.T

            # Auto-detect normalized data: raw UMI counts are non-negative integers
            # (or very close). Floats like 0.449, 1.752 indicate log-normalized data.
            if is_raw_counts is None:
                sample_vals = df.values.flatten()
                # Sample up to 100k values for performance
                if len(sample_vals) > 100_000:
                    rng = np.random.default_rng(42)
                    indices = rng.choice(len(sample_vals), size=100_000, replace=False)
                    sample_vals = sample_vals[indices]
                nonzero = sample_vals[sample_vals > 0]
                if len(nonzero) > 0:
                    frac_non_integer = np.mean(~np.isclose(nonzero, np.round(nonzero)))
                    if frac_non_integer > 0.5:
                        is_raw_counts = False
                    elif frac_non_integer < 0.05:
                        is_raw_counts = True

            adata = anndata.AnnData(
                X=df.values.astype(np.float32),
                obs=pd.DataFrame(index=df.index),
                var=pd.DataFrame(index=df.columns),
            )
            if sample_id:
                adata.obs["sample_id"] = sample_id

            warnings: list[str] = []
            if is_raw_counts is True:
                adata.layers["counts"] = adata.X.copy()
            elif is_raw_counts is False:
                adata.layers["normalized"] = adata.X.copy()
                adata.uns["limited_use"] = True
                adata.uns["limited_use_reason"] = (
                    "Auto-detected as normalized (non-integer values). "
                    "Raw counts not available."
                )
                warnings.append(
                    "Auto-detected normalized data (non-integer values) — "
                    "marked limited-use, counts NOT preserved in layers['counts']"
                )
            else:
                # Unknown: preserve as counts but warn
                adata.layers["counts"] = adata.X.copy()
                warnings.append("Raw counts status unknown — verify before using for count-based methods")

            wash_msgs = wash_gene_symbols(adata)
            if wash_msgs:
                warnings.extend(wash_msgs)

            output_path = self._output_path(csv_path, "csv_matrix", sample_id)
            adata.write_h5ad(output_path)

            return ConversionResult(
                success=True,
                output_path=output_path,
                provenance={
                    "method": "csv_matrix_to_h5ad",
                    "input": str(csv_path),
                    "delimiter": delimiter,
                    "is_raw_counts": is_raw_counts,
                    "n_cells": adata.n_obs,
                    "n_genes": adata.n_vars,
                },
                warnings=warnings if warnings else None,
            )
        except Exception as exc:
            self._logger.error("CSV matrix conversion failed", error=str(exc))
            return ConversionResult(success=False, error=str(exc))

    # ------------------------------------------------------------------
    # RDS -> h5ad (via anndata2ri / SeuratDisk)
    # ------------------------------------------------------------------

    def convert_rds_to_h5ad(
        self,
        rds_path: Path,
        sample_id: str | None = None,
        object_type: str = "seurat",
    ) -> ConversionResult:
        """Convert an RDS file (Seurat or SingleCellExperiment) to h5ad.

        Calls an R script via subprocess (Rscript) to perform the conversion.
        Supports Seurat (.rds), SingleCellExperiment (.rds), and h5Seurat
        (.h5seurat) objects.

        Prerequisites:
            R must be installed with Rscript on PATH, and the corresponding
            R packages (Seurat, SeuratDisk, zellkonverter) must be available.

        Args:
            rds_path: Path to the .rds or .h5seurat file.
            sample_id: Optional sample identifier added to obs.
            object_type: 'seurat', 'sce', or 'h5seurat'.

        Returns:
            ConversionResult with path to output .h5ad.
        """
        self._logger.info("Converting RDS to h5ad", path=str(rds_path), object_type=object_type)

        # --- Determine which R script to use ---
        script_path = R_SCRIPTS.get(object_type)
        if script_path is None:
            return ConversionResult(
                success=False,
                error=f"Unknown object_type '{object_type}' — expected 'seurat', 'sce', or 'h5seurat'",
            )

        # --- Check script exists ---
        if not script_path.is_file():
            return ConversionResult(
                success=False,
                error=f"R conversion script not found: '{script_path}'",
            )

        # --- Check Rscript is available ---
        try:
            subprocess.run(
                ["Rscript", "--version"],
                capture_output=True,
                text=True,
                timeout=30,
                check=True,
            )
        except FileNotFoundError:
            return ConversionResult(
                success=False,
                error="Rscript not found on PATH — install R and ensure Rscript is available",
            )
        except subprocess.TimeoutExpired:
            return ConversionResult(
                success=False,
                error="Rscript --version timed out — R installation may be broken",
            )
        except Exception as exc:
            return ConversionResult(
                success=False,
                error=f"Failed to run Rscript: {exc}",
            )

        # --- Build output path ---
        output_path = self._output_path(rds_path, f"rds_{object_type}", sample_id)

        # --- Build the Rscript command ---
        cmd = [
            "Rscript", str(script_path),
            "--input", str(rds_path),
            "--output", str(output_path),
        ]
        if sample_id:
            cmd.extend(["--sample-id", sample_id])
        if object_type == "seurat":
            cmd.extend(["--assay", "RNA"])

        self._logger.info("Running R conversion", cmd=" ".join(cmd))

        # --- Execute ---
        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=1200,  # 20 minutes
            )

            stdout = result.stdout.strip()
            stderr = result.stderr.strip()

            if result.returncode != 0:
                self._logger.error(
                    "R conversion failed",
                    returncode=result.returncode,
                    stderr=stderr,
                )
                return ConversionResult(
                    success=False,
                    error=stderr or f"R script exited with return code {result.returncode}",
                )

            return ConversionResult(
                success=True,
                output_path=output_path,
                provenance={
                    "method": "convert_rds_to_h5ad",
                    "r_script": str(script_path),
                    "object_type": object_type,
                    "input": str(rds_path),
                    "cmd": " ".join(cmd),
                },
            )

        except subprocess.TimeoutExpired:
            self._logger.error("R conversion timed out", path=str(rds_path))
            return ConversionResult(
                success=False,
                error=f"R conversion timed out after 20 minutes: '{rds_path}'",
            )
        except FileNotFoundError:
            return ConversionResult(
                success=False,
                error="Rscript not found on PATH — install R and ensure Rscript is available",
            )
        except Exception as exc:
            self._logger.error("RDS conversion failed", error=str(exc))
            return ConversionResult(success=False, error=str(exc))

    # ------------------------------------------------------------------
    # Merge per-sample h5ad -> combined h5ad
    # ------------------------------------------------------------------

    def merge_h5ad_outputs(
        self,
        h5ad_paths: list[Path],
        dataset_id: str,
        join_method: str = "inner",
    ) -> ConversionResult:
        """Merge multiple per-sample h5ad files into a single combined h5ad.

        Uses anndata.concat() with join='inner' by default.
        Adds dataset_id to adata.obs.
        Preserves layers['counts'] if all inputs have it.

        Args:
            h5ad_paths: List of paths to individual .h5ad files.
            dataset_id: GEO accession for labeling.
            join_method: 'inner' or 'outer' for gene intersection.

        Returns:
            ConversionResult with path to merged .h5ad.
        """
        self._logger.info(
            "Merging h5ad outputs",
            dataset_id=dataset_id,
            n_files=len(h5ad_paths),
            join_method=join_method,
        )
        try:
            import anndata

            # --- Edge case: empty list ---
            if not h5ad_paths:
                return ConversionResult(
                    success=False,
                    error="No h5ad files provided to merge",
                )

            # --- Edge case: single file — just copy and label ---
            if len(h5ad_paths) == 1:
                output_path = self.output_dir / f"{dataset_id}_merged.h5ad"
                self.output_dir.mkdir(parents=True, exist_ok=True)
                shutil.copy2(h5ad_paths[0], output_path)

                adata = anndata.read_h5ad(output_path)
                adata.obs["dataset_id"] = dataset_id
                if "sample_id" not in adata.obs.columns:
                    adata.obs["sample_id"] = f"{dataset_id}_sample_0"
                wash_msgs = wash_gene_symbols(adata)
                adata.write_h5ad(output_path)

                return ConversionResult(
                    success=True,
                    output_path=output_path,
                    provenance={
                        "method": "merge_h5ad_outputs",
                        "n_input_files": 1,
                        "join_method": join_method,
                        "n_cells": adata.n_obs,
                        "n_genes": adata.n_vars,
                    },
                    warnings=wash_msgs if wash_msgs else None,
                )

            # --- Multi-file merge ---
            h5ads: list = []
            all_have_counts = True

            for i, p in enumerate(h5ad_paths):
                adata = anndata.read_h5ad(p)
                if "sample_id" not in adata.obs.columns:
                    adata.obs["sample_id"] = f"{dataset_id}_sample_{i}"
                if "counts" not in adata.layers:
                    all_have_counts = False
                h5ads.append(adata)

            merged = anndata.concat(
                h5ads,
                join=join_method,
                index_unique="-",
            )
            merged.obs["dataset_id"] = dataset_id

            wash_msgs = wash_gene_symbols(merged)

            output_path = self.output_dir / f"{dataset_id}_merged.h5ad"
            merged.write_h5ad(output_path)

            warnings: list[str] = []
            if wash_msgs:
                warnings.extend(wash_msgs)
            if not all_have_counts:
                warnings.append(
                    "Not all inputs have a 'counts' layer — merged object "
                    "may have incomplete or missing counts layer"
                )

            return ConversionResult(
                success=True,
                output_path=output_path,
                provenance={
                    "method": "merge_h5ad_outputs",
                    "n_input_files": len(h5ad_paths),
                    "join_method": join_method,
                    "n_cells": merged.n_obs,
                    "n_genes": merged.n_vars,
                    "all_have_counts": all_have_counts,
                },
                warnings=warnings if warnings else None,
            )
        except Exception as exc:
            self._logger.error("h5ad merge failed", error=str(exc))
            return ConversionResult(success=False, error=str(exc))

    # ------------------------------------------------------------------
    # MATLAB .mat -> h5ad
    # ------------------------------------------------------------------

    def convert_mat_to_h5ad(
        self,
        mat_path: Path,
        sample_id: str | None = None,
        variable_name: str | None = None,
    ) -> ConversionResult:
        """Convert a MATLAB .mat file to .h5ad.

        Uses scipy.io.loadmat to read the matrix, then builds an AnnData object.
        Expects a variable named 'data' or the largest matrix in the file.

        This handles SPLiT-seq and other protocols that distribute count
        matrices as MATLAB files (.mat or .mat.gz).

        Args:
            mat_path: Path to the .mat or .mat.gz file.
            sample_id: Optional sample identifier.
            variable_name: Name of the MATLAB variable containing the count
                matrix. Auto-detected if None (finds largest numeric array).

        Returns:
            ConversionResult with path to output .h5ad.
        """
        self._logger.info("Converting MATLAB .mat to h5ad", path=str(mat_path))
        try:
            import gzip

            import anndata
            import numpy as np
            from scipy.io import loadmat
        except ImportError as exc:
            return ConversionResult(
                success=False,
                error=f"Required package not installed: {exc}. "
                       "Install scipy: pip install scipy",
            )

        try:
            # Handle .mat.gz — scipy can't read gzipped .mat directly
            if mat_path.suffix.lower() == ".gz":
                import tempfile
                with tempfile.NamedTemporaryFile(suffix=".mat", delete=False) as tmp:
                    with gzip.open(mat_path, "rb") as gz_fh:
                        tmp.write(gz_fh.read())
                    tmp_path = Path(tmp.name)
                mat_data = loadmat(str(tmp_path))
                tmp_path.unlink()  # Clean up temp file
            else:
                mat_data = loadmat(str(mat_path))

            # Find the count matrix: use specified variable, or auto-detect
            # by finding the largest numeric array
            numeric_arrays: dict[str, np.ndarray] = {}
            for key, value in mat_data.items():
                if key.startswith("__"):
                    continue  # Skip MATLAB metadata keys
                if isinstance(value, np.ndarray) and np.issubdtype(value.dtype, np.number):
                    numeric_arrays[key] = value

            if not numeric_arrays:
                return ConversionResult(
                    success=False,
                    error=f"No numeric arrays found in .mat file: '{mat_path}'. "
                           f"Available keys: {list(mat_data.keys())}",
                )

            if variable_name and variable_name in numeric_arrays:
                count_matrix = numeric_arrays[variable_name]
            else:
                # Auto-detect: pick the largest numeric array
                largest_key = max(numeric_arrays, key=lambda k: numeric_arrays[k].size)
                count_matrix = numeric_arrays[largest_key]
                self._logger.info(
                    "Auto-detected count matrix variable",
                    variable=largest_key,
                    shape=count_matrix.shape,
                )

            # Build AnnData
            # MATLAB matrices may be genes x cells or cells x genes
            # Heuristic: if the first dimension is > second, it's likely genes x cells
            if count_matrix.shape[0] > count_matrix.shape[1]:
                # Genes x cells — transpose to cells x genes (AnnData convention)
                adata = anndata.AnnData(X=count_matrix.T)
            else:
                adata = anndata.AnnData(X=count_matrix)

            # Preserve raw counts in layers
            adata.layers["counts"] = adata.X.copy()

            if sample_id:
                adata.obs["sample_id"] = sample_id

            wash_msgs = wash_gene_symbols(adata)

            output_path = self._output_path(mat_path, "mat", sample_id)
            adata.write_h5ad(output_path)

            return ConversionResult(
                success=True,
                output_path=output_path,
                provenance={
                    "method": "convert_mat_to_h5ad",
                    "input": str(mat_path),
                    "matrix_shape": list(count_matrix.shape),
                    "variable_name": variable_name or largest_key,
                },
                warnings=wash_msgs if wash_msgs else None,
            )

        except Exception as exc:
            self._logger.error("MATLAB .mat conversion failed", error=str(exc))
            return ConversionResult(success=False, error=str(exc))

    # ------------------------------------------------------------------
    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _check_file_exists(path: Path, description: str) -> tuple[bool, str | None]:
        """Check that a file exists and return (True, None) or (False, error_msg).

        Args:
            path: File path to check.
            description: Human-readable description for error messages.

        Returns:
            Tuple of (exists: bool, error_message: str | None).
        """
        if not path.exists():
            return False, f"File not found: {description} at '{path}'"
        if not path.is_file():
            return False, f"Path is not a regular file: {description} at '{path}'"
        return True, None

    def _output_path(self, input_path: Path, method: str, sample_id: str | None) -> Path:
        """Generate a deterministic output path for a conversion result."""
        stem = input_path.stem
        suffix = f"_{sample_id}" if sample_id else ""
        filename = f"{stem}{suffix}_{method}.h5ad"
        self.output_dir.mkdir(parents=True, exist_ok=True)
        return self.output_dir / filename


# ------------------------------------------------------------------
# CLI: argparse-based entry points for MCP tool dispatch
# ------------------------------------------------------------------

def _result_to_json(result: ConversionResult) -> str:
    """Serialize a ConversionResult to a JSON string for stdout."""
    import json

    payload: dict = {
        "success": result.success,
        "output_path": str(result.output_path) if result.output_path else None,
        "provenance": result.provenance,
        "error": result.error,
        "warnings": result.warnings,
    }
    return json.dumps(payload, indent=2, default=str)


def _add_common_args(parser: object) -> None:
    """Add --output-dir and --sample-id to a subcommand parser."""
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Directory for output .h5ad files (default: ./output)",
    )
    parser.add_argument(
        "--sample-id",
        type=str,
        default=None,
        help="Sample identifier added to AnnData.obs['sample_id']",
    )


def _build_parser() -> object:
    """Build the full argparse parser with all subcommands."""
    import argparse

    parser = argparse.ArgumentParser(
        prog="gse_router.tools.converters",
        description="Convert single-cell RNA-seq data formats to standardised h5ad.",
    )
    sub = parser.add_subparsers(dest="subcommand", required=True)

    # --- convert-10x-mtx ---
    p_10x_mtx = sub.add_parser(
        "convert-10x-mtx",
        help="Convert 10x MEX triplet (matrix.mtx + barcodes.tsv + features.tsv) to h5ad",
    )
    _add_common_args(p_10x_mtx)
    p_10x_mtx.add_argument(
        "--mtx", type=Path, required=True, metavar="PATH",
        help="Path to the MatrixMarket .mtx file",
    )
    p_10x_mtx.add_argument(
        "--barcodes", type=Path, required=True, metavar="PATH",
        help="Path to barcodes.tsv(.gz)",
    )
    p_10x_mtx.add_argument(
        "--features", type=Path, required=True, metavar="PATH",
        help="Path to features/genes.tsv(.gz)",
    )
    p_10x_mtx.add_argument(
        "--var-names", type=str, default="gene_symbols",
        choices=["gene_symbols", "gene_ids"],
        help="How to name variables (default: gene_symbols)",
    )

    # --- convert-10x-h5 ---
    p_10x_h5 = sub.add_parser(
        "convert-10x-h5",
        help="Convert 10x HDF5 (.h5) to h5ad via scanpy.read_10x_h5",
    )
    _add_common_args(p_10x_h5)
    p_10x_h5.add_argument(
        "--h5-path", type=Path, required=True, metavar="PATH",
        help="Path to the 10x .h5 file",
    )
    p_10x_h5.add_argument(
        "--genome", type=str, default=None,
        help="Genome string for multi-genome H5 files",
    )

    # --- convert-loom ---
    p_loom = sub.add_parser(
        "convert-loom",
        help="Convert .loom / .loom.gz to h5ad via scanpy.read_loom",
    )
    _add_common_args(p_loom)
    p_loom.add_argument(
        "--loom-path", type=Path, required=True, metavar="PATH",
        help="Path to the .loom file",
    )

    # --- read-h5ad ---
    p_read = sub.add_parser(
        "read-h5ad",
        help="Read existing .h5ad, validate schema, standardise structure",
    )
    _add_common_args(p_read)
    p_read.add_argument(
        "--h5ad-path", type=Path, required=True, metavar="PATH",
        help="Path to the .h5ad file to validate and standardise",
    )

    # --- convert-raw-umi ---
    p_raw = sub.add_parser(
        "convert-raw-umi",
        help="Parse raw UMI count matrix (TXT/CSV/TSV) and convert to h5ad",
    )
    _add_common_args(p_raw)
    p_raw.add_argument(
        "--matrix", type=Path, required=True, metavar="PATH",
        help="Path to the raw UMI count matrix file",
    )
    p_raw.add_argument(
        "--delimiter", type=str, default=None,
        help="Column delimiter (auto-detected if omitted)",
    )
    p_raw.add_argument(
        "--transpose", action="store_true", default=False,
        help="Transpose the matrix before building AnnData",
    )
    p_raw.add_argument(
        "--file-format", type=str, default="txt",
        choices=["txt", "csv", "tsv", "rds"],
        help="File format hint (default: txt)",
    )

    # --- convert-normalized ---
    p_norm = sub.add_parser(
        "convert-normalized",
        help="Parse normalized matrix (TPM/FPKM/log) — stores in layers['normalized'], marks limited-use",
    )
    _add_common_args(p_norm)
    p_norm.add_argument(
        "--matrix", type=Path, required=True, metavar="PATH",
        help="Path to the normalized expression matrix file",
    )
    p_norm.add_argument(
        "--delimiter", type=str, default=None,
        help="Column delimiter (auto-detected if omitted)",
    )
    p_norm.add_argument(
        "--normalization-method", type=str, default="unknown",
        help="Normalization method label, e.g. log2TPM, FPKM, log-normalized (default: unknown)",
    )

    # --- convert-csv ---
    p_csv = sub.add_parser(
        "convert-csv",
        help="Parse ambiguous CSV/TSV count matrix with gene/cell orientation inference",
    )
    _add_common_args(p_csv)
    p_csv.add_argument(
        "--csv", type=Path, required=True, metavar="PATH",
        help="Path to the CSV/TSV count matrix file",
    )
    p_csv.add_argument(
        "--delimiter", type=str, default=None,
        help="Column delimiter (auto-detected if omitted)",
    )
    p_csv.add_argument(
        "--transpose", action="store_true", default=False,
        help="Transpose the matrix before building AnnData",
    )
    p_csv.add_argument(
        "--is-raw-counts", type=str, default="unknown",
        choices=["true", "false", "unknown"],
        help="Are these raw counts? 'true', 'false', or 'unknown' (default: unknown)",
    )

    return parser


def _resolve_bool_tristate(value: str | None) -> bool | None:
    """Map CLI string choices to Python bool | None for is_raw_counts."""
    if value is None or value == "unknown":
        return None
    if value == "true":
        return True
    if value == "false":
        return False
    return None


def main(args: list[str] | None = None) -> None:
    """Entry point for CLI subcommands dispatched by the MCP server.

    Each subcommand constructs a FormatConverters instance, calls the
    corresponding method, prints JSON result to stdout, and exits with
    code 0 (success) or 1 (failure).
    """
    import sys

    parser = _build_parser()
    ns = parser.parse_args(args)

    converters = FormatConverters(output_dir=ns.output_dir)

    # Dispatch map: subcommand -> (method_name, kwargs builder)
    if ns.subcommand == "convert-10x-mtx":
        result = converters.convert_10x_mtx_to_h5ad(
            mtx_path=ns.mtx,
            barcodes_path=ns.barcodes,
            features_path=ns.features,
            sample_id=ns.sample_id,
            var_names=ns.var_names,
        )
    elif ns.subcommand == "convert-10x-h5":
        result = converters.convert_10x_h5_to_h5ad(
            h5_path=getattr(ns, "h5_path"),
            sample_id=ns.sample_id,
            genome=ns.genome,
        )
    elif ns.subcommand == "convert-loom":
        result = converters.convert_loom_to_h5ad(
            loom_path=getattr(ns, "loom_path"),
            sample_id=ns.sample_id,
        )
    elif ns.subcommand == "read-h5ad":
        result = converters.validate_and_standardize_h5ad(
            h5ad_path=getattr(ns, "h5ad_path"),
            sample_id=ns.sample_id,
        )
    elif ns.subcommand == "convert-raw-umi":
        result = converters.convert_raw_umi_matrix_to_h5ad(
            matrix_path=ns.matrix,
            sample_id=ns.sample_id,
            delimiter=ns.delimiter,
            transpose=ns.transpose,
            file_format=ns.file_format,
        )
    elif ns.subcommand == "convert-normalized":
        result = converters.convert_normalized_matrix_to_h5ad(
            matrix_path=ns.matrix,
            sample_id=ns.sample_id,
            delimiter=ns.delimiter,
            normalization_method=ns.normalization_method,
        )
    elif ns.subcommand == "convert-csv":
        result = converters.convert_csv_matrix_to_h5ad(
            csv_path=getattr(ns, "csv"),
            sample_id=ns.sample_id,
            delimiter=ns.delimiter,
            transpose=ns.transpose,
            is_raw_counts=_resolve_bool_tristate(ns.is_raw_counts),
        )
    else:
        # argparse should prevent this, but guard anyway
        print(_result_to_json(ConversionResult(success=False, error=f"Unknown subcommand: {ns.subcommand}")))
        sys.exit(1)

    print(_result_to_json(result))
    sys.exit(0 if result.success else 1)


if __name__ == "__main__":
    main()
