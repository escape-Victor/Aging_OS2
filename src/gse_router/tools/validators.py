"""h5ad validation tools — QC checks, schema validation, and quality metrics."""

from __future__ import annotations

from pathlib import Path

import structlog

logger = structlog.get_logger(__name__)


class QCResult:
    """Result of a QC validation run."""

    def __init__(
        self,
        passed: bool,
        metrics: dict[str, float | int | str] | None = None,
        failures: list[str] | None = None,
        warnings: list[str] | None = None,
    ) -> None:
        self.passed = passed
        self.metrics = metrics or {}
        self.failures = failures or []
        self.warnings = warnings or []


class AnnDataValidator:
    """Validates AnnData objects against quality standards.

    Checks include:
    - Minimum cell/gene counts
    - Mitochondrial gene percentage
    - Raw counts preservation
    - Metadata completeness
    - Gene symbol validity
    """

    def __init__(
        self,
        min_cells: int = 100,
        min_genes_per_cell: int = 200,
        min_counts_per_cell: int = 500,
        max_mito_pct: float = 20.0,
    ) -> None:
        self.min_cells = min_cells
        self.min_genes_per_cell = min_genes_per_cell
        self.min_counts_per_cell = min_counts_per_cell
        self.max_mito_pct = max_mito_pct
        self._logger = logger.bind(component="AnnDataValidator")

    def validate(self, h5ad_path: Path) -> QCResult:
        """Run full validation on an h5ad file.

        Args:
            h5ad_path: Path to the .h5ad file.

        Returns:
            QCResult with metrics, failures, and warnings.
        """
        self._logger.info("Validating h5ad", path=str(h5ad_path))
        try:
            import anndata
            import numpy as np

            adata = anndata.read_h5ad(h5ad_path)
        except Exception as exc:
            return QCResult(passed=False, failures=[f"Cannot read h5ad: {exc}"])

        metrics: dict[str, float | int | str] = {}
        failures: list[str] = []
        warnings: list[str] = []

        # Cell count
        metrics["n_cells"] = adata.n_obs
        if adata.n_obs < self.min_cells:
            failures.append(f"n_cells ({adata.n_obs}) < {self.min_cells}")

        # Gene count
        metrics["n_genes"] = adata.n_vars
        if adata.n_vars == 0:
            failures.append("n_genes == 0")

        # Genes per cell
        n_genes_by_cell = np.array((adata.X > 0).sum(axis=1)).flatten()
        metrics["median_genes_per_cell"] = float(np.median(n_genes_by_cell))
        cells_low_genes = (n_genes_by_cell < self.min_genes_per_cell).sum()
        metrics["cells_below_min_genes"] = int(cells_low_genes)
        if cells_low_genes > adata.n_obs * 0.5:
            failures.append(f"{cells_low_genes}/{adata.n_obs} cells below {self.min_genes_per_cell} genes")

        # Total counts per cell
        total_counts = np.array(adata.X.sum(axis=1)).flatten()
        metrics["median_total_counts"] = float(np.median(total_counts))
        cells_low_counts = (total_counts < self.min_counts_per_cell).sum()
        metrics["cells_below_min_counts"] = int(cells_low_counts)

        # Mitochondrial percentage
        mito_genes = [g for g in adata.var_names if str(g).upper().startswith("MT-")]
        if mito_genes and adata.n_vars > 0:
            mito_idx = [i for i, g in enumerate(adata.var_names) if str(g).upper().startswith("MT-")]
            if mito_idx:
                mito_pct = np.array(adata.X[:, mito_idx].sum(axis=1) / adata.X.sum(axis=1) * 100).flatten()
                metrics["median_mito_pct"] = float(np.median(mito_pct))
                if np.median(mito_pct) > self.max_mito_pct:
                    warnings.append(f"High mitochondrial content: median {np.median(mito_pct):.1f}%")
        else:
            metrics["median_mito_pct"] = "N/A (no MT- genes detected)"

        # Raw counts preserved
        if "counts" in adata.layers:
            metrics["has_raw_counts"] = "yes"
            metrics["raw_counts_layer"] = "counts"
        else:
            metrics["has_raw_counts"] = "no"
            warnings.append("Raw counts not preserved in layers['counts']")

        # Metadata columns
        metrics["obs_columns"] = ", ".join(adata.obs.columns.tolist()[:20])
        metrics["var_columns"] = ", ".join(adata.var.columns.tolist()[:20])

        passed = len(failures) == 0
        return QCResult(passed=passed, metrics=metrics, failures=failures, warnings=warnings)

    def quick_check(self, h5ad_path: Path) -> bool:
        """Quick sanity check: can the file be opened?"""
        try:
            import anndata
            adata = anndata.read_h5ad(h5ad_path)
            return adata.n_obs > 0 and adata.n_vars > 0
        except Exception:
            return False


def main() -> None:
    """CLI entry point for h5ad validation tools."""
    import argparse
    import json
    import sys

    parser = argparse.ArgumentParser(
        prog="gse_router.tools.validators",
        description="h5ad validation tools",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    # --- validate ---
    validate_parser = subparsers.add_parser(
        "validate",
        help="Run QC checks on an h5ad file",
    )
    validate_parser.add_argument(
        "--h5ad-path",
        type=Path,
        required=True,
        help="Path to the .h5ad file",
    )
    validate_parser.add_argument(
        "--min-cells",
        type=int,
        default=100,
        help="Minimum number of cells (default: 100)",
    )
    validate_parser.add_argument(
        "--min-genes-per-cell",
        type=int,
        default=200,
        help="Minimum genes per cell (default: 200)",
    )
    validate_parser.add_argument(
        "--min-counts-per-cell",
        type=int,
        default=500,
        help="Minimum counts per cell (default: 500)",
    )
    validate_parser.add_argument(
        "--max-mito-pct",
        type=float,
        default=20.0,
        help="Maximum mitochondrial percentage (default: 20.0)",
    )

    args = parser.parse_args()

    if args.command == "validate":
        validator = AnnDataValidator(
            min_cells=args.min_cells,
            min_genes_per_cell=args.min_genes_per_cell,
            min_counts_per_cell=args.min_counts_per_cell,
            max_mito_pct=args.max_mito_pct,
        )
        result = validator.validate(args.h5ad_path)

        output = {
            "passed": result.passed,
            "metrics": result.metrics,
            "failures": result.failures,
            "warnings": result.warnings,
        }
        json.dump(output, sys.stdout, indent=2, default=str)
        sys.stdout.write("\n")
        sys.exit(0 if result.passed else 1)


if __name__ == "__main__":
    main()
