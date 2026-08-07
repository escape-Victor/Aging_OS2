"""Metadata tools — merge external metadata into AnnData objects.

Handles:
- Cell-level metadata (obs) — barcode-to-metadata joins
- Sample-level metadata — sample-level annotations
- Gene annotations (var) — gene symbol/ID mapping
- GEO-derived study context — condition groups, treatments
"""

from __future__ import annotations

from pathlib import Path

import structlog

logger = structlog.get_logger(__name__)


class MetadataMergeResult:
    """Result of a metadata merge operation."""

    def __init__(
        self,
        success: bool,
        output_path: Path | None = None,
        n_cells_matched: int = 0,
        n_cells_total: int = 0,
        n_columns_added: int = 0,
        error: str | None = None,
        provenance: dict | None = None,
    ) -> None:
        self.success = success
        self.output_path = output_path
        self.n_cells_matched = n_cells_matched
        self.n_cells_total = n_cells_total
        self.n_columns_added = n_columns_added
        self.error = error
        self.provenance = provenance or {}


def _read_first_line(path: Path) -> str:
    """Read the first line of a file, handling gzip transparently."""
    if path.suffix.lower() == ".gz":
        import gzip
        with gzip.open(path, "rt", encoding="utf-8", errors="replace") as fh:
            return fh.readline()
    else:
        with open(path, encoding="utf-8", errors="replace") as fh:
            return fh.readline()


def _auto_detect_join_column(meta_df: "pd.DataFrame", obs_index: "pd.Index") -> str:
    """Auto-detect which column in metadata best matches the h5ad obs index.

    Tries common column names first (barcode, cell, cell_name, cellname),
    then falls back to finding the column with the highest overlap with obs_index.
    """
    # Priority-ordered candidate column names
    candidates = [
        "barcode", "Barcode", "BARCODE",
        "cell", "Cell", "CELL",
        "cell_name", "CellName", "cellname",
        "cell_id", "CellID", "cellid",
    ]
    for col in candidates:
        if col in meta_df.columns:
            return col

    # Fallback: find column with highest overlap with obs_index
    obs_set = set(str(i) for i in obs_index)
    best_col = None
    best_overlap = 0
    for col in meta_df.columns:
        col_vals = set(str(v) for v in meta_df[col])
        overlap = len(obs_set & col_vals)
        if overlap > best_overlap:
            best_overlap = overlap
            best_col = col
    if best_col and best_overlap > 0:
        return best_col
    return None  # No suitable column found — agent must decide


class MetadataTools:
    """Tools for merging metadata into AnnData objects.

    These are idempotent merge operations that take an existing h5ad
    and enrich its obs/var with external annotations.
    """

    def __init__(self) -> None:
        self._logger = logger.bind(component="MetadataTools")

    # ------------------------------------------------------------------
    # Cell metadata merge
    # ------------------------------------------------------------------

    def merge_cell_metadata(
        self,
        h5ad_path: Path,
        metadata_path: Path,
        join_column: str | None = None,
        delimiter: str | None = None,
        output_path: Path | None = None,
    ) -> MetadataMergeResult:
        """Merge cell-level metadata from a CSV/TSV file into adata.obs.

        Matches on barcode by default. Can join on any shared column.

        Args:
            h5ad_path: Path to the h5ad to enrich.
            metadata_path: Path to the CSV/TSV metadata file.
            join_column: Column in metadata to join on. If None, auto-detected
                from common names (barcode, cell, etc.) and overlap analysis.
            delimiter: Column delimiter (auto-detected if None).
            output_path: Where to write the enriched h5ad. Overwrites input if None.

        Returns:
            MetadataMergeResult with match statistics.
        """
        self._logger.info(
            "Merging cell metadata",
            h5ad=str(h5ad_path),
            metadata=str(metadata_path),
            join_column=join_column,
        )
        try:
            import anndata
            import pandas as pd

            adata = anndata.read_h5ad(h5ad_path)

            if delimiter is None:
                first_line = _read_first_line(metadata_path)
                delimiter = "\t" if "\t" in first_line else ","

            meta_df = pd.read_csv(metadata_path, sep=delimiter, index_col=None)

            # Auto-detect join column if the requested one doesn't exist
            if join_column not in meta_df.columns and join_column not in adata.obs.columns:
                detected = _auto_detect_join_column(meta_df, adata.obs.index)
                if detected:
                    join_column = detected
                    self._logger.info(
                        "Auto-detected join column",
                        requested=join_column,
                        detected=detected,
                    )
                else:
                    return MetadataMergeResult(
                        success=False,
                        error=(
                            f"No join column found in metadata. "
                            f"Available columns: {list(meta_df.columns)}. "
                            f"adata.obs columns: {list(adata.obs.columns)}. "
                            f"Tried common names: barcode, cell, cell_name, cell_id, "
                            f"and overlap-based detection found no match. "
                            f"Specify join_column explicitly."
                        ),
                    )

            # Determine if we join on index (barcodes as index) or column
            existing_cols = set(adata.obs.columns)
            if join_column in existing_cols:
                # Join on column
                merged = adata.obs.merge(
                    meta_df,
                    left_on=join_column,
                    right_on=join_column,
                    how="left",
                    suffixes=("", "_meta"),
                )
            else:
                # Join on index (barcode in index)
                merged = adata.obs.merge(
                    meta_df,
                    left_index=True,
                    right_on=join_column,
                    how="left",
                    suffixes=("", "_meta"),
                )

            n_new_cols = len(merged.columns) - len(adata.obs.columns)
            n_matched = merged.notna().any(axis=1).sum()

            # Replace obs with merged version (keeping index)
            merged.index = adata.obs.index
            adata.obs = merged

            out = output_path or h5ad_path
            adata.write_h5ad(out)

            return MetadataMergeResult(
                success=True,
                output_path=out,
                n_cells_matched=int(n_matched),
                n_cells_total=adata.n_obs,
                n_columns_added=n_new_cols,
                provenance={
                    "method": "merge_cell_metadata",
                    "metadata_file": str(metadata_path),
                    "join_column": join_column,
                    "new_columns": list(meta_df.columns),
                },
            )
        except Exception as exc:
            self._logger.error("Cell metadata merge failed", error=str(exc))
            # Include column context so the agent can decide next steps
            try:
                col_info = f" (metadata columns: {list(meta_df.columns)})"
            except Exception:
                col_info = ""
            return MetadataMergeResult(success=False, error=f"{exc}{col_info}")

    # ------------------------------------------------------------------
    # Sample metadata merge
    # ------------------------------------------------------------------

    def merge_sample_metadata(
        self,
        h5ad_path: Path,
        sample_metadata: dict[str, dict[str, str]],
        sample_column: str = "sample_id",
        output_path: Path | None = None,
    ) -> MetadataMergeResult:
        """Merge sample-level metadata into adata.obs.

        Takes a dict mapping sample_id -> {key: value, ...} and broadcasts
        values to all cells belonging to each sample.

        Args:
            h5ad_path: Path to the h5ad to enrich.
            sample_metadata: Dict of sample_id -> attribute dict.
            sample_column: Column in adata.obs that contains sample_id.
            output_path: Where to write the enriched h5ad.

        Returns:
            MetadataMergeResult with match statistics.
        """
        self._logger.info("Merging sample metadata", h5ad=str(h5ad_path), samples=list(sample_metadata.keys()))
        try:
            import anndata
            import pandas as pd

            adata = anndata.read_h5ad(h5ad_path)

            if sample_column not in adata.obs.columns:
                return MetadataMergeResult(
                    success=False,
                    error=f"Sample column '{sample_column}' not found in adata.obs",
                )

            meta_df = pd.DataFrame(sample_metadata).T
            meta_df.index.name = sample_column
            meta_df = meta_df.reset_index()

            merged = adata.obs.merge(meta_df, on=sample_column, how="left", suffixes=("", "_sample"))
            n_new_cols = len(merged.columns) - len(adata.obs.columns)
            n_matched = merged[[c for c in meta_df.columns if c != sample_column]].notna().any(axis=1).sum()

            merged.index = adata.obs.index
            adata.obs = merged

            out = output_path or h5ad_path
            adata.write_h5ad(out)

            return MetadataMergeResult(
                success=True,
                output_path=out,
                n_cells_matched=int(n_matched),
                n_cells_total=adata.n_obs,
                n_columns_added=n_new_cols,
                provenance={"method": "merge_sample_metadata", "samples": list(sample_metadata.keys())},
            )
        except Exception as exc:
            self._logger.error("Sample metadata merge failed", error=str(exc))
            return MetadataMergeResult(success=False, error=str(exc))

    # ------------------------------------------------------------------
    # Gene annotation merge
    # ------------------------------------------------------------------

    def merge_gene_annotation(
        self,
        h5ad_path: Path,
        annotation_path: Path,
        join_column: str = "gene_id",
        output_path: Path | None = None,
    ) -> MetadataMergeResult:
        """Merge gene-level annotations into adata.var.

        Args:
            h5ad_path: Path to the h5ad to enrich.
            annotation_path: Path to the gene annotation file (CSV/TSV/GTF).
            join_column: Column in annotation to join on.
            output_path: Where to write the enriched h5ad.

        Returns:
            MetadataMergeResult with match statistics.
        """
        self._logger.info("Merging gene annotation", h5ad=str(h5ad_path), annotation=str(annotation_path))
        try:
            import anndata
            import pandas as pd

            adata = anndata.read_h5ad(h5ad_path)

            if annotation_path.suffix.lower() in (".gtf", ".gff", ".gff3", ".bed"):
                # Parse GTF/GFF (simplified — full parsing would be more robust)
                anno_df = self._parse_gtf_simple(annotation_path)
            else:
                anno_df = pd.read_csv(annotation_path, sep=None, engine="python")

            merged = adata.var.merge(
                anno_df,
                left_index=True,
                right_on=join_column,
                how="left",
                suffixes=("", "_anno"),
            )
            n_new_cols = len(merged.columns) - len(adata.var.columns)
            merged.index = adata.var.index
            adata.var = merged

            out = output_path or h5ad_path
            adata.write_h5ad(out)

            return MetadataMergeResult(
                success=True,
                output_path=out,
                n_cells_matched=adata.n_obs,
                n_cells_total=adata.n_obs,
                n_columns_added=n_new_cols,
                provenance={"method": "merge_gene_annotation", "annotation_file": str(annotation_path)},
            )
        except Exception as exc:
            self._logger.error("Gene annotation merge failed", error=str(exc))
            return MetadataMergeResult(success=False, error=str(exc))

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _parse_gtf_simple(gtf_path: Path):
        """Parse a GTF file into a simple gene annotation DataFrame (returns pandas.DataFrame)."""
        import pandas as pd  # noqa: F811

        columns = ["seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"]
        df = pd.read_csv(gtf_path, sep="\t", comment="#", names=columns, header=None)

        # Extract gene_id and gene_name from attributes column
        def _extract_attr(attr_str: str) -> dict[str, str]:
            attrs = {}
            for part in attr_str.split(";"):
                part = part.strip()
                if " " in part:
                    key, val = part.split(" ", 1)
                    attrs[key.strip()] = val.strip().strip('"')
            return attrs

        attrs_df = df["attribute"].apply(_extract_attr).apply(pd.Series)
        result = pd.concat([df, attrs_df], axis=1)

        # Keep only gene features
        if "feature" in result.columns:
            result = result[result["feature"] == "gene"]

        return result.drop_duplicates(subset="gene_id") if "gene_id" in result.columns else result


# ------------------------------------------------------------------
# CLI entry point
# ------------------------------------------------------------------


def main() -> None:
    import argparse
    import json
    import sys

    parser = argparse.ArgumentParser(description="Metadata merge tools for AnnData objects")
    subparsers = parser.add_subparsers(dest="command", required=True)

    # merge-cells subcommand
    cells_parser = subparsers.add_parser("merge-cells", help="Merge cell-level metadata into AnnData.obs")
    cells_parser.add_argument("--h5ad-path", required=True, type=Path, help="Path to the h5ad file")
    cells_parser.add_argument("--metadata-path", required=True, type=Path, help="Path to the metadata CSV/TSV file")
    cells_parser.add_argument("--join-column", default="barcode", help="Column to join on (default: barcode)")
    cells_parser.add_argument("--output-path", type=Path, default=None, help="Output path (overwrites input if not set)")

    # merge-genes subcommand
    genes_parser = subparsers.add_parser("merge-genes", help="Merge gene annotation into AnnData.var")
    genes_parser.add_argument("--h5ad-path", required=True, type=Path, help="Path to the h5ad file")
    genes_parser.add_argument("--metadata-path", required=True, type=Path, help="Path to the gene annotation file")
    genes_parser.add_argument("--join-column", default="gene_id", help="Column to join on (default: gene_id)")
    genes_parser.add_argument("--output-path", type=Path, default=None, help="Output path (overwrites input if not set)")

    args = parser.parse_args()
    tools = MetadataTools()

    if args.command == "merge-cells":
        result = tools.merge_cell_metadata(
            h5ad_path=args.h5ad_path,
            metadata_path=args.metadata_path,
            join_column=args.join_column,
            output_path=args.output_path,
        )
    elif args.command == "merge-genes":
        result = tools.merge_gene_annotation(
            h5ad_path=args.h5ad_path,
            annotation_path=args.metadata_path,
            join_column=args.join_column,
            output_path=args.output_path,
        )
    else:
        parser.print_help()
        sys.exit(1)

    output = {
        "success": result.success,
        "output_path": str(result.output_path) if result.output_path else None,
        "n_cells_matched": result.n_cells_matched,
        "n_cells_total": result.n_cells_total,
        "n_columns_added": result.n_columns_added,
        "error": result.error,
        "provenance": result.provenance,
    }
    print(json.dumps(output, indent=2))
    sys.exit(0 if result.success else 1)


if __name__ == "__main__":
    main()
