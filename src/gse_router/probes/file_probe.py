"""Local file probe — inspects downloaded files and writes per-file JSON sidecars.

Scans directories, decompresses archives, fingerprints content, and detects:
FASTQ, 10x MEX, 10x H5, h5ad, loom, RDS, CSV/TSV matrices, metadata, BAM, etc.
"""

from __future__ import annotations

import gzip as _gzip_module
import hashlib
import re as _re
import tarfile
import zipfile
from pathlib import Path

import structlog

from ..schemas.sidecar import (
    FileDetectedType,
    FileSidecar,
    ProcessingStage,
    RouteRecommendation,
)

logger = structlog.get_logger(__name__)

# Signature-based detection constants
MATRIX_MARKET_HEADER = b"%%MatrixMarket"
HDF5_SIGNATURE = b"\x89HDF\r\n\x1a\n"
LOOM_SIGNATURES = [b"LOOM", b"/matrix", b"/layers"]  # Loom is HDF5 with specific groups
FASTQ_QUAL_HEADER = b"@"
GZIP_MAGIC = b"\x1f\x8b"
ZLIB_MAGIC = b"\x78\x9c"
RDS_MAGIC = b"RDX2\n"  # R serialized object
BAM_MAGIC = b"BAM\x01"
CRAM_MAGIC = b"CRAM"


class FileProbe:
    """Inspects local files and writes per-file JSON sidecars.

    This is the second intelligence layer: it examines every file in the
    dataset, detects its format by content (not filename), and documents
    findings in structured sidecars so the route planner never reasons
    from raw directory listings.
    """

    # File patterns for common single-cell formats
    MEX_MATRIX_PATTERNS = ["matrix.mtx", "matrix.mtx.gz"]
    MEX_BARCODES_PATTERNS = ["barcodes.tsv", "barcodes.tsv.gz"]
    MEX_FEATURES_PATTERNS = ["features.tsv", "features.tsv.gz", "genes.tsv", "genes.tsv.gz"]
    H5_PATTERNS = ["filtered_feature_bc_matrix.h5", "raw_feature_bc_matrix.h5"]

    def __init__(self, max_file_size_mb: int = 5000) -> None:
        self.max_file_size_mb = max_file_size_mb
        self._sha256_cache: dict[str, str] = {}
        self._logger = logger.bind(component="FileProbe")

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    @staticmethod
    def _is_in_extracted_dir(file_path: Path, root: Path) -> bool:
        """Check whether *file_path* lives inside an ``_extracted_*`` directory.

        These directories are artifacts from prior archive extraction runs.
        They should be excluded from the initial probe pass because:
        - Their contents will be re-extracted and probed fresh during execution
        - Including them causes duplicate sidecars that confuse the route planner
          (both MEX_10X from extracted files AND ARCHIVE_10X from the archive)
        """
        try:
            rel = file_path.relative_to(root)
        except ValueError:
            return False
        for part in rel.parts:
            if part.startswith("_extracted_"):
                return True
        return False

    def probe_directory(self, directory: Path) -> list[FileSidecar]:
        """Recursively probe all files in a directory.

        Entry point for scanning a downloaded dataset. Decompresses archives
        in-place and probes extracted files as well.

        Args:
            directory: Root directory of the dataset.

        Returns:
            List of FileSidecar, one per file found.
        """
        self._logger.info("Probing directory", path=str(directory))
        sidecars: list[FileSidecar] = []
        for file_path in directory.rglob("*"):
            # Skip files inside _extracted_* directories — these are artifacts
            # from prior extraction runs. Their contents will be probed fresh
            # when the archive is extracted again during execution.
            if self._is_in_extracted_dir(file_path, directory):
                continue
            if file_path.is_file():
                sc = self.probe_file(file_path, relative_root=directory)
                if sc:
                    sidecars.append(sc)
                    # If this is an archive, extract and probe contents
                    if sc.detected_type == FileDetectedType.ARCHIVE:
                        extracted = self._extract_and_probe(file_path, sc.file_id, directory)
                        sidecars.extend(extracted)
                    # If this is a gzip-compressed data file (csv/tsv/txt), decompress
                    # and inspect the actual content for better classification
                    elif (sc.extension in (".gz",) and sc.detected_type in (
                            FileDetectedType.CSV_MATRIX,
                            FileDetectedType.TSV_MATRIX,
                            FileDetectedType.TXT_MATRIX,
                            FileDetectedType.CELL_METADATA,
                            FileDetectedType.UNKNOWN,
                    )):
                        enriched = self._decompress_gzip_and_reprobe(
                            file_path, sc, directory
                        )
                        if enriched:
                            # Replace the original sidecar with the enriched one
                            for i, s in enumerate(sidecars):
                                if s.file_id == sc.file_id:
                                    sidecars[i] = enriched
                                    break
        return sidecars

    def probe_file(self, file_path: Path, relative_root: Path | None = None) -> FileSidecar | None:
        """Probe a single file and return its sidecar.

        Args:
            file_path: Absolute path to the file.
            relative_root: Root for computing relative path. Defaults to parent.

        Returns:
            FileSidecar, or None if the file should be skipped (too large, unreadable).
        """
        if not file_path.is_file():
            return None

        size = file_path.stat().st_size
        if size > self.max_file_size_mb * 1024 * 1024:
            self._logger.warning("Skipping large file", path=str(file_path), size_mb=size / 1e6)
            return None

        root = relative_root or file_path.parent
        rel_path = str(file_path.relative_to(root))

        sha256 = self._compute_sha256(file_path)
        detected_type, confidence, evidence = self._detect_format(file_path)

        sidecar = FileSidecar(
            file_id=f"sha256:{sha256}",
            path=rel_path,
            filename=file_path.name,
            size_bytes=size,
            extension=file_path.suffix.lower(),
            detected_type=detected_type,
            processing_stage=self._infer_stage(detected_type),
            confidence=confidence,
            evidence=evidence,
        )

        # Add recommended routes
        sidecar.recommended_routes = self._recommend_routes(sidecar)

        self._logger.debug(
            "Probed file",
            path=rel_path,
            type=detected_type.value,
            confidence=confidence,
        )
        return sidecar

    # ------------------------------------------------------------------
    # Format detection
    # ------------------------------------------------------------------

    def _detect_format(self, file_path: Path) -> tuple[FileDetectedType, float, list[str]]:
        """Detect file format by reading header bytes and checking patterns.

        Returns:
            (detected_type, confidence, evidence_list)
        """
        try:
            with open(file_path, "rb") as fh:
                header = fh.read(4096)
        except (OSError, PermissionError) as exc:
            return FileDetectedType.UNKNOWN, 0.0, [f"Cannot read file: {exc}"]

        if not header:
            return FileDetectedType.UNKNOWN, 0.3, ["Empty file"]

        # GEO supplementary manifests — not biological data
        if file_path.name.lower() in ("filelist.txt", "filelist.txt.gz", "checksum_md5.txt"):
            return FileDetectedType.METADATA, 0.95, [
                "GEO supplementary checksum manifest — exclude from processing"
            ]

        # Check against signatures in priority order
        detector = self._match_signature(file_path, header)
        if detector:
            return detector

        # Fallback: detect by filename pattern
        detector = self._match_filename(file_path)
        if detector:
            # For CSV/TSV/TXT files, the filename alone is ambiguous —
            # prefer text content inspection for more accurate classification.
            dtype, conf, ev = detector
            if dtype in (FileDetectedType.CSV_MATRIX, FileDetectedType.TSV_MATRIX, FileDetectedType.TXT_MATRIX):
                if self._looks_like_text(header):
                    return self._inspect_text_content(file_path, header)
            return detector

        # If nothing matched but it looks like text, try deeper inspection
        if self._looks_like_text(header):
            return self._inspect_text_content(file_path, header)

        return FileDetectedType.UNKNOWN, 0.2, [f"Unrecognized format, header: {header[:20]!r}"]

    def _match_signature(self, file_path: Path, header: bytes) -> tuple[FileDetectedType, float, list[str]] | None:
        """Match file content against known format signatures."""
        ext = file_path.suffix.lower()

        # MatrixMarket
        if header[:14] == MATRIX_MARKET_HEADER:
            return FileDetectedType.MEX_MATRIX_10X, 0.96, ["MatrixMarket header detected"]

        # HDF5 — needs deeper inspection to distinguish 10x H5 from h5ad from loom
        if header[:8] == HDF5_SIGNATURE:
            return self._classify_hdf5(file_path)

        # Gzip-compressed files
        if header[:2] == GZIP_MAGIC:
            return self._classify_gzip(file_path)

        # RDS (R serialized object)
        if header[:5] == RDS_MAGIC:
            return FileDetectedType.RDS_SEURAT, 0.85, ["RDS magic bytes detected"]
        if header[:2] == b"R\n":
            return FileDetectedType.RDS_SEURAT, 0.60, ["Possible R serialized object"]

        # BAM/CRAM
        if header[:4] == BAM_MAGIC:
            return FileDetectedType.BAM, 0.99, ["BAM magic bytes detected"]
        if header[:4] == CRAM_MAGIC:
            return FileDetectedType.CRAM, 0.99, ["CRAM magic bytes detected"]

        # FASTQ (starts with @ followed by printable ASCII)
        if ext in (".fastq", ".fq") and len(header) > 0 and header[0] == 0x40:
            return FileDetectedType.FASTQ, 0.85, ["FASTQ extension + @ start character"]

        # ZIP
        if header[:4] == b"PK\x03\x04":
            evidence = ["ZIP signature detected"]
            # Peek at archive contents for better classification
            contents = self._peek_archive_contents(file_path)
            if contents:
                _, _, content_evidence = self._classify_archive_by_contents(contents, file_path)
                evidence.extend(content_evidence)
            return FileDetectedType.ARCHIVE, 0.99, evidence
        # TAR
        if header[257:262] == b"ustar":
            evidence = ["TAR signature detected"]
            # Peek at archive contents for better classification
            contents = self._peek_archive_contents(file_path)
            if contents:
                _, _, content_evidence = self._classify_archive_by_contents(contents, file_path)
                evidence.extend(content_evidence)
            return FileDetectedType.ARCHIVE, 0.99, evidence

        return None

    def _classify_hdf5(self, file_path: Path) -> tuple[FileDetectedType, float, list[str]]:
        """Classify an HDF5 file by inspecting its internal groups with h5py.

        Content-first approach: opens the file with h5py and checks for
        the canonical group structure of each format. Filename/extension
        only adjusts confidence ±0.03 — it is never the primary signal.
        """
        ext = file_path.suffix.lower()
        name = file_path.name.lower()

        # Inspect internal HDF5 groups for definitive classification
        groups = self._inspect_hdf5_groups(file_path)

        if groups is not None:
            has_X = groups.get("X", False)
            has_obs = groups.get("obs", False)
            has_matrix = groups.get("matrix", False)
            has_layers = groups.get("layers", False)

            # h5ad: canonical structure has /X and /obs at minimum
            if has_X and has_obs:
                confidence = 0.95
                evidence = ["HDF5 internal groups: /X + /obs → AnnData (h5ad)"]
                if ext == ".h5ad":
                    confidence = 0.97
                    evidence.append("Extension .h5ad confirms")
                return FileDetectedType.H5AD, confidence, evidence

            # Loom: has /matrix and /layers (and typically /row_attrs, /col_attrs)
            if has_matrix and has_layers:
                confidence = 0.93
                evidence = ["HDF5 internal groups: /matrix + /layers → Loom"]
                if ext == ".loom":
                    confidence = 0.95
                    evidence.append("Extension .loom confirms")
                return FileDetectedType.LOOM, confidence, evidence

            # 10x H5: has /matrix but NOT /layers
            if has_matrix and not has_layers:
                confidence = 0.90
                evidence = ["HDF5 internal groups: /matrix (no /layers) → 10x H5"]
                if ext in (".h5", ".hdf5") and any(
                    p in name for p in ["feature", "bc_matrix", "filtered", "raw"]
                ):
                    confidence = 0.93
                    evidence.append("Filename consistent with 10x output")
                return FileDetectedType.H5_10X, confidence, evidence

            # HDF5 with unrecognized internal structure — do NOT guess
            return FileDetectedType.UNKNOWN, 0.40, [
                f"HDF5 with unrecognized internal structure. "
                f"Groups found: {sorted(groups.keys())}"
            ]

        # h5py failed to open — fall back to extension-based classification
        # with sharply reduced confidence
        if ext == ".h5ad":
            return FileDetectedType.H5AD, 0.60, [
                "HDF5 signature + .h5ad extension (internal inspection failed)"
            ]
        if ext in (".h5seurat",):
            return FileDetectedType.H5SEURAT, 0.55, [
                "HDF5 signature + .h5seurat extension (internal inspection failed)"
            ]
        if ext in (".loom",):
            return FileDetectedType.LOOM, 0.60, [
                "HDF5 signature + .loom extension (internal inspection failed)"
            ]
        if ext in (".h5", ".hdf5"):
            if any(p in name for p in ["feature", "bc_matrix", "filtered", "raw"]):
                return FileDetectedType.H5_10X, 0.55, [
                    "HDF5 signature + 10x-like filename (internal inspection failed)"
                ]

        return FileDetectedType.UNKNOWN, 0.35, [
            "HDF5 signature — could not inspect internal structure, ambiguous subtype"
        ]

    @staticmethod
    def _inspect_hdf5_groups(file_path: Path) -> dict[str, bool] | None:
        """Inspect the top-level groups of an HDF5 file using h5py.

        Returns a dict mapping group name → present (True), or None if the
        file cannot be opened with h5py.
        """
        try:
            import h5py
        except ImportError:
            return None

        try:
            with h5py.File(str(file_path), "r") as f:
                groups: dict[str, bool] = {}
                # Check for canonical single-cell HDF5 groups
                for key in ("X", "obs", "var", "matrix", "layers", "raw"):
                    groups[key] = key in f
                # Also collect any other top-level groups for evidence
                for key in f.keys():
                    if key not in groups:
                        groups[key] = True
                return groups
        except Exception:
            return None

    def _classify_gzip(self, file_path: Path) -> tuple[FileDetectedType, float, list[str]]:
        """Classify a gzip-compressed file — content-first, filename as tiebreaker.

        Three ordered stages:
        1. Content-first detection — decompress header, inspect structure, sample values
        2. Specialized non-count detection — BED/peak/GFF/V(D)J structure by content
        3. Filename hints — only when content is ambiguous, with reduced confidence
        """
        name = file_path.name.lower().replace(".gz", "")
        full_name = file_path.name.lower()
        inner_ext = Path(name).suffix.lower() if name else ""

        # ==================================================================
        # STAGE 1: Content-first detection (runs FIRST, before any filename check)
        # ==================================================================

        # 1a. Check if this is actually a compressed archive (tar.gz or zip within gz)
        if self._is_gzip_wrapping_archive(file_path):
            # Peek inside the archive to classify its contents
            archive_files = self._peek_archive_contents(file_path)
            if archive_files:
                dtype, conf, ev = self._classify_archive_by_contents(archive_files, file_path)
                return dtype, conf, ["Gzip archive, content-peeked: " + e for e in ev]

        # 1b. Decompress and inspect content for MEX components / matrices
        content_result = self._inspect_gzip_mex_content(file_path)
        if content_result:
            return content_result

        # 1c. For .txt/.tsv/.csv files, do deeper content inspection
        if name.endswith((".txt", ".tsv", ".csv")) or inner_ext in (".txt", ".tsv", ".csv"):
            deeper = self._inspect_gzip_text_content(file_path)
            if deeper:
                return deeper

        # ==================================================================
        # STAGE 2: Specialized non-count format detection (by content structure)
        # ==================================================================

        # 2a. GFF/GTF gene annotations — check for annotation column structure
        if inner_ext in (".gff", ".gff3", ".gtf") or name.endswith((".gff", ".gff3", ".gtf")):
            if self._verify_gff_content(file_path):
                return FileDetectedType.GENE_ANNOTATION, 0.90, [
                    "Gzip + GFF/GTF annotation structure verified by content"
                ]

        # 2b. Peak data — check for BED/narrowPeak column structure
        if inner_ext in (".bed", ".narrowpeak", ".broadpeak"):
            if self._verify_peak_content(file_path):
                return FileDetectedType.PEAK_DATA, 0.90, [
                    "Gzip + peak file structure verified by content"
                ]

        # 2c. Spatial transcriptomics — check for image coordinate data
        spatial_kw = [
            "spatial", "tissue_positions", "scalefactors",
            "tissue_hires_image", "tissue_lowres_image",
            "aligned_fiducials", "detected_tissue_image",
        ]
        if any(kw in name for kw in spatial_kw):
            spatial_verified = self._verify_spatial_content(file_path)
            if spatial_verified:
                return FileDetectedType.UNKNOWN, 0.88, [
                    "Gzip + spatial transcriptomics structure verified by content"
                ]

        # ==================================================================
        # STAGE 3: Filename hints — tiebreaker only, reduced confidence
        # ==================================================================

        # 3a. MEX components — filename confirmation only (content already checked above)
        if "matrix.mtx" in name:
            return FileDetectedType.MEX_MATRIX_10X, 0.75, [
                "Gzip + MatrixMarket filename (content ambiguous, filename suggests matrix)"
            ]
        if "barcodes.tsv" in name:
            return FileDetectedType.MEX_BARCODES_10X, 0.75, [
                "Gzip + barcodes filename (content ambiguous, filename suggests barcodes)"
            ]
        if "features.tsv" in name or "genes.tsv" in name:
            return FileDetectedType.MEX_FEATURES_10X, 0.75, [
                "Gzip + features/genes filename (content ambiguous, filename suggests features)"
            ]

        # 3b. Non-count filename patterns (with content verification)
        non_count_result = self._classify_gzip_noncount_by_filename(name, inner_ext)
        if non_count_result:
            return non_count_result

        # 3c. FASTQ
        if name.endswith(".fastq"):
            return FileDetectedType.FASTQ, 0.65, ["Gzip + FASTQ filename (content not verified)"]

        # 3d. MATLAB .mat
        if name.endswith(".mat") or inner_ext == ".mat":
            return FileDetectedType.MATRIX_MAT, 0.60, ["Gzip + MATLAB .mat filename (content not verified)"]

        # 3e. Generic delimited matrix files — lowest confidence
        if name.endswith(".csv"):
            return FileDetectedType.CSV_MATRIX, 0.55, ["Gzip + CSV extension (content ambiguous)"]
        if name.endswith(".tsv"):
            return FileDetectedType.TSV_MATRIX, 0.55, ["Gzip + TSV extension (content ambiguous)"]
        if name.endswith(".txt"):
            return FileDetectedType.TXT_MATRIX, 0.50, ["Gzip + TXT extension (content ambiguous)"]

        return FileDetectedType.UNKNOWN, 0.35, ["Gzip-compressed, unknown content"]

    def _match_filename(self, file_path: Path) -> tuple[FileDetectedType, float, list[str]] | None:
        """Detect format by filename pattern — low-confidence fallback only.

        Filename-based detection is the weakest signal. All matches are capped
        at confidence 0.65 and carry a "content not verified" warning.
        Content-based detection (magic bytes, HDF5 groups, text inspection)
        is always preferred and runs before this method is called.
        """
        name = file_path.name.lower()
        ext = file_path.suffix.lower()

        if ext == ".h5ad":
            return FileDetectedType.H5AD, 0.60, [
                "Extension .h5ad — filename-based only, content not verified"
            ]
        if ext == ".loom":
            return FileDetectedType.LOOM, 0.60, [
                "Extension .loom — filename-based only, content not verified"
            ]
        if ext in (".rds", ".rdata", ".rda"):
            if any(kw in name for kw in ["seurat", "seu"]):
                return FileDetectedType.RDS_SEURAT, 0.58, [
                    "Extension + Seurat keyword — filename-based only, content not verified"
                ]
            if any(kw in name for kw in ["sce", "singlecell", "single_cell_experiment"]):
                return FileDetectedType.RDS_SCE, 0.58, [
                    "Extension + SCE keyword — filename-based only, content not verified"
                ]
            return FileDetectedType.RDS_SEURAT, 0.35, [
                "Extension .rds (ambiguous R object) — content not verified"
            ]
        if ext == ".h5seurat":
            return FileDetectedType.H5SEURAT, 0.60, [
                "Extension .h5seurat — filename-based only, content not verified"
            ]

        # --- Non-count supplementary file detection ---
        # Filename-based quick-reject filters. Content must be verified in
        # _classify_gzip / _inspect_text_content to confirm; filename alone
        # is insufficient for high-confidence exclusion.

        # Immune receptor data (BCR/TCR contig annotations)
        if any(kw in name for kw in [
            "filtered_contig_annotations_bcr",
            "filtered_contig_annotations_tcr",
            "contig_annotations",
            "vdj",
            "clonotypes",
        ]):
            return FileDetectedType.IMMUNE_RECEPTOR_DATA, 0.55, [
                "Filename suggests BCR/TCR immune receptor data — content not verified"
            ]

        # Peak data (ATAC-seq, ChIP-seq)
        if ext in (".bed", ".bed.gz", ".narrowpeak", ".broadpeak"):
            return FileDetectedType.PEAK_DATA, 0.62, [
                "Peak file extension — content not verified"
            ]
        if any(kw in name for kw in [
            "_peaks_", "_peak.", "chip_peaks", "atac_peaks",
            "narrowpeak", "broadpeak", "macs2_peaks",
        ]):
            return FileDetectedType.PEAK_DATA, 0.52, [
                "Filename suggests peak data — content not verified"
            ]

        # Differential expression/results
        if any(kw in name for kw in [
            "deseq2", "diffbind", "diffbind_result",
            "voom_weights", "limma_result", "edger_result",
            "differential_expression", "deg_",
        ]):
            return FileDetectedType.DIFF_RESULTS, 0.52, [
                "Filename suggests differential expression results — content not verified"
            ]

        # Spatial transcriptomics data
        if any(kw in name for kw in [
            "spatial", "tissue_positions_list", "scalefactors_json",
            "tissue_hires_image", "tissue_lowres_image",
            "aligned_fiducials", "detected_tissue_image",
        ]):
            return FileDetectedType.UNKNOWN, 0.55, [
                "Filename suggests spatial transcriptomics data — content not verified"
            ]

        # Gene annotation extensions
        if ext in (".gtf", ".gff", ".gff3"):
            return FileDetectedType.GENE_ANNOTATION, 0.65, [
                "Annotation file extension — filename-based only"
            ]

        # Other specialized formats
        if ext in (".fastq", ".fq", ".fastq.gz", ".fq.gz"):
            return FileDetectedType.FASTQ, 0.55, ["FASTQ extension — content not verified"]
        if ext in (".sra",):
            return FileDetectedType.SRA_ARCHIVE, 0.65, ["SRA extension"]
        if ext in (".mat",):
            return FileDetectedType.MATRIX_MAT, 0.55, [
                "MATLAB .mat file — may contain count matrix, content not verified"
            ]
        if ext in (".bam",):
            return FileDetectedType.BAM, 0.65, ["BAM extension — content not verified"]
        if ext in (".cram",):
            return FileDetectedType.CRAM, 0.65, ["CRAM extension — content not verified"]

        # Generic delimited text — last resort (AFTER all specific checks)
        if ext in (".csv", ".csv.gz"):
            return FileDetectedType.CSV_MATRIX, 0.50, [
                "Extension .csv — filename-based only, content not verified"
            ]
        if ext in (".tsv", ".tsv.gz", ".txt", ".txt.gz"):
            return FileDetectedType.TSV_MATRIX, 0.45, [
                "Delimited text extension — filename-based only, content not verified"
            ]

    def _inspect_text_content(
        self, file_path: Path, header: bytes
    ) -> tuple[FileDetectedType, float, list[str]]:
        """Inspect text file content for matrix or metadata patterns."""
        try:
            text = header.decode("utf-8", errors="replace")[:4096]
        except UnicodeDecodeError:
            return FileDetectedType.UNKNOWN, 0.1, ["Binary content"]

        lines = text.split("\n")
        if not lines:
            return FileDetectedType.UNKNOWN, 0.1, ["Empty text file"]

        first_line = lines[0].strip()

        # Comma-separated with gene-like first column
        if "\t" in first_line or "," in first_line:
            sep = "\t" if "\t" in first_line else ","
            cols = first_line.split(sep)
            if len(cols) > 10:
                # Sample data values to determine raw vs normalized
                data_stats = self._sample_data_values(file_path)
                is_raw, raw_evidence = self._classify_raw_vs_normalized(data_stats)

                dtype = FileDetectedType.TSV_MATRIX if sep == "\t" else FileDetectedType.CSV_MATRIX
                evidence = [f"Matrix: {len(cols)} columns detected"] + raw_evidence

                if is_raw is True:
                    # Upgrade to confirmed raw UMI matrix
                    dtype = FileDetectedType.RAW_UMI_MATRIX
                    confidence = 0.82
                elif is_raw is False:
                    # Downgrade to normalized matrix
                    dtype = FileDetectedType.NORMALIZED_MATRIX
                    confidence = 0.80
                else:
                    confidence = 0.75

                return dtype, confidence, evidence
            # Check for metadata headers
            if any(
                kw in first_line.lower()
                for kw in ["cell", "barcode", "sample", "condition", "cluster", "cell_type"]
            ):
                return FileDetectedType.CELL_METADATA, 0.78, ["Cell metadata headers detected"]
            if any(
                kw in first_line.lower()
                for kw in ["gene", "symbol", "ensembl", "feature"]
            ):
                return FileDetectedType.GENE_ANNOTATION, 0.72, ["Gene annotation headers detected"]

            # Fallback: few columns (2-10) with data values → likely metadata
            if 2 <= len(cols) <= 10:
                return FileDetectedType.CELL_METADATA, 0.55, [
                    f"Few columns ({len(cols)}) — likely cell/sample metadata"
                ]

        return FileDetectedType.UNKNOWN, 0.2, ["Text file, unknown structure"]

    # ------------------------------------------------------------------
    # Stage inference & route recommendation
    # ------------------------------------------------------------------

    def _infer_stage(self, detected_type: FileDetectedType) -> ProcessingStage:
        """Infer the processing stage from the detected file type."""
        raw_types = {FileDetectedType.FASTQ, FileDetectedType.SRA_ARCHIVE}
        count_types = {
            FileDetectedType.MEX_MATRIX_10X,
            FileDetectedType.H5_10X,
            FileDetectedType.CSV_MATRIX,
            FileDetectedType.TSV_MATRIX,
            FileDetectedType.TXT_MATRIX,
            FileDetectedType.RAW_UMI_MATRIX,
            FileDetectedType.MATRIX_MAT,
            FileDetectedType.BAM,
            FileDetectedType.CRAM,
        }
        obj_types = {FileDetectedType.H5AD, FileDetectedType.LOOM,
                     FileDetectedType.RDS_SEURAT, FileDetectedType.RDS_SCE,
                     FileDetectedType.H5SEURAT}
        meta_types = {FileDetectedType.CELL_METADATA, FileDetectedType.SAMPLE_METADATA}
        anno_types = {FileDetectedType.GENE_ANNOTATION}
        non_count_types = {
            FileDetectedType.IMMUNE_RECEPTOR_DATA,
            FileDetectedType.PEAK_DATA,
            FileDetectedType.DIFF_RESULTS,
        }

        if detected_type in raw_types:
            return ProcessingStage.RAW_SEQUENCING
        if detected_type in count_types:
            return ProcessingStage.RAW_COUNTS
        if detected_type in obj_types:
            return ProcessingStage.PROCESSED_OBJECT
        if detected_type in meta_types:
            return ProcessingStage.METADATA
        if detected_type in anno_types:
            return ProcessingStage.ANNOTATION
        if detected_type in non_count_types:
            return ProcessingStage.UNKNOWN  # Not scRNA-seq count data — skip
        return ProcessingStage.UNKNOWN

    def _recommend_routes(self, sidecar: FileSidecar) -> list[RouteRecommendation]:
        """Recommend processing routes based on detected file type."""
        routes: dict[FileDetectedType, list[RouteRecommendation]] = {
            FileDetectedType.FASTQ: [
                RouteRecommendation(route="cellranger_count", priority=1,
                                    description="Cell Ranger count (10x)"),
                RouteRecommendation(route="starsolo", priority=2,
                                    description="STARsolo alignment"),
                RouteRecommendation(route="kallisto_bustools", priority=3,
                                    description="kallisto/bustools pseudoalignment"),
            ],
            FileDetectedType.MEX_MATRIX_10X: [
                RouteRecommendation(route="read_10x_mtx_to_h5ad", priority=1,
                                    description="Read 10x MEX with Scanpy"),
            ],
            FileDetectedType.H5_10X: [
                RouteRecommendation(route="read_10x_h5_to_h5ad", priority=1,
                                    description="Read 10x H5 with Scanpy"),
            ],
            FileDetectedType.H5AD: [
                RouteRecommendation(route="read_h5ad_validate", priority=1,
                                    description="Validate and standardize h5ad"),
            ],
            FileDetectedType.LOOM: [
                RouteRecommendation(route="read_loom_to_h5ad", priority=1,
                                    description="Convert Loom to h5ad"),
            ],
            FileDetectedType.RDS_SEURAT: [
                RouteRecommendation(route="read_rds_seurat_to_h5ad", priority=1,
                                    description="Convert Seurat RDS to h5ad"),
            ],
            FileDetectedType.RDS_SCE: [
                RouteRecommendation(route="read_rds_sce_to_h5ad", priority=1,
                                    description="Convert SCE RDS to h5ad"),
            ],
            FileDetectedType.CSV_MATRIX: [
                RouteRecommendation(route="read_csv_matrix_to_h5ad", priority=1,
                                    description="Parse CSV count matrix to h5ad"),
            ],
            FileDetectedType.TSV_MATRIX: [
                RouteRecommendation(route="read_csv_matrix_to_h5ad", priority=1,
                                    description="Parse TSV count matrix to h5ad"),
            ],
        }
        return routes.get(sidecar.detected_type, [
            RouteRecommendation(route="manual_review", priority=1,
                                description="Unknown format — requires manual review"),
        ])

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------

    def _compute_sha256(self, file_path: Path) -> str:
        """Compute SHA-256 hash of a file."""
        sha256 = hashlib.sha256()
        with open(file_path, "rb") as fh:
            for chunk in iter(lambda: fh.read(8192), b""):
                sha256.update(chunk)
        return sha256.hexdigest()

    @staticmethod
    def _inspect_gzip_mex_content(
        file_path: Path,
    ) -> tuple[FileDetectedType, float, list[str]] | None:
        """Inspect a gzip-compressed text file for MEX component patterns.

        Handles non-standard MEX filenames (e.g., *.cellBarcodes.txt.gz,
        *.genes.tsv.gz, *.dem.txt.gz) by inspecting the actual content.

        Returns None if the content doesn't look like a MEX component.
        """
        import gzip

        try:
            with gzip.open(file_path, "rt", encoding="utf-8", errors="replace") as fh:
                header = fh.read(8192)
        except Exception:
            return None

        if not header:
            return None

        lines = header.split("\n")
        if not lines:
            return None

        first_line = lines[0].strip()

        # MatrixMarket header → definitely a matrix
        if first_line.startswith("%%MatrixMarket"):
            return FileDetectedType.MEX_MATRIX_10X, 0.94, [
                "Gzip + MatrixMarket header (non-standard filename)"
            ]

        # Single-column file: could be barcodes
        # Count the number of tab-separated columns in the first few lines
        n_cols_first = len(first_line.split("\t"))
        if n_cols_first == 1:
            # Check a few more lines for single-column consistency
            single_col = all(
                len(line.split("\t")) == 1
                for line in lines[1:6] if line.strip()
            )
            if single_col and len(lines) > 1 and lines[1].strip():
                # Check if it looks like barcodes (ACGT sequences or long alphanumeric)
                sample = [line.strip() for line in lines[1:6] if line.strip()]
                barcode_like = sum(
                    1 for s in sample
                    if len(s) >= 14 and all(c in "ACGTacgt" for c in s)
                )
                if barcode_like >= len(sample) * 0.6:
                    return FileDetectedType.MEX_BARCODES_10X, 0.78, [
                        "Gzip + single-column barcode-like content (non-standard filename)"
                    ]
                # Could also be generic barcodes (not pure ACGT)
                if all(10 <= len(s) <= 30 for s in sample):
                    return FileDetectedType.MEX_BARCODES_10X, 0.70, [
                        "Gzip + single-column content, possible barcodes (non-standard filename)"
                    ]

        # Two-column file: could be features/genes (ID + symbol)
        if n_cols_first == 2:
            # Check for gene-like content in the second column
            gene_like = 0
            total = 0
            for line in lines[1:21]:
                if not line.strip():
                    continue
                parts = line.split("\t")
                if len(parts) >= 2:
                    total += 1
                    if _looks_like_gene_symbol(parts[1].strip().strip('"')):
                        gene_like += 1
            if total > 0 and gene_like / total >= 0.6:
                return FileDetectedType.MEX_FEATURES_10X, 0.78, [
                    "Gzip + two-column gene symbol content (non-standard filename)"
                ]
            # Check for Ensembl IDs in first column
            ensembl_like = sum(
                1 for line in lines[1:21] if line.strip()
                and bool(_re.match(r"^ENS[GT]\d{11}", line.split("\t")[0].strip().strip('"')))
            )
            if ensembl_like >= 3:
                return FileDetectedType.MEX_FEATURES_10X, 0.78, [
                    "Gzip + two-column Ensembl ID content (non-standard filename)"
                ]

        return None

    # ------------------------------------------------------------------
    # Content-first gzip helpers (Stage 1 & 2 of _classify_gzip)
    # ------------------------------------------------------------------

    @staticmethod
    def _is_gzip_wrapping_archive(file_path: Path) -> bool:
        """Check if a gzip file wraps a tar or zip archive.

        Reads through gzip decompression and checks the inner magic bytes
        for tar (ustar at offset 257) or zip (PK at offset 0).
        """
        import gzip
        try:
            with gzip.open(file_path, "rb") as fh:
                inner = fh.read(512)
        except Exception:
            return False
        if len(inner) < 4:
            return False
        if inner[:4] == b"PK\x03\x04":
            return True
        if len(inner) >= 262 and inner[257:262] == b"ustar":
            return True
        return False

    @staticmethod
    def _inspect_gzip_text_content(
        file_path: Path,
    ) -> tuple[FileDetectedType, float, list[str]] | None:
        """Deeper content inspection for gzip text files.

        Decompresses the header and classifies based on column count,
        header keywords, and sampled data values. Used when
        _inspect_gzip_mex_content() does not find MEX component patterns.
        """
        import csv
        import gzip

        try:
            with gzip.open(file_path, "rt", encoding="utf-8", errors="replace") as fh:
                header_text = fh.read(16384)
        except Exception:
            return None

        if not header_text:
            return None

        lines = header_text.split("\n")
        if not lines or not lines[0].strip():
            return None

        first_line = lines[0].strip()
        sep = "\t" if "\t" in first_line else ("," if "," in first_line else None)
        if sep is None:
            return None

        cols = first_line.split(sep)
        n_cols = len(cols)
        first_lower = first_line.lower()

        # Metadata header keywords
        meta_kw = [
            "barcode", "cell", "cluster", "cell_type", "celltype",
            "sample", "condition", "treatment", "tissue", "orig.ident",
            "model", "type", "group", "batch", "patient", "donor",
        ]
        # Gene annotation keywords
        gene_kw = ["gene", "symbol", "ensembl", "feature", "gene_id", "gene_name"]
        # Immune receptor keywords
        immune_kw = [
            "contig_id", "cdr3", "v_gene", "d_gene", "j_gene", "c_gene",
            "chain", "productive", "raw_clonotype_id",
        ]

        meta_score = sum(1 for kw in meta_kw if kw in first_lower)
        gene_score = sum(1 for kw in gene_kw if kw in first_lower)
        immune_score = sum(1 for kw in immune_kw if kw in first_lower)

        # Immune receptor detection by content
        if immune_score >= 3:
            return FileDetectedType.IMMUNE_RECEPTOR_DATA, 0.85, [
                "Gzip decompressed: V(D)J/immune receptor column headers detected"
            ]

        # Cell metadata by content
        if meta_score >= 2:
            return FileDetectedType.CELL_METADATA, 0.82, [
                f"Gzip decompressed: metadata headers ({meta_score} keywords)"
            ]

        # Gene annotation by content
        if gene_score >= 2:
            return FileDetectedType.GENE_ANNOTATION, 0.78, [
                f"Gzip decompressed: gene annotation headers ({gene_score} keywords)"
            ]

        # Matrix-like: many columns or numeric data
        is_matrix_like = n_cols > 10 or (n_cols >= 3 and meta_score == 0 and gene_score == 0)
        if is_matrix_like:
            # Sample data values for raw vs normalized
            data_stats = FileProbe._sample_data_values_gzip(file_path)
            is_raw, raw_evidence = FileProbe._classify_raw_vs_normalized(data_stats)

            if data_stats and data_stats.get("n_values", 0) > n_cols:
                if is_raw is True:
                    return FileDetectedType.RAW_UMI_MATRIX, 0.85, [
                        f"Gzip decompressed: {n_cols} cols, raw UMI counts "
                        f"(integer, sparse, {data_stats['sparsity']:.0%} zeros)"
                    ] + raw_evidence
                elif is_raw is False:
                    return FileDetectedType.NORMALIZED_MATRIX, 0.83, [
                        f"Gzip decompressed: {n_cols} cols, normalized data "
                        f"(float, {data_stats['sparsity']:.0%} zeros)"
                    ] + raw_evidence
                else:
                    dtype = FileDetectedType.CSV_MATRIX if sep == "," else FileDetectedType.TSV_MATRIX
                    return dtype, 0.70, [
                        f"Gzip decompressed: {n_cols} cols, numeric data (raw status uncertain)"
                    ] + raw_evidence

        # Few columns with no clear metadata/gene headers → could still be matrix
        if 2 <= n_cols <= 10 and meta_score == 0:
            data_stats = FileProbe._sample_data_values_gzip(file_path)
            if data_stats and data_stats.get("n_values", 0) > n_cols:
                is_raw, raw_evidence = FileProbe._classify_raw_vs_normalized(data_stats)
                dtype = FileDetectedType.RAW_UMI_MATRIX if is_raw is True else (
                    FileDetectedType.NORMALIZED_MATRIX if is_raw is False else
                    (FileDetectedType.CSV_MATRIX if sep == "," else FileDetectedType.TSV_MATRIX)
                )
                return dtype, 0.65, [
                    f"Gzip decompressed: {n_cols} cols with numeric data (small matrix)"
                ] + raw_evidence

        return None

    @staticmethod
    def _verify_gff_content(file_path: Path) -> bool:
        """Verify a gzip file contains GFF/GTF annotation by content structure."""
        import gzip
        try:
            with gzip.open(file_path, "rt", encoding="utf-8", errors="replace") as fh:
                # Read first 20 non-comment lines
                lines_read = 0
                valid_lines = 0
                for line in fh:
                    line = line.strip()
                    if not line or line.startswith("#"):
                        continue
                    parts = line.split("\t")
                    # GFF/GTF has 9 tab-separated fields
                    if len(parts) == 9:
                        # Check if col 4 and 5 are numeric (start, end)
                        try:
                            int(parts[3])
                            int(parts[4])
                            valid_lines += 1
                        except ValueError:
                            pass
                    lines_read += 1
                    if lines_read >= 20:
                        break
                return valid_lines >= 1
        except Exception:
            return False

    @staticmethod
    def _verify_peak_content(file_path: Path) -> bool:
        """Verify a gzip file contains BED/narrowPeak peak data by column structure."""
        import gzip
        try:
            with gzip.open(file_path, "rt", encoding="utf-8", errors="replace") as fh:
                lines_read = 0
                valid_lines = 0
                for line in fh:
                    line = line.strip()
                    if not line or line.startswith("#") or line.startswith("track"):
                        continue
                    parts = line.split("\t")
                    if len(parts) >= 3:
                        # BED: chr, start, end (col 2 and 3 must be numeric)
                        try:
                            int(parts[1])
                            int(parts[2])
                            valid_lines += 1
                        except ValueError:
                            pass
                    lines_read += 1
                    if lines_read >= 20:
                        break
                return valid_lines >= 1
        except Exception:
            return False

    @staticmethod
    def _verify_spatial_content(file_path: Path) -> bool:
        """Verify a gzip file contains spatial transcriptomics data."""
        import gzip
        try:
            with gzip.open(file_path, "rt", encoding="utf-8", errors="replace") as fh:
                header = fh.readline().strip()
                # Spatial position files have specific column patterns
                spatial_headers = [
                    "barcode", "tissue", "row", "col", "imagerow", "imagecol",
                    "in_tissue", "array_row", "array_col", "pxl_row_in_fullres",
                    "pxl_col_in_fullres",
                ]
                h_lower = header.lower()
                matches = sum(1 for kw in spatial_headers if kw in h_lower)
                return matches >= 3
        except Exception:
            return False

    @staticmethod
    def _sample_data_values_gzip(file_path: Path, n_samples: int = 100) -> dict | None:
        """Sample data values from a gzip-compressed delimited text file.

        Same logic as _sample_data_values() but for gzip files.
        """
        import csv
        import gzip

        try:
            with gzip.open(file_path, "rt", encoding="utf-8", errors="replace") as fh:
                first_line = fh.readline()
                if not first_line:
                    return None
                sep = "\t" if "\t" in first_line else ","
                reader = csv.reader(fh, delimiter=sep)
                values: list[float] = []
                for row_idx, row in enumerate(reader):
                    if row_idx >= n_samples:
                        break
                    if not row:
                        continue
                    for cell in row[1:]:
                        cell = cell.strip().strip('"')
                        if not cell:
                            continue
                        try:
                            values.append(float(cell))
                        except ValueError:
                            continue
                if not values:
                    return None
                integers = [v for v in values if v == int(v)]
                zeros = [v for v in values if v == 0.0]
                negatives = [v for v in values if v < 0]
                return {
                    "is_integer": len(integers) == len(values),
                    "has_negatives": len(negatives) > 0,
                    "sparsity": len(zeros) / len(values) if values else 0.0,
                    "value_range": (min(values), max(values)),
                    "n_values": len(values),
                }
        except Exception:
            return None

    def _classify_gzip_noncount_by_filename(
        self, name: str, inner_ext: str
    ) -> tuple[FileDetectedType, float, list[str]] | None:
        """Filename-based non-count filtering with content verification.

        This is Stage 3b of _classify_gzip — it only fires when Stage 1
        (content-first) and Stage 2 (specialized structure) did not produce
        a definitive answer. Confidence is reduced accordingly.
        """
        # Immune receptor data (BCR/TCR contig annotations) — NOT count matrices
        if any(kw in name for kw in [
            "filtered_contig_annotations_bcr",
            "filtered_contig_annotations_tcr",
            "contig_annotations",
            "vdj",
            "clonotypes",
        ]):
            return FileDetectedType.IMMUNE_RECEPTOR_DATA, 0.60, [
                "Gzip + BCR/TCR immune receptor filename (content ambiguous, confidence reduced)"
            ]

        # Peak data — NOT count matrices
        if inner_ext in (".bed", ".narrowpeak", ".broadpeak"):
            return FileDetectedType.PEAK_DATA, 0.65, [
                "Gzip + peak file extension (content ambiguous, confidence reduced)"
            ]
        if any(kw in name for kw in [
            "_peaks_", "_peak.", "chip_peaks", "atac_peaks",
            "narrowpeak", "broadpeak", "macs2_peaks",
        ]):
            return FileDetectedType.PEAK_DATA, 0.55, [
                "Gzip + peak filename pattern (content ambiguous, confidence reduced)"
            ]

        # Differential expression results — NOT count matrices
        if any(kw in name for kw in [
            "deseq2", "diffbind", "diffbind_result",
            "voom_weights", "limma_result", "edger_result",
            "differential_expression", "deg_",
        ]):
            return FileDetectedType.DIFF_RESULTS, 0.55, [
                "Gzip + differential expression filename (content ambiguous, confidence reduced)"
            ]

        return None

    # ------------------------------------------------------------------
    # Archive content peeking (Phase B1)
    # ------------------------------------------------------------------

    @staticmethod
    def _peek_archive_contents(archive_path: Path) -> list[str] | None:
        """List files inside an archive without full extraction.

        Supports ZIP, uncompressed TAR, and gzip-compressed TAR.
        Returns a list of filenames, or None if the archive cannot be read.

        This enables content-based archive classification without the
        overhead of extracting potentially large archives.
        """
        try:
            with open(archive_path, "rb") as fh:
                magic = fh.read(512)
        except (OSError, PermissionError):
            return None

        # ZIP archive
        if magic[:4] == b"PK\x03\x04":
            try:
                import zipfile
                with zipfile.ZipFile(archive_path, "r") as zf:
                    return zf.namelist()
            except Exception:
                return None

        # Gzip-compressed (possibly tar.gz)
        if magic[:2] == b"\x1f\x8b":
            try:
                import gzip
                import tarfile
                with gzip.open(archive_path, "rb") as gz_fh:
                    inner = gz_fh.read(512)
                if len(inner) >= 262 and inner[257:262] == b"ustar":
                    # gzip wrapping a tar — list tar contents
                    with tarfile.open(archive_path, "r:gz") as tf:
                        return [m.name for m in tf.getmembers() if m.isfile()]
                # gzip wrapping something else — just note it's compressed
                return ["<gzip-compressed-single-file>"]
            except Exception:
                return None

        # Uncompressed TAR (ustar magic at bytes 257-261)
        if len(magic) >= 262 and magic[257:262] == b"ustar":
            try:
                import tarfile
                # Try r:* first (auto-detect compression, seeks to end for GNU tar)
                with tarfile.open(archive_path, "r:*") as tf:
                    return [m.name for m in tf.getmembers() if m.isfile()]
            except Exception:
                # r:* can fail on truncated downloads — try streaming mode r:| which
                # doesn't seek to end; just reads entries sequentially
                try:
                    with tarfile.open(archive_path, "r:|") as tf:
                        return [m.name for m in tf.getmembers() if m.isfile()]
                except Exception:
                    return None

        # Also check for "ustar" at position 0 for tar files read from certain tools
        if magic[:5] == b"ustar":
            try:
                import tarfile
                with tarfile.open(archive_path, "r:*") as tf:
                    return [m.name for m in tf.getmembers() if m.isfile()]
            except Exception:
                return None

        return None

    @staticmethod
    def _classify_archive_by_contents(
        filenames: list[str], archive_path: Path | None = None
    ) -> tuple[FileDetectedType, float, list[str]]:
        """Classify an archive based on its internal file listing.

        This uses content-based archive peeking (_peek_archive_contents)
        instead of filename-based guessing. The archive type is determined by
        what files are inside, not by the archive's own name.
        """
        if not filenames:
            return FileDetectedType.ARCHIVE, 0.50, ["Empty archive — cannot classify"]

        names_lower = [f.lower() for f in filenames]
        all_names = " ".join(names_lower)
        evidence: list[str] = []

        # 10x MEX triplet detection
        has_mtx = any("matrix.mtx" in n for n in names_lower)
        has_barcodes = any("barcodes.tsv" in n for n in names_lower)
        has_features = any("features.tsv" in n for n in names_lower) or any(
            "genes.tsv" in n for n in names_lower
        )
        if has_mtx and has_barcodes and has_features:
            return FileDetectedType.ARCHIVE, 0.92, [
                "Archive contains 10x MEX triplet (matrix.mtx + barcodes.tsv + features/genes.tsv)"
            ]
        if has_mtx:
            evidence.append("Contains .mtx file(s)")

        # Object file detection
        has_h5ad = any(n.endswith(".h5ad") for n in names_lower)
        has_loom = any(n.endswith((".loom", ".loom.gz")) for n in names_lower)
        has_rds = any(n.endswith((".rds", ".rdata", ".rda")) for n in names_lower)
        has_h5seurat = any(n.endswith(".h5seurat") for n in names_lower)
        has_h5 = any(n.endswith((".h5", ".hdf5")) for n in names_lower)

        if has_h5ad:
            evidence.append("Contains .h5ad file(s)")
            return FileDetectedType.ARCHIVE, 0.88, evidence
        if has_loom:
            evidence.append("Contains .loom file(s)")
            return FileDetectedType.ARCHIVE, 0.85, evidence
        if has_rds:
            evidence.append("Contains .rds file(s)")
            return FileDetectedType.ARCHIVE, 0.82, evidence
        if has_h5seurat:
            evidence.append("Contains .h5seurat file(s)")
            return FileDetectedType.ARCHIVE, 0.82, evidence
        if has_h5:
            evidence.append("Contains .h5 file(s) — may be 10x HDF5")
            return FileDetectedType.ARCHIVE, 0.78, evidence

        # FASTQ detection
        has_fastq = any(
            n.endswith((".fastq", ".fastq.gz", ".fq", ".fq.gz"))
            for n in names_lower
        )
        if has_fastq:
            evidence.append("Contains FASTQ file(s)")
            return FileDetectedType.ARCHIVE, 0.80, evidence

        # CSV/TSV/TXT matrix detection
        has_delimited = any(
            n.endswith((".csv", ".csv.gz", ".tsv", ".tsv.gz", ".txt", ".txt.gz"))
            for n in names_lower
        )
        if has_delimited:
            evidence.append("Contains delimited text files — may be count matrices")
            return FileDetectedType.ARCHIVE, 0.70, evidence

        # Generic archive with recognized data files
        if evidence:
            return FileDetectedType.ARCHIVE, 0.65, evidence

        return FileDetectedType.ARCHIVE, 0.55, [
            f"Archive with {len(filenames)} file(s) — contents not recognized as scRNA-seq data"
        ]

    def _sample_data_values(
        self, file_path: Path, n_samples: int = 100
    ) -> dict | None:
        """Sample data values from a delimited text file to determine
        whether it contains raw counts or normalized data.

        Reads the first N data rows (skipping the header) and computes
        statistics on the numeric columns.

        Returns dict with:
        - is_integer: bool — all sampled values are integers
        - has_negatives: bool — any negative values present
        - sparsity: float — fraction of zero values
        - value_range: tuple — (min, max) of sampled values
        - n_values: int — number of values sampled

        Returns None if the file cannot be read or has no numeric data.
        """
        import csv
        import gzip

        try:
            if file_path.suffix.lower() == ".gz":
                opener = gzip.open
            else:
                opener = open

            with opener(file_path, "rt", encoding="utf-8", errors="replace") as fh:
                # Auto-detect delimiter from first line
                first_line = fh.readline()
                if not first_line:
                    return None
                sep = "\t" if "\t" in first_line else ","
                reader = csv.reader(fh, delimiter=sep)

                values: list[float] = []
                for row_idx, row in enumerate(reader):
                    if row_idx >= n_samples:
                        break
                    if not row:
                        continue
                    # Skip the first column (gene names / barcodes) and collect numeric values
                    for cell in row[1:]:
                        cell = cell.strip().strip('"')
                        if not cell:
                            continue
                        try:
                            val = float(cell)
                            values.append(val)
                        except ValueError:
                            continue

                if not values:
                    return None

                integers = [v for v in values if v == int(v)]
                zeros = [v for v in values if v == 0.0]
                negatives = [v for v in values if v < 0]

                return {
                    "is_integer": len(integers) == len(values),
                    "has_negatives": len(negatives) > 0,
                    "sparsity": len(zeros) / len(values) if values else 0.0,
                    "value_range": (min(values), max(values)),
                    "n_values": len(values),
                }
        except Exception:
            return None

    @staticmethod
    def _classify_raw_vs_normalized(data_stats: dict | None) -> tuple[bool | None, list[str]]:
        """Classify whether sampled data values indicate raw counts or normalized.

        Returns (is_raw_counts, evidence_list).
        is_raw_counts: True = confirmed raw, False = confirmed normalized, None = uncertain.
        """
        if data_stats is None:
            return None, ["Could not sample data values"]

        evidence: list[str] = []
        is_integer = data_stats["is_integer"]
        has_negatives = data_stats["has_negatives"]
        sparsity = data_stats["sparsity"]
        value_range = data_stats["value_range"]
        n_values = data_stats["n_values"]

        # Classic UMI count signature: integers, highly sparse, no negatives
        if is_integer and sparsity > 0.5 and not has_negatives:
            return True, [
                f"Raw UMI counts: integer values, {sparsity:.1%} zeros, "
                f"range [{value_range[0]:.0f}, {value_range[1]:.0f}], "
                f"no negative values ({n_values} values sampled)"
            ]

        # Log-normalized: negative values present (log of <1)
        if has_negatives:
            return False, [
                f"Normalized data: negative values detected "
                f"(log-transformed, {n_values} values sampled)"
            ]

        # Float values with bounded range → likely TPM/FPKM or other normalization
        if not is_integer and value_range[1] < 100:
            return False, [
                f"Normalized data: float values, range [{value_range[0]:.1f}, "
                f"{value_range[1]:.1f}] suggests TPM/FPKM or similar ({n_values} values sampled)"
            ]

        # Integer but low sparsity → could be normalized or bulk RNA-seq
        if is_integer and sparsity < 0.3:
            return None, [
                f"Ambiguous: integer values but low sparsity ({sparsity:.1%} zeros). "
                f"Could be normalized counts or bulk RNA-seq ({n_values} values sampled)"
            ]

        # Float but wide range → could be raw counts in float format or normalized
        if not is_integer and value_range[1] >= 100:
            return None, [
                f"Ambiguous: float values with wide range [{value_range[0]:.1f}, "
                f"{value_range[1]:.1f}]. Could be raw counts stored as float "
                f"or normalized data ({n_values} values sampled)"
            ]

        return None, [f"Uncertain: {n_values} values sampled, sparsity={sparsity:.1%}"]

    def _looks_like_text(self, header: bytes) -> bool:
        """Heuristic: does the header look like text (not binary)?"""
        if b"\x00" in header[:512]:
            return False
        text_bytes = sum(1 for b in header[:512] if 0x20 <= b <= 0x7E or b in (0x09, 0x0A, 0x0D))
        return text_bytes / min(len(header[:512]), 1) > 0.85

    def _decompress_gzip_and_reprobe(
        self, gz_path: Path, original_sidecar: FileSidecar, root: Path
    ) -> FileSidecar | None:
        """Decompress a single .gz file and re-classify based on its contents.

        Many GEO datasets distribute count matrices as .csv.gz or .txt.gz.
        Without decompression the probe can only guess from the filename;
        with decompression we can inspect the actual header rows to determine
        whether it is a count matrix, metadata table, or gene annotation.
        """
        import gzip

        try:
            # Read up to 16 KiB of decompressed content for header inspection
            with gzip.open(gz_path, "rb") as fh:
                header = fh.read(16384)
        except Exception as exc:
            self._logger.debug(
                "Failed to decompress gzip for re-probe",
                path=str(gz_path),
                error=str(exc),
            )
            return None

        if not header:
            return None

        # Try to classify from the decompressed content
        try:
            text = header.decode("utf-8", errors="replace")
        except UnicodeDecodeError:
            return None  # Binary content — keep original classification

        lines = text.split("\n")
        if not lines or not lines[0].strip():
            return None

        first_line = lines[0].strip()
        sep = "\t" if "\t" in first_line else ("," if "," in first_line else None)

        if sep is None:
            return None

        cols = first_line.split(sep)
        n_cols = len(cols)

        # Heuristic classification from column headers
        first_lower = first_line.lower()

        # Cell metadata: has cell/barcode/cluster/condition headers
        metadata_keywords = [
            "barcode", "cell", "cluster", "cell_type", "celltype",
            "sample", "condition", "treatment", "tissue", "orig.ident",
            "model", "type", "group", "batch", "patient", "donor",
        ]
        metadata_score = sum(1 for kw in metadata_keywords if kw in first_lower)

        # Gene annotation: has gene/symbol/ensembl/feature headers
        gene_keywords = ["gene", "symbol", "ensembl", "feature", "gene_id", "gene_name"]
        gene_score = sum(1 for kw in gene_keywords if kw in first_lower)

        # Count matrix: many columns, OR few columns with numeric data values.
        # The original threshold was n_cols > 10, but many real GEO datasets have
        # ≤10 samples. When the data values are numeric (especially integer/sparse),
        # it is far more likely a count matrix than metadata.
        is_matrix_like = n_cols > 10 or (
            n_cols >= 3
            and metadata_score == 0
            and gene_score == 0
        )
        if is_matrix_like:
            # Sample data values to determine raw vs normalized
            data_stats = self._sample_data_values(gz_path)
            is_raw, raw_evidence = self._classify_raw_vs_normalized(data_stats)

            # If data values are numeric (not metadata strings), treat as matrix
            if data_stats and data_stats.get("n_values", 0) > n_cols:
                if first_col := cols[0].strip().strip('"'):
                    # Check if first column looks like gene symbols
                    if _looks_like_gene_symbol(first_col):
                        # Gene symbols with many sample columns = count matrix
                        if is_raw is True:
                            dtype = FileDetectedType.RAW_UMI_MATRIX
                            confidence = 0.85
                        elif is_raw is False:
                            dtype = FileDetectedType.NORMALIZED_MATRIX
                            confidence = 0.83
                        else:
                            dtype = FileDetectedType.CSV_MATRIX if sep == "," else FileDetectedType.TSV_MATRIX
                            confidence = 0.78

                        return self._build_enriched_sidecar(
                            original_sidecar, dtype, confidence,
                            [f"Gzip decompressed: {n_cols} columns, gene-like first column, "
                             f"{len(lines)} rows visible"] + raw_evidence,
                            is_raw_counts=is_raw,
                        )
                    else:
                        # Could be a transposed matrix (cells as rows)
                        if is_raw is True:
                            dtype = FileDetectedType.RAW_UMI_MATRIX
                            confidence = 0.78
                        elif is_raw is False:
                            dtype = FileDetectedType.NORMALIZED_MATRIX
                            confidence = 0.76
                        else:
                            dtype = FileDetectedType.CSV_MATRIX if sep == "," else FileDetectedType.TSV_MATRIX
                            confidence = 0.70

                        return self._build_enriched_sidecar(
                            original_sidecar, dtype, confidence,
                            [f"Gzip decompressed: {n_cols} columns, {len(lines)} rows visible"] + raw_evidence,
                            is_raw_counts=is_raw,
                        )
            elif n_cols > 10:
                # No numeric data found but many columns — still likely a matrix
                # (could be all-zero first rows or unusual formatting)
                if is_raw is True:
                    dtype = FileDetectedType.RAW_UMI_MATRIX
                elif is_raw is False:
                    dtype = FileDetectedType.NORMALIZED_MATRIX
                else:
                    dtype = FileDetectedType.CSV_MATRIX if sep == "," else FileDetectedType.TSV_MATRIX
                return self._build_enriched_sidecar(
                    original_sidecar, dtype, 0.65,
                    [f"Gzip decompressed: {n_cols} columns, {len(lines)} rows visible"] + raw_evidence,
                    is_raw_counts=is_raw,
                )

        # Cell metadata
        if metadata_score >= 2:
            return self._build_enriched_sidecar(
                original_sidecar, FileDetectedType.CELL_METADATA, 0.82,
                [f"Gzip decompressed: metadata headers ({metadata_score} keywords)"],
            )

        # Gene annotation
        if gene_score >= 2:
            return self._build_enriched_sidecar(
                original_sidecar, FileDetectedType.GENE_ANNOTATION, 0.78,
                [f"Gzip decompressed: gene annotation headers ({gene_score} keywords)"],
            )

        # Fallback: delimited text with few columns — check data values first
        if n_cols <= 10 and n_cols >= 2:
            # Before concluding it's metadata, check if data values are numeric
            data_stats = self._sample_data_values(gz_path)
            if data_stats and data_stats.get("n_values", 0) > n_cols:
                # Has numeric data — more likely a small matrix than metadata
                is_raw, raw_evidence = self._classify_raw_vs_normalized(data_stats)
                dtype = FileDetectedType.RAW_UMI_MATRIX if is_raw is True else (
                    FileDetectedType.NORMALIZED_MATRIX if is_raw is False else
                    (FileDetectedType.CSV_MATRIX if sep == "," else FileDetectedType.TSV_MATRIX)
                )
                return self._build_enriched_sidecar(
                    original_sidecar, dtype, 0.70,
                    [f"Gzip decompressed: {n_cols} columns with numeric data, "
                     f"classifying as matrix ({data_stats['n_values']} values sampled)"] + raw_evidence,
                    is_raw_counts=is_raw,
                )
            return self._build_enriched_sidecar(
                original_sidecar, FileDetectedType.CELL_METADATA, 0.65,
                [f"Gzip decompressed: {n_cols} columns, looks like metadata"],
            )

        return None

    def _build_enriched_sidecar(
        self,
        original: FileSidecar,
        new_type: FileDetectedType,
        confidence: float,
        extra_evidence: list[str],
        is_raw_counts: bool | None = None,
    ) -> FileSidecar:
        """Build a new sidecar with updated classification, keeping original data."""
        evidence = list(original.evidence) + extra_evidence
        sc = FileSidecar(
            file_id=original.file_id,
            path=original.path,
            filename=original.filename,
            size_bytes=original.size_bytes,
            extension=original.extension,
            mime_type=original.mime_type,
            detected_type=new_type,
            processing_stage=self._infer_stage(new_type),
            modality=original.modality,
            organism=original.organism,
            sample_id_guess=original.sample_id_guess,
            gse=original.gse,
            gsm=original.gsm,
            srr=original.srr,
            paired_files=list(original.paired_files),
            confidence=confidence,
            evidence=evidence,
            is_raw_counts=is_raw_counts,
            limited_use=original.limited_use,
            parent_file_ids=list(original.parent_file_ids),
            extracted_files=list(original.extracted_files),
            probe_version=original.probe_version,
            checksum_sha256=original.checksum_sha256,
        )
        sc.recommended_routes = self._recommend_routes(sc)
        return sc

    # File extensions that look like archives but are actually 10x Genomics
    # compressed data files. These MUST NOT be extracted recursively.
    _FORMAT_EXEMPT_EXTENSIONS: set[str] = {
        ".mtx.gz", ".tsv.gz", ".csv.gz", ".txt.gz", ".h5ad.gz",
    }

    # Maximum nesting depth for recursive archive extraction to prevent
    # infinite loops from circular or pathological archive structures.
    _MAX_NESTING_DEPTH: int = 5

    def _extract_and_probe(
        self, archive_path: Path, parent_id: str, root: Path
    ) -> list[FileSidecar]:
        """Extract an archive and recursively probe its contents.

        Handles nested archives (Russian Doll structure) by recursing into
        extracted directories. Uses format-aware exemption to avoid
        extracting 10x compressed data files (.mtx.gz, .tsv.gz, .csv.gz).

        Ported from v1 ``TarBlindBoxUnpacker._recursive_unpack()``.
        """
        extract_dir = root / f"_extracted_{archive_path.stem}"
        try:
            self._extract_archive(archive_path, extract_dir)
        except Exception as exc:
            self._logger.warning("Failed to extract archive", path=str(archive_path), error=str(exc))
            return []

        return self._probe_extracted_recursive(
            extract_dir, parent_id, root, depth=0
        )

    def _probe_extracted_recursive(
        self,
        directory: Path,
        parent_id: str,
        root: Path,
        depth: int = 0,
    ) -> list[FileSidecar]:
        """Recursively probe an extracted directory, handling nested archives.

        For each file found:
        - If it is a 10x compressed data file (format-exempt) → probe as data.
        - If it is a known archive format → extract recursively (respecting
          max nesting depth).
        - Otherwise → probe as normal file.

        Args:
            directory: The extracted directory to scan.
            parent_id: File ID of the parent archive.
            root: Root directory for relative path computation.
            depth: Current nesting depth (0-based).

        Returns:
            List of FileSidecar for all discovered files.
        """
        if depth >= self._MAX_NESTING_DEPTH:
            self._logger.warning(
                "Max archive nesting depth reached — skipping nested extraction",
                directory=str(directory), depth=depth,
            )
            return []

        sidecars: list[FileSidecar] = []

        for file_path in directory.rglob("*"):
            if not file_path.is_file():
                continue

            filename_lower = file_path.name.lower()

            # Format-aware exemption: 10x compressed data files are NOT archives
            if any(filename_lower.endswith(ext) for ext in self._FORMAT_EXEMPT_EXTENSIONS):
                sc = self.probe_file(file_path, relative_root=root)
                if sc:
                    sc.parent_file_ids.append(parent_id)
                    sidecars.append(sc)
                continue

            # Detect and recursively extract nested archives
            if self._is_archive_file(file_path):
                self._logger.info("Cracking nested archive", path=str(file_path), depth=depth + 1)
                nested_sidecars = self._extract_and_probe(
                    file_path, parent_id, root
                )
                sidecars.extend(nested_sidecars)
                continue

            # Regular file — probe normally
            sc = self.probe_file(file_path, relative_root=root)
            if sc:
                sc.parent_file_ids.append(parent_id)
                sidecars.append(sc)

        return sidecars

    @staticmethod
    def _is_archive_file(file_path: Path) -> bool:
        """Check if a file is an extractable archive (multi-file container).

        Only matches archives that can contain multiple files: tar, tar.gz, tgz, zip.
        Plain ``.gz`` files (gzip-compressed single files like ``.csv.gz``,
        ``.txt.gz``) are NOT archives — they are decompressed and probed as
        data by ``_decompress_gzip_and_reprobe``.
        """
        name = file_path.name.lower()
        # Extension-based check
        if name.endswith((".tar.gz", ".tgz", ".tar", ".zip")):
            return True
        # Magic byte check for misnamed zip files only.
        # We deliberately do NOT check gzip magic here: a bare .gz file
        # is a single compressed file, not a multi-file archive.  Treating
        # it as an archive causes _extract_archive to fail with "Unsupported
        # archive format", which kills the entire extraction branch for that
        # file — losing data that should have been decompressed and probed.
        try:
            with open(file_path, "rb") as fh:
                magic = fh.read(4)
            if magic == b"PK\x03\x04":      # ZIP magic (4 bytes)
                return True
        except (OSError, PermissionError):
            pass
        return False

    def _extract_archive(self, archive_path: Path, dest: Path) -> None:
        """Extract zip or tar archive."""
        dest.mkdir(parents=True, exist_ok=True)
        if zipfile.is_zipfile(archive_path):
            with zipfile.ZipFile(archive_path, "r") as zf:
                zf.extractall(dest)
        elif tarfile.is_tarfile(archive_path):
            with tarfile.open(archive_path, "r:*") as tf:
                tf.extractall(dest)
        else:
            raise ValueError(f"Unsupported archive format: {archive_path}")


def _looks_like_gene_symbol(text: str) -> bool:
    """Check if text looks like a gene symbol (not a barcode or sample ID).

    Gene symbols in real scRNA-seq data follow species-specific casing:
    - Human: all UPPERCASE with digits (BRCA1, TP53, GAPDH, ACTB)
    - Mouse: Titlecase (Foxp3, Gapdh, Actb)
    - Common: short alphanumeric, may have hyphens/dots (H2-K1, RP11-34P13.7)

    All-lowercase words like "value", "header", "sample" are NOT gene symbols.
    """
    import re
    if not text:
        return False
    # Barcode check: long pure-ACGT string (case-insensitive — barcodes are ACGT only)
    if re.match(r"^[ACGT]{14,}$", text, re.IGNORECASE):
        return False
    # Numeric-only check (sample index)
    if re.match(r"^\d+$", text):
        return False
    # Human-style UPPERCASE genes: TP53, BRCA1, GAPDH, ACTB, H2-K1
    if re.match(r"^[A-Z][A-Z0-9\-\.]{1,14}$", text):
        return True
    # Mouse-style Titlecase genes: Foxp3, Gapdh, Actb
    if re.match(r"^[A-Z][a-z][A-Za-z0-9\-\.]{0,13}$", text):
        return True
    # Ensembl IDs: ENSG00000141510, ENSMUSG00000051951
    if re.match(r"^ENS[A-Z]{0,4}G\d{11}$", text, re.IGNORECASE):
        return True
    return False
