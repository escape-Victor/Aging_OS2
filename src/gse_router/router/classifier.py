"""File type classifier — groups related files into logical input groups.

Takes a collection of file sidecars and organizes them into biologically
meaningful groups: 10x MEX triplets, paired FASTQ runs, standalone objects, etc.
"""

from __future__ import annotations

import gzip
import re
from pathlib import Path

import structlog

from ..schemas.manifest import FileGroup, FileGroupType, Strategy
from ..schemas.sidecar import FileDetectedType, FileSidecar

logger = structlog.get_logger(__name__)


class FileClassifier:
    """Groups individual file sidecars into logical input groups.

    A "group" represents one logical input to a processing route:
    a 10x MEX triplet (matrix + barcodes + features), a single h5ad file,
    a pair of FASTQ files, etc.
    """

    def __init__(self) -> None:
        self._logger = logger.bind(component="FileClassifier")

    def group_files(
        self,
        sidecars: dict[str, FileSidecar],
        data_dir: str | None = None,
        multi_triplet_strategy: str = "first",
    ) -> list[FileGroup]:
        """Group sidecar'd files into logical processing units.

        Args:
            sidecars: Dict mapping file_id -> FileSidecar.
            data_dir: Optional root directory for the dataset. When provided,
                enables shape-matching fallback for MEX triplets with
                non-standard filenames.
            multi_triplet_strategy: How to resolve multiple MEX triplets in the
                same directory. Options:
                - ``"first"`` — first candidate (backward compatible, may
                  mismatch across samples)
                - ``"prefix_match"`` — match by longest common filename prefix
                  (same sample shares prefix)
                - ``"shape_only"`` — only group via shape-matching, skip
                  name-based grouping

        Returns:
            List of FileGroup, each representing one logical input unit.
        """
        groups: list[FileGroup] = []
        assigned: set[str] = set()

        # 1. Group 10x MEX triplets
        mex_groups = self._group_mex_triplets(
            sidecars, assigned, data_dir=data_dir,
            multi_triplet_strategy=multi_triplet_strategy,
        )
        groups.extend(mex_groups)

        # 2. Group 10x H5 files (each is standalone)
        h5_groups = self._group_h5_files(sidecars, assigned)
        groups.extend(h5_groups)

        # 3. Standalone processed objects
        obj_groups = self._group_standalone_objects(sidecars, assigned)
        groups.extend(obj_groups)

        # 4. FASTQ files (paired or single)
        fq_groups = self._group_fastq_files(sidecars, assigned)
        groups.extend(fq_groups)

        # 5. Remaining files
        remaining_groups = self._group_remaining(sidecars, assigned)
        groups.extend(remaining_groups)

        return groups

    def select_strategy(
        self,
        groups: list[FileGroup],
        prefer_processed: bool = True,
    ) -> Strategy:
        """Select the overall processing strategy based on available groups.

        Args:
            groups: All file groups.
            prefer_processed: Whether to prefer processed counts over raw FASTQ.

        Returns:
            The recommended Strategy.
        """
        has_processed = any(
            g.group_type
            in {
                FileGroupType.MEX_10X,
                FileGroupType.H5_10X,
                FileGroupType.H5AD,
                FileGroupType.LOOM,
                FileGroupType.SEURAT_RDS,
                FileGroupType.SCE_RDS,
                FileGroupType.H5SEURAT,
                FileGroupType.CSV_MATRIX,
                FileGroupType.RAW_UMI_MATRIX,
                FileGroupType.ARCHIVE_10X,  # Contains 10x MEX — extractable processed data
            }
            for g in groups
        )
        has_normalized_only = any(
            g.group_type == FileGroupType.NORMALIZED_MATRIX
            for g in groups
        ) and not has_processed
        has_raw = any(
            g.group_type in {FileGroupType.FASTQ_PAIRED, FileGroupType.FASTQ_SINGLE}
            for g in groups
        )
        has_unknown = any(g.group_type == FileGroupType.UNKNOWN for g in groups)
        has_metadata = any(
            g.group_type in (FileGroupType.METADATA, FileGroupType.ANNOTATION)
            for g in groups
        )
        no_usable = not has_processed and not has_raw and not has_normalized_only

        if no_usable:
            if has_metadata:
                return Strategy.METADATA_ONLY
            return Strategy.UNSOLVABLE
        if has_normalized_only:
            return Strategy.NORMALIZED_ONLY
        if has_processed and prefer_processed:
            if has_raw and has_unknown:
                return Strategy.HYBRID
            return Strategy.PROCESSED_COUNTS_FIRST
        if has_raw:
            return Strategy.RAW_FASTQ_ONLY
        return Strategy.METADATA_ONLY

    # ------------------------------------------------------------------
    # Grouping strategies
    # ------------------------------------------------------------------

    def _group_mex_triplets(
        self,
        sidecars: dict[str, FileSidecar],
        assigned: set[str],
        data_dir: str | None = None,
        multi_triplet_strategy: str = "first",
    ) -> list[FileGroup]:
        """Find and group 10x MEX matrix/barcodes/features triplets.

        Uses filename-based matching first. Falls back to mathematical
        shape-matching (line-count dimension verification) when standard
        filename patterns don't find both barcode and feature companions.

        Args:
            sidecars: Dict of file_id -> FileSidecar.
            assigned: Set of already-assigned file IDs.
            data_dir: Optional root directory for shape-matching fallback.
            multi_triplet_strategy: How to pick among multiple candidates:
                ``"first"``, ``"prefix_match"``, or ``"shape_only"``.
        """
        groups: list[FileGroup] = []
        matrices = [
            (fid, sc)
            for fid, sc in sidecars.items()
            if sc.detected_type == FileDetectedType.MEX_MATRIX_10X and fid not in assigned
        ]

        for mid, mtx in matrices:
            # Find barcodes and features in the same directory
            mtx_dir = self._dir_of(mtx.path)
            barcode_candidates = [
                (fid, sc)
                for fid, sc in sidecars.items()
                if sc.detected_type == FileDetectedType.MEX_BARCODES_10X
                and self._dir_of(sc.path) == mtx_dir
                and fid not in assigned
            ]
            feature_candidates = [
                (fid, sc)
                for fid, sc in sidecars.items()
                if sc.detected_type == FileDetectedType.MEX_FEATURES_10X
                and self._dir_of(sc.path) == mtx_dir
                and fid not in assigned
            ]

            if barcode_candidates and feature_candidates:
                if multi_triplet_strategy == "shape_only":
                    # Skip name-based grouping; only use shape-matching
                    pass
                elif multi_triplet_strategy == "prefix_match" and (len(barcode_candidates) > 1 or len(feature_candidates) > 1):
                    # Match by longest common filename prefix (same sample)
                    bc_fid, bc = self._pick_best_candidate(mtx, barcode_candidates)
                    ft_fid, ft = self._pick_best_candidate(mtx, feature_candidates)
                else:
                    # "first" strategy — backward compatible
                    bc_fid, bc = barcode_candidates[0]
                    ft_fid, ft = feature_candidates[0]

                group = FileGroup(
                    group_id=f"mex_{mid[:12]}",
                    group_type=FileGroupType.MEX_10X,
                    sample_id=mtx.sample_id_guess,
                    gsm=mtx.gsm,
                    files={
                        "matrix": mid,
                        "barcodes": bc_fid,
                        "features": ft_fid,
                    },
                    route="read_10x_mtx_to_h5ad",
                    route_description="Convert 10x MEX triplet to h5ad via Scanpy",
                    confidence=min(mtx.confidence, bc.confidence, ft.confidence),
                    notes=[
                        f"Matrix: {mtx.filename}",
                        f"Barcodes: {bc.filename}",
                        f"Features: {ft.filename}",
                    ],
                )
                groups.append(group)
                assigned.update({mid, bc_fid, ft_fid})
            elif data_dir is not None:
                # Fallback: try mathematical shape-matching for non-standard filenames
                shape_group = self._shape_match_mex(
                    mid, mtx, sidecars, assigned, data_dir
                )
                if shape_group:
                    groups.append(shape_group)

        return groups

    def _group_h5_files(
        self, sidecars: dict[str, FileSidecar], assigned: set[str]
    ) -> list[FileGroup]:
        """Group 10x H5 files (each is a standalone group)."""
        groups: list[FileGroup] = []
        for fid, sc in sidecars.items():
            if sc.detected_type == FileDetectedType.H5_10X and fid not in assigned:
                group = FileGroup(
                    group_id=f"h5_{fid[:12]}",
                    group_type=FileGroupType.H5_10X,
                    sample_id=sc.sample_id_guess,
                    gsm=sc.gsm,
                    files={"h5": fid},
                    route="read_10x_h5_to_h5ad",
                    route_description="Convert 10x H5 to h5ad via Scanpy",
                    confidence=sc.confidence,
                )
                groups.append(group)
                assigned.add(fid)
        return groups

    def _group_standalone_objects(
        self, sidecars: dict[str, FileSidecar], assigned: set[str]
    ) -> list[FileGroup]:
        """Group standalone processed objects (h5ad, loom, rds, raw UMI, normalized)."""
        type_map = {
            FileDetectedType.H5AD: (FileGroupType.H5AD, "read_h5ad_validate"),
            FileDetectedType.LOOM: (FileGroupType.LOOM, "read_loom_to_h5ad"),
            FileDetectedType.RDS_SEURAT: (FileGroupType.SEURAT_RDS, "read_rds_seurat_to_h5ad"),
            FileDetectedType.RDS_SCE: (FileGroupType.SCE_RDS, "read_rds_sce_to_h5ad"),
            FileDetectedType.H5SEURAT: (FileGroupType.H5SEURAT, "read_rds_seurat_to_h5ad"),
            FileDetectedType.RAW_UMI_MATRIX: (FileGroupType.RAW_UMI_MATRIX, "read_raw_umi_matrix_to_h5ad"),
            FileDetectedType.NORMALIZED_MATRIX: (FileGroupType.NORMALIZED_MATRIX, "read_normalized_matrix_to_h5ad"),
            FileDetectedType.CSV_MATRIX: (FileGroupType.CSV_MATRIX, "read_csv_matrix_to_h5ad"),
            FileDetectedType.TSV_MATRIX: (FileGroupType.CSV_MATRIX, "read_csv_matrix_to_h5ad"),
            FileDetectedType.TXT_MATRIX: (FileGroupType.CSV_MATRIX, "read_csv_matrix_to_h5ad"),
        }
        groups: list[FileGroup] = []
        for fid, sc in sidecars.items():
            if sc.detected_type in type_map and fid not in assigned:
                gtype, route = type_map[sc.detected_type]
                group = FileGroup(
                    group_id=f"{gtype.value}_{fid[:12]}",
                    group_type=gtype,
                    sample_id=sc.sample_id_guess,
                    gsm=sc.gsm,
                    files={"data": fid},
                    route=route,
                    route_description=f"Convert {sc.detected_type.value} to h5ad",
                    confidence=sc.confidence,
                )
                groups.append(group)
                assigned.add(fid)
        return groups

    def _group_fastq_files(
        self, sidecars: dict[str, FileSidecar], assigned: set[str]
    ) -> list[FileGroup]:
        """Group FASTQ files — detect paired-end vs single-end."""
        groups: list[FileGroup] = []
        fq_files = [
            (fid, sc)
            for fid, sc in sidecars.items()
            if sc.detected_type == FileDetectedType.FASTQ and fid not in assigned
        ]

        # Simple paired-end detection: look for R1/R2 pattern in same directory
        paired: set[str] = set()
        for i, (fid1, sc1) in enumerate(fq_files):
            if fid1 in paired:
                continue
            name1 = sc1.filename.lower()
            dir1 = self._dir_of(sc1.path)

            # Try to find R2 mate
            mate = None
            for fid2, sc2 in fq_files[i + 1 :]:
                if fid2 in paired:
                    continue
                name2 = sc2.filename.lower()
                if (
                    self._dir_of(sc2.path) == dir1
                    and self._are_paired(name1, name2)
                ):
                    mate = (fid2, sc2)
                    break

            if mate:
                fid2, sc2 = mate
                group = FileGroup(
                    group_id=f"fastq_paired_{fid1[:12]}",
                    group_type=FileGroupType.FASTQ_PAIRED,
                    sample_id=sc1.sample_id_guess,
                    gsm=sc1.gsm or sc2.gsm,
                    files={"read1": fid1, "read2": fid2},
                    route="cellranger_count",
                    route_description="Cell Ranger count from paired FASTQ",
                    confidence=0.85,
                )
                groups.append(group)
                assigned.update({fid1, fid2})
                paired.update({fid1, fid2})
            else:
                group = FileGroup(
                    group_id=f"fastq_single_{fid1[:12]}",
                    group_type=FileGroupType.FASTQ_SINGLE,
                    sample_id=sc1.sample_id_guess,
                    gsm=sc1.gsm,
                    files={"reads": fid1},
                    route="starsolo",
                    route_description="STARsolo from single-end FASTQ",
                    confidence=0.80,
                )
                groups.append(group)
                assigned.add(fid1)

        return groups

    def _group_remaining(
        self, sidecars: dict[str, FileSidecar], assigned: set[str]
    ) -> list[FileGroup]:
        """Group remaining unassigned files (metadata, annotations, unknowns)."""
        groups: list[FileGroup] = []
        for fid, sc in sidecars.items():
            if fid in assigned:
                continue
            gtype = self._map_detected_to_group(sc.detected_type, sidecar=sc)

            # Determine route based on group type and detected type
            if gtype == FileGroupType.ARCHIVE_10X:
                route = "extract_10x_archive_to_h5ad"
            elif sc.detected_type in (FileDetectedType.CELL_METADATA, FileDetectedType.SAMPLE_METADATA):
                route = "merge_cell_metadata"
            elif sc.detected_type == FileDetectedType.METADATA:
                # GEO manifest/checksum files — not meaningful processing targets, skip
                assigned.add(fid)
                continue
            elif sc.detected_type == FileDetectedType.GENE_ANNOTATION:
                route = "merge_gene_annotation"
            else:
                route = "manual_review"

            group = FileGroup(
                group_id=f"{gtype.value}_{fid[:12]}",
                group_type=gtype,
                sample_id=sc.sample_id_guess,
                gsm=sc.gsm,
                files={"file": fid},
                route=route,
                confidence=sc.confidence,
            )
            groups.append(group)
            assigned.add(fid)
        return groups

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _dir_of(path: str) -> str:
        return str(path).rsplit("/", 1)[0] if "/" in str(path) else "."

    @staticmethod
    def _pick_best_candidate(
        matrix: FileSidecar,
        candidates: list[tuple[str, FileSidecar]],
    ) -> tuple[str, FileSidecar]:
        """Pick the best companion file for a matrix from multiple candidates.

        When multiple barcode or feature files exist in the same directory,
        this matches by longest common filename prefix — files from the same
        sample typically share a prefix (e.g., ``GSM123_sampleA_barcodes.tsv``
        and ``GSM123_sampleA_matrix.mtx``).

        Falls back to the first candidate when prefix matching finds no clear
        winner (ties or zero-length common prefix).

        Args:
            matrix: The matrix sidecar (reference for prefix matching).
            candidates: List of (file_id, FileSidecar) tuples.

        Returns:
            The best-matching (file_id, FileSidecar) tuple.
        """
        if len(candidates) == 1:
            return candidates[0]

        import os

        def _common_prefix_len(a: str, b: str) -> int:
            """Length of longest common prefix between two strings."""
            n = 0
            for ca, cb in zip(a, b):
                if ca != cb:
                    break
                n += 1
            return n

        matrix_name = os.path.basename(str(matrix.path))
        best = candidates[0]
        best_score = _common_prefix_len(matrix_name, os.path.basename(str(candidates[0][1].path)))

        for candidate in candidates[1:]:
            cand_name = os.path.basename(str(candidate[1].path))
            score = _common_prefix_len(matrix_name, cand_name)
            if score > best_score:
                best_score = score
                best = candidate

        return best

    @staticmethod
    def _are_paired(name1: str, name2: str) -> bool:
        """Detect paired-end FASTQ filenames (R1/R2 pattern)."""
        # Replace R1/R2 with placeholder and compare
        n1 = name1.replace("_r1_", "_rx_").replace("_r2_", "_rx_").replace("_R1_", "_RX_").replace("_R2_", "_RX_")
        n1 = n1.replace(".r1.", ".rx.").replace(".r2.", ".rx.")
        n2 = name2.replace("_r1_", "_rx_").replace("_r2_", "_rx_").replace("_R1_", "_RX_").replace("_R2_", "_RX_")
        n2 = n2.replace(".r1.", ".rx.").replace(".r2.", ".rx.")
        return n1 == n2

    # Evidence phrases from content-based archive peeking that indicate
    # the archive contains 10x Genomics data (MEX triplet or H5 files).
    # These are produced by FileProbe._classify_archive_by_contents().
    TENX_ARCHIVE_EVIDENCE_KEYWORDS = [
        "10x MEX triplet",
        "10x MEX bundle",
        ".h5 file(s) — may be 10x",
        "contains .h5 file",
        "matrix.mtx",
    ]

    @staticmethod
    def _map_detected_to_group(
        detected_type: FileDetectedType,
        sidecar: FileSidecar | None = None,
    ) -> FileGroupType:
        """Map detected file type to group type.

        When a sidecar is provided, checks evidence for 10x archive
        hints (from content-peeking) to upgrade generic ARCHIVE → ARCHIVE_10X.
        No longer uses filename patterns — content-based evidence only.
        """
        # Upgrade ARCHIVE → ARCHIVE_10X when content-peeking evidence
        # indicates 10x data inside the archive
        if detected_type == FileDetectedType.ARCHIVE and sidecar is not None:
            evidence_text = " ".join(sidecar.evidence).lower()
            if any(kw.lower() in evidence_text for kw in FileClassifier.TENX_ARCHIVE_EVIDENCE_KEYWORDS):
                return FileGroupType.ARCHIVE_10X

        mapping: dict[FileDetectedType, FileGroupType] = {
            FileDetectedType.CELL_METADATA: FileGroupType.METADATA,
            FileDetectedType.SAMPLE_METADATA: FileGroupType.METADATA,
            FileDetectedType.GENE_ANNOTATION: FileGroupType.ANNOTATION,
            FileDetectedType.CSV_MATRIX: FileGroupType.CSV_MATRIX,
            FileDetectedType.TSV_MATRIX: FileGroupType.CSV_MATRIX,
            FileDetectedType.TXT_MATRIX: FileGroupType.CSV_MATRIX,
            FileDetectedType.RAW_UMI_MATRIX: FileGroupType.RAW_UMI_MATRIX,
            FileDetectedType.NORMALIZED_MATRIX: FileGroupType.NORMALIZED_MATRIX,
            FileDetectedType.MATRIX_MAT: FileGroupType.CSV_MATRIX,  # Treat as generic matrix — converter handles .mat specifics
            FileDetectedType.BAM: FileGroupType.BAM,
            FileDetectedType.CRAM: FileGroupType.BAM,
            FileDetectedType.ARCHIVE: FileGroupType.ARCHIVE,
            FileDetectedType.RAW_TAR: FileGroupType.ARCHIVE,
            FileDetectedType.FILELIST_TXT: FileGroupType.METADATA,
            FileDetectedType.METADATA: FileGroupType.METADATA,  # GEO manifest / checksum files
            # Non-count supplementary data — explicitly map to UNKNOWN so they
            # are never routed as count matrices (BCR/TCR, peaks, diff results)
            FileDetectedType.IMMUNE_RECEPTOR_DATA: FileGroupType.UNKNOWN,
            FileDetectedType.PEAK_DATA: FileGroupType.UNKNOWN,
            FileDetectedType.DIFF_RESULTS: FileGroupType.UNKNOWN,
        }
        return mapping.get(detected_type, FileGroupType.UNKNOWN)

    # ------------------------------------------------------------------
    # Shape-matching helpers (ported from v1 UniversalOmicsIntegrator)
    # ------------------------------------------------------------------

    @staticmethod
    def _read_mtx_dimensions(mtx_abs_path: str) -> tuple[int, int] | None:
        """Read matrix dimensions from a MatrixMarket header.

        Handles both plain .mtx and gzipped .mtx.gz files.
        Reads only the header portion (up to the size line) to avoid
        loading the full matrix into memory.

        Args:
            mtx_abs_path: Absolute path to the .mtx or .mtx.gz file.

        Returns:
            (rows, cols) tuple, or None if the file cannot be read.
        """
        try:
            if mtx_abs_path.endswith(".gz"):
                with gzip.open(mtx_abs_path, "rt", encoding="utf-8", errors="ignore") as fh:
                    return FileClassifier._parse_mtx_header(fh)
            else:
                with open(mtx_abs_path, "r", encoding="utf-8", errors="ignore") as fh:
                    return FileClassifier._parse_mtx_header(fh)
        except (OSError, ValueError, gzip.BadGzipFile) as exc:
            logger.warning("Failed to read MTX dimensions", path=mtx_abs_path, error=str(exc))
            return None

    @staticmethod
    def _parse_mtx_header(fh) -> tuple[int, int] | None:
        """Parse the MatrixMarket header from an open file handle.

        The MatrixMarket format is:
            %%MatrixMarket matrix coordinate <type> <format>
            % comments (optional, start with %)
            <rows> <cols> <nnz>

        Returns (rows, cols) or None.
        """
        for line in fh:
            line = line.strip()
            if not line or line.startswith("%"):
                continue
            parts = line.split()
            if len(parts) >= 2:
                return int(parts[0]), int(parts[1])
        return None

    @staticmethod
    def _count_lines_in_file(filepath: str) -> int:
        """Count lines in a text file, handling gzip transparently.

        Args:
            filepath: Absolute path to the file.

        Returns:
            Number of lines, or 0 on any error.
        """
        try:
            if filepath.endswith(".gz"):
                with gzip.open(filepath, "rt", encoding="utf-8", errors="ignore") as fh:
                    return sum(1 for _ in fh)
            else:
                with open(filepath, "r", encoding="utf-8", errors="ignore") as fh:
                    return sum(1 for _ in fh)
        except (OSError, gzip.BadGzipFile):
            return 0

    def _shape_match_mex(
        self,
        matrix_fid: str,
        mtx_sidecar: FileSidecar,
        sidecars: dict[str, FileSidecar],
        assigned: set[str],
        data_dir: str,
    ) -> FileGroup | None:
        """Fallback: use mathematical dimension matching to find MEX barcodes/features.

        When standard filename-based detection fails to identify barcode and
        feature companion files, this method reads the MatrixMarket dimensions
        and finds companion files by matching their line counts to the matrix
        shape. This handles non-standard filenames (e.g., ``*cellBarcodes.txt.gz``,
        ``SAM*.genes.tsv.gz``) that the probe layer cannot identify by name alone.

        Ported from v1 ``UniversalOmicsIntegrator.process_and_get_tracks()``
        lines 177–235.

        Algorithm:
        1. Read matrix dimensions (A rows × B cols) from the MTX header.
        2. For each unassigned file in the same directory, count its lines.
        3. Candidates_A = files with line_count == A
        4. Candidates_B = files with line_count == B
        5. Use filename heuristics (barc/feat/gene) to identify which is which:
           - Files with "barc" / "cell" / "meta" in name → likely barcodes
           - Files with "feat" / "gene" in name → likely features
        6. If A != B and heuristics allow identification → create group.
        7. If A != B and heuristics fail: larger dim → features, smaller → barcodes.
        8. If no valid pair found → return None.

        Args:
            matrix_fid: File ID of the MTX matrix.
            mtx_sidecar: Sidecar of the MTX matrix.
            sidecars: All file sidecars (for finding candidate files).
            assigned: Set of already-assigned file IDs.
            data_dir: Root directory of the dataset for constructing absolute paths.

        Returns:
            A FileGroup if matching succeeds, or None.
        """
        mtx_dir = self._dir_of(mtx_sidecar.path)
        mtx_abs = str(Path(data_dir) / mtx_sidecar.path)
        dims = self._read_mtx_dimensions(mtx_abs)
        if dims is None:
            self._logger.debug("Could not read MTX dimensions for shape-matching",
                               file=mtx_sidecar.filename)
            return None

        A, B = dims
        if A == 0 or B == 0:
            return None

        # Collect line counts for all unassigned TSV/CSV/TXT files in same directory
        dir_candidates: list[tuple[str, FileSidecar, int]] = []
        for fid, sc in sidecars.items():
            if fid == matrix_fid or fid in assigned:
                continue
            if self._dir_of(sc.path) != mtx_dir:
                continue
            # Only consider text-like files
            fname = sc.filename.lower()
            ext = Path(fname).suffix
            is_text_ext = any(
                fname.endswith(sfx) for sfx in (
                    ".tsv", ".tsv.gz", ".csv", ".csv.gz",
                    ".txt", ".txt.gz", ".tsv", ".tsv.gz",
                )
            )
            if not is_text_ext and ext not in ("", ".gz"):
                continue

            abs_path = str(Path(data_dir) / sc.path)
            line_count = self._count_lines_in_file(abs_path)
            if line_count > 0:
                dir_candidates.append((fid, sc, line_count))

        if not dir_candidates:
            return None

        # Match by dimension
        candidates_A = [(fid, sc, lc) for fid, sc, lc in dir_candidates if lc == A]
        candidates_B = [(fid, sc, lc) for fid, sc, lc in dir_candidates if lc == B]

        def _name_hints_barcodes(name: str) -> bool:
            n = name.lower()
            return any(kw in n for kw in ("barc", "cell", "meta"))

        def _name_hints_features(name: str) -> bool:
            n = name.lower()
            return any(kw in n for kw in ("feat", "gene"))

        barcodes_fid: str | None = None
        features_fid: str | None = None
        barcodes_sc: FileSidecar | None = None
        features_sc: FileSidecar | None = None

        if A != B:
            # Try heuristic-based identification
            barc_A = [(f, s, _) for f, s, _ in candidates_A if _name_hints_barcodes(s.filename)]
            feat_B = [(f, s, _) for f, s, _ in candidates_B if _name_hints_features(s.filename)]
            feat_A = [(f, s, _) for f, s, _ in candidates_A if _name_hints_features(s.filename)]
            barc_B = [(f, s, _) for f, s, _ in candidates_B if _name_hints_barcodes(s.filename)]

            if barc_A and feat_B:
                # A-dimension files are barcodes, B-dimension files are features
                barcodes_fid, barcodes_sc, _ = barc_A[0]
                features_fid, features_sc, _ = feat_B[0]
            elif feat_A and barc_B:
                # A-dimension files are features, B-dimension files are barcodes
                features_fid, features_sc, _ = feat_A[0]
                barcodes_fid, barcodes_sc, _ = barc_B[0]
            elif candidates_A and candidates_B:
                # Fallback: larger dimension is typically genes (features),
                # smaller dimension is typically cells (barcodes)
                if A > B:
                    features_fid, features_sc, _ = candidates_A[0]
                    barcodes_fid, barcodes_sc, _ = candidates_B[0]
                else:
                    barcodes_fid, barcodes_sc, _ = candidates_A[0]
                    features_fid, features_sc, _ = candidates_B[0]

        if not barcodes_fid or not features_fid:
            self._logger.debug("Shape-matching could not resolve MEX companions",
                               matrix=mtx_sidecar.filename, dims=f"{A}x{B}")
            return None

        # Build the group
        confidence = min(
            mtx_sidecar.confidence,
            barcodes_sc.confidence,
            features_sc.confidence,
        )
        # Slightly reduce confidence for shape-matched groups vs filename-matched
        confidence = round(confidence * 0.95, 2)

        group = FileGroup(
            group_id=f"mex_shape_{matrix_fid[:12]}",
            group_type=FileGroupType.MEX_10X,
            sample_id=mtx_sidecar.sample_id_guess,
            gsm=mtx_sidecar.gsm,
            files={
                "matrix": matrix_fid,
                "barcodes": barcodes_fid,
                "features": features_fid,
            },
            route="read_10x_mtx_to_h5ad",
            route_description="Convert 10x MEX triplet to h5ad (shape-matched)",
            confidence=confidence,
            notes=[
                f"Matrix: {mtx_sidecar.filename} ({A}x{B})",
                f"Barcodes: {barcodes_sc.filename} (shape-matched, {B if A > B else A} lines)",
                f"Features: {features_sc.filename} (shape-matched, {A if A > B else B} lines)",
                "Grouped via mathematical shape-matching (non-standard filenames)",
            ],
        )
        assigned.update({matrix_fid, barcodes_fid, features_fid})
        self._logger.info("Shape-matched MEX triplet",
                          matrix=mtx_sidecar.filename,
                          barcodes=barcodes_sc.filename,
                          features=features_sc.filename)
        return group
