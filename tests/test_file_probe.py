"""Tests for the file probe — format detection, archive handling, non-count filtering.

Tests cover:
- 10x archive filename patterns (ARCHIVE_10X detection)
- Non-count supplementary file filtering (BCR/TCR, peaks, diff results, GFF)
- Spatial transcriptomics detection
- Gzip content inspection for MEX components
- Gene annotation detection inside gzip
"""

from __future__ import annotations

from pathlib import Path

import pytest

from gse_router.probes.file_probe import FileProbe
from gse_router.schemas.sidecar import FileDetectedType


class TestArchive10xDetection:
    """Test detection of archived 10x Genomics bundles via content peeking.

    Creates real tar.gz archives containing MEX triplets so the
    content-peeking logic (_peek_archive_contents + _classify_archive_by_contents)
    can detect the 10x MEX bundle inside.
    """

    def _create_tar_gz_with_mex(self, tmp_path, archive_name):
        """Create a tar.gz archive containing MEX triplet files."""
        import gzip
        import io
        import tarfile

        archive_path = tmp_path / archive_name

        # Build the tar in memory
        tar_buffer = io.BytesIO()
        with tarfile.open(fileobj=tar_buffer, mode="w") as tf:
            # matrix.mtx
            mtx_info = tarfile.TarInfo(name="matrix.mtx")
            mtx_content = b"%%MatrixMarket matrix coordinate integer general\n%comment\n10 5 50\n1 1 10\n"
            mtx_info.size = len(mtx_content)
            tf.addfile(mtx_info, io.BytesIO(mtx_content))

            # barcodes.tsv
            bc_info = tarfile.TarInfo(name="barcodes.tsv")
            bc_content = b"AAACCTGAGCCTA\nAAACCTGAGGTCA\n"
            bc_info.size = len(bc_content)
            tf.addfile(bc_info, io.BytesIO(bc_content))

            # features.tsv
            ft_info = tarfile.TarInfo(name="features.tsv")
            ft_content = b"ENSG000001\tTP53\tGene Expression\n"
            ft_info.size = len(ft_content)
            tf.addfile(ft_info, io.BytesIO(ft_content))

        # Gzip-compress the tar
        tar_bytes = tar_buffer.getvalue()
        with gzip.open(archive_path, "wb") as gz_fh:
            gz_fh.write(tar_bytes)

        return archive_path

    def test_tar_gz_with_mex_triplet_detected(self, tmp_path):
        """A tar.gz containing a MEX triplet should be classified as ARCHIVE with 10x evidence."""
        probe = FileProbe()
        file_path = self._create_tar_gz_with_mex(tmp_path, "sample_data.tar.gz")

        sidecar = probe.probe_file(file_path, relative_root=tmp_path)
        assert sidecar is not None
        assert sidecar.detected_type == FileDetectedType.ARCHIVE
        # Content peeking should find the MEX triplet
        assert any("10x MEX triplet" in ev for ev in sidecar.evidence), \
            f"Expected '10x MEX triplet' evidence, got: {sidecar.evidence}"

    def test_tar_gz_with_mex_nonstandard_name(self, tmp_path):
        """Even with non-standard names, MEX contents should be detected by peeking."""
        probe = FileProbe()
        file_path = self._create_tar_gz_with_mex(tmp_path, "GSE123_RAW.tar.gz")

        sidecar = probe.probe_file(file_path, relative_root=tmp_path)
        assert sidecar is not None
        assert sidecar.detected_type == FileDetectedType.ARCHIVE
        assert any("10x MEX triplet" in ev for ev in sidecar.evidence), \
            f"Content peeking should detect MEX regardless of archive name. Evidence: {sidecar.evidence}"

    def test_zip_archive_with_mex_files(self, tmp_path):
        """ZIP archive containing MEX files should be detected via content peeking."""
        import zipfile
        probe = FileProbe()

        file_path = tmp_path / "sample_data.zip"
        with zipfile.ZipFile(file_path, "w") as zf:
            zf.writestr("matrix.mtx", "%%MatrixMarket\n10 5 50\n")
            zf.writestr("barcodes.tsv", "AAACCTGAGCCTA\n")
            zf.writestr("features.tsv", "TP53\n")

        sidecar = probe.probe_file(file_path, relative_root=tmp_path)
        assert sidecar is not None
        assert sidecar.detected_type == FileDetectedType.ARCHIVE
        assert any("10x MEX triplet" in ev for ev in sidecar.evidence), \
            f"Expected '10x MEX triplet' evidence in ZIP, got: {sidecar.evidence}"

    def test_archive_with_h5ad_contents(self, tmp_path):
        """Archive containing .h5ad files should be detected and noted in evidence."""
        import zipfile
        probe = FileProbe()

        file_path = tmp_path / "processed_data.zip"
        with zipfile.ZipFile(file_path, "w") as zf:
            zf.writestr("sample1.h5ad", "dummy h5ad content")
            zf.writestr("sample2.h5ad", "dummy h5ad content")

        sidecar = probe.probe_file(file_path, relative_root=tmp_path)
        assert sidecar is not None
        assert sidecar.detected_type == FileDetectedType.ARCHIVE
        assert any(".h5ad" in ev for ev in sidecar.evidence), \
            f"Expected h5ad evidence, got: {sidecar.evidence}"


class TestNonCountSupplementaryFiltering:
    """Test that non-count supplementary files are NOT routed as count matrices."""

    @pytest.mark.parametrize("filename,expected_type", [
        # Immune receptor data (BCR/TCR)
        ("filtered_contig_annotations_bcr.csv.gz", FileDetectedType.IMMUNE_RECEPTOR_DATA),
        ("filtered_contig_annotations_tcr.csv.gz", FileDetectedType.IMMUNE_RECEPTOR_DATA),
        ("contig_annotations.csv.gz", FileDetectedType.IMMUNE_RECEPTOR_DATA),
        ("vdj_b_contig_annotations.csv", FileDetectedType.IMMUNE_RECEPTOR_DATA),
        ("clonotypes.csv", FileDetectedType.IMMUNE_RECEPTOR_DATA),
        # Peak data (ATAC/ChIP-seq)
        ("peaks.bed.gz", FileDetectedType.PEAK_DATA),
        ("macs2_peaks.xls.gz", FileDetectedType.PEAK_DATA),
        ("sample_peaks.bed", FileDetectedType.PEAK_DATA),
        ("chip_peaks.narrowPeak", FileDetectedType.PEAK_DATA),
        ("atac_peaks.broadPeak", FileDetectedType.PEAK_DATA),
        # Diff results
        ("deseq2_results.csv.gz", FileDetectedType.DIFF_RESULTS),
        ("diffbind_result.csv", FileDetectedType.DIFF_RESULTS),
        ("differential_expression.csv.gz", FileDetectedType.DIFF_RESULTS),
        ("deg_wt_vs_ko.csv", FileDetectedType.DIFF_RESULTS),
    ])
    def test_non_count_files_detected(self, filename, expected_type, tmp_path):
        """Verify non-count supplementary files are detected correctly."""
        probe = FileProbe()
        file_path = tmp_path / filename

        # Create a file with appropriate content
        if filename.endswith(".gz"):
            import gzip
            with gzip.open(file_path, "wb") as fh:
                fh.write(b"header1\theader2\nvalue1\tvalue2\n")
        else:
            file_path.write_text("header1\theader2\nvalue1\tvalue2\n")

        sidecar = probe.probe_file(file_path, relative_root=tmp_path)
        assert sidecar is not None
        assert sidecar.detected_type == expected_type, \
            f"Expected {expected_type} for '{filename}', got {sidecar.detected_type}"

    def test_non_count_files_skip_count_routing(self, tmp_path):
        """Verify non-count files are mapped to UNKNOWN group type (never count matrix)."""
        from gse_router.router.classifier import FileClassifier
        from gse_router.schemas.manifest import FileGroupType

        probe = FileProbe()
        classifier = FileClassifier()

        non_count_files = [
            "filtered_contig_annotations_bcr.csv",
            "peaks.bed",
            "deseq2_results.csv",
        ]

        sidecars = {}
        for filename in non_count_files:
            file_path = tmp_path / filename
            file_path.write_text("col1\tcol2\nval1\tval2\n")
            sc = probe.probe_file(file_path, relative_root=tmp_path)
            sidecars[sc.file_id] = sc

        groups = classifier.group_files(sidecars)
        for g in groups:
            assert g.group_type != FileGroupType.CSV_MATRIX, \
                f"Non-count file mapped to CSV_MATRIX: {g}"
            assert g.group_type != FileGroupType.RAW_UMI_MATRIX, \
                f"Non-count file mapped to RAW_UMI_MATRIX: {g}"


class TestGFFAnnotationInGzip:
    """Test detection of gene annotations inside gzip."""

    def test_gff_gz_detected_as_annotation(self, tmp_path):
        """GFF files inside .gz should be detected as gene annotations, not generic matrices."""
        probe = FileProbe()
        file_path = tmp_path / "annotations.gff.gz"

        import gzip
        with gzip.open(file_path, "wb") as fh:
            fh.write(b"##gff-version 3\nchr1\tsource\tgene\t1\t1000\t.\t+\t.\tID=gene1\n")

        sidecar = probe.probe_file(file_path, relative_root=tmp_path)
        assert sidecar is not None
        assert sidecar.detected_type == FileDetectedType.GENE_ANNOTATION, \
            f"Expected GENE_ANNOTATION, got {sidecar.detected_type}"

    def test_gtf_gz_detected_as_annotation(self, tmp_path):
        """GTF files inside .gz should be detected as gene annotations."""
        probe = FileProbe()
        file_path = tmp_path / "genes.gtf.gz"

        import gzip
        with gzip.open(file_path, "wb") as fh:
            fh.write(b"chr1\tensembl\tgene\t1\t1000\t.\t+\t.\tgene_id \"ENSG1\";\n")

        sidecar = probe.probe_file(file_path, relative_root=tmp_path)
        assert sidecar is not None
        assert sidecar.detected_type == FileDetectedType.GENE_ANNOTATION, \
            f"Expected GENE_ANNOTATION, got {sidecar.detected_type}"

    def test_gff3_gz_detected_as_annotation(self, tmp_path):
        """GFF3 files inside .gz should be detected as gene annotations."""
        probe = FileProbe()
        file_path = tmp_path / "features.gff3.gz"

        import gzip
        with gzip.open(file_path, "wb") as fh:
            fh.write(b"##gff-version 3\nchr1\t.\tgene\t1\t500\t.\t+\t.\tID=g1\n")

        sidecar = probe.probe_file(file_path, relative_root=tmp_path)
        assert sidecar is not None
        assert sidecar.detected_type == FileDetectedType.GENE_ANNOTATION, \
            f"Expected GENE_ANNOTATION, got {sidecar.detected_type}"

    def test_plain_gff_detected(self, tmp_path):
        """Plain .gff files should be detected as gene annotations."""
        probe = FileProbe()
        file_path = tmp_path / "genome.gff"
        file_path.write_text("##gff-version 3\nchr1\t.\tgene\t1\t1000\t.\t+\t.\tID=g1\n")

        sidecar = probe.probe_file(file_path, relative_root=tmp_path)
        assert sidecar is not None
        assert sidecar.detected_type == FileDetectedType.GENE_ANNOTATION


class TestSpatialTranscriptomics:
    """Test detection of spatial transcriptomics data."""

    @pytest.mark.parametrize("filename", [
        "tissue_positions_list.csv",
        "scalefactors_json.json",
        "tissue_hires_image.png",
        "tissue_lowres_image.png",
        "aligned_fiducials.jpg",
        "detected_tissue_image.jpg",
        "spatial_barcodes.tsv",
    ])
    def test_spatial_files_not_count_matrix(self, filename, tmp_path):
        """Spatial transcriptomics files should not be classified as count matrices."""
        probe = FileProbe()
        file_path = tmp_path / filename

        if filename.endswith((".csv", ".tsv", ".json")):
            file_path.write_text("dummy\tcontent\n")
        else:
            file_path.write_bytes(b"\x89PNG\r\n\x1a\n" + b"\x00" * 100)

        sidecar = probe.probe_file(file_path, relative_root=tmp_path)
        assert sidecar is not None
        # Should NOT be classified as a count matrix type
        assert sidecar.detected_type not in (
            FileDetectedType.CSV_MATRIX,
            FileDetectedType.TSV_MATRIX,
            FileDetectedType.RAW_UMI_MATRIX,
            FileDetectedType.NORMALIZED_MATRIX,
        ), f"Spatial file '{filename}' should not be classified as count matrix"


class TestContentBasedArchiveClassification:
    """Test content-based archive classification (replaces filename heuristics)."""

    def test_empty_archive_fallback(self, tmp_path):
        """Empty or unreadable archives return None from peek."""
        # A file with no archive magic bytes should return None
        probe = FileProbe()

        file_path = tmp_path / "not_an_archive.txt"
        file_path.write_text("just some text, not an archive")

        contents = FileProbe._peek_archive_contents(file_path)
        assert contents is None, \
            "Non-archive files should return None from _peek_archive_contents"

    def test_unrecognized_contents_confidence(self, tmp_path):
        """Archives with only text files get moderate confidence from content peeking."""
        import zipfile
        probe = FileProbe()

        file_path = tmp_path / "mystery.zip"
        with zipfile.ZipFile(file_path, "w") as zf:
            zf.writestr("README.txt", "no data here")
            zf.writestr("notes.md", "just notes")

        contents = FileProbe._peek_archive_contents(file_path)
        assert contents is not None
        assert "README.txt" in contents

        dtype, conf, evidence = FileProbe._classify_archive_by_contents(contents)
        assert dtype == FileDetectedType.ARCHIVE
        # .txt files trigger the delimited-text heuristic → moderate confidence
        assert 0.65 <= conf <= 0.80, \
            f"Expected confidence 0.65-0.80 for delimited text archive, got {conf}"

    def test_unknown_contents_fallback(self, tmp_path):
        """Archives with only non-data files get low confidence."""
        import zipfile
        probe = FileProbe()

        file_path = tmp_path / "docs.zip"
        with zipfile.ZipFile(file_path, "w") as zf:
            zf.writestr("methods.pdf", "pdf content")
            zf.writestr("supplementary.docx", "docx content")

        contents = FileProbe._peek_archive_contents(file_path)
        assert contents is not None
        dtype, conf, evidence = FileProbe._classify_archive_by_contents(contents)
        assert dtype == FileDetectedType.ARCHIVE
        # No recognizable scRNA-seq data files → lower confidence
        assert conf <= 0.65, \
            f"Expected low confidence for unrecognized archive, got {conf}"
