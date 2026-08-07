"""Enhanced batch tester — downloads files, runs content-first probing, deletes data.

For each GSE dataset:
1. Fetch metadata + list supplementary files
2. Download the smallest informative file (or peek inside archives)
3. Run content-first FileProbe on it
4. Compare: what would filename-only say vs content-first?
5. Record results, delete data, move to next

This actually exercises our content-first improvements (h5py HDF5 inspection,
gzip content-first restructuring, archive content peeking).
"""
from __future__ import annotations

import asyncio
import json
import shutil
import sys
import traceback
from datetime import datetime, timezone
from pathlib import Path

# Add project root to path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

# Load .env before anything else (needed for DEEPSEEK_API_KEY, NCBI_API_KEY, etc.)
from gse_router.config import get_config
get_config()  # triggers _load_dotenv()

from gse_router.probes.file_probe import FileProbe
from gse_router.probes.metadata import MetadataProbe
from gse_router.tools.geo_downloader import GeoSupplementaryLister, GeoFileDownloader

DATA_DIR = Path("E:/project/Aging_OS2/data")
RESULTS_FILE = Path("E:/project/Aging_OS2/content_probe_results.jsonl")
ERROR_LOG = Path("E:/project/Aging_OS2/content_probe_errors.log")

GSE_LIST = [
    "GSE337655", "GSE337289", "GSE327907", "GSE327389", "GSE325198",
    "GSE301966", "GSE276559", "GSE303713", "GSE337368", "GSE326370",
    "GSE309120", "GSE301717", "GSE336843", "GSE327186", "GSE319066",
    "GSE318800", "GSE318657", "GSE311869", "GSE308872", "GSE301594",
    "GSE248386", "GSE248385", "GSE337384", "GSE336861", "GSE336668",
    "GSE327808", "GSE319693", "GSE319592", "GSE300504", "GSE299068",
    "GSE271348", "GSE337230", "GSE337208", "GSE337162", "GSE336973",
    "GSE336400", "GSE335941", "GSE335811", "GSE332958", "GSE329941",
    "GSE329924", "GSE328525", "GSE328097", "GSE327732", "GSE322748",
    "GSE322527", "GSE320454", "GSE317299", "GSE314969", "GSE314267",
]

# Files that are never useful for format probing
SKIP_PATTERNS = ["filelist.txt", "README", "readme", "_feature_reference",
                 "_feature_README", "metadata.tsv", "metadata.csv",
                 "cell_metadata", "sample_metadata", "demuxalot"]


def log_result(entry: dict) -> None:
    with open(RESULTS_FILE, "a") as f:
        f.write(json.dumps(entry) + "\n")
    status = entry.get("status", "???")
    gse = entry.get("gse_id", "???")
    detected = entry.get("content_detected_type", "????")
    fn_type = entry.get("filename_detected_type", "????")
    match = "MATCH" if entry.get("content_matches_filename", False) else "DIFF"
    confidence = entry.get("content_confidence", 0)
    print(f"  [{status:>14s}] {gse}  content={detected:>25s}({confidence:.0f})  "
          f"filename={fn_type:>25s}  {match}")


def log_error(gse_id: str, stage: str, error: str) -> None:
    ts = datetime.now(timezone.utc).isoformat()
    with open(ERROR_LOG, "a") as f:
        f.write(f"[{ts}] {gse_id} @ {stage}\n{error}\n---\n")


def cleanup(gse_id: str) -> None:
    """Remove all downloaded data for this GSE."""
    paths = [
        DATA_DIR / gse_id,
        DATA_DIR / "supplementary" / gse_id,
    ]
    for p in paths:
        if p.exists():
            try:
                if p.is_dir():
                    shutil.rmtree(p)
                else:
                    p.unlink()
            except Exception:
                pass
    # Also clean any _extracted_ dirs
    if DATA_DIR.exists():
        for d in DATA_DIR.iterdir():
            if d.is_dir() and d.name.startswith("_extracted_"):
                try:
                    shutil.rmtree(d)
                except Exception:
                    pass


def _should_skip_file(filename: str) -> bool:
    """Check if this file is a non-data file we should skip for probing."""
    name_lower = filename.lower()
    for pat in SKIP_PATTERNS:
        if pat.lower() in name_lower:
            return True
    return False


def _pick_best_file(files: list, max_disk_mb: int = 500) -> tuple | None:
    """Pick the best file to download for content probing.

    Priority (within size constraints):
    1. Small .h5ad files — tests HDF5 content inspection (h5py groups)
    2. Small .h5 files — tests HDF5 content inspection for 10x/loom
    3. Small .mtx.gz files — tests gzip MEX detection
    4. Small .csv.gz/.tsv.gz files — tests delimited text probing
    5. Small .loom files — tests HDF5 with /matrix + /layers
    6. .tar.gz / .zip archives — tests archive content peeking
    7. .rds files — tests RDS magic byte detection
    8. Unknown .gz — probe to determine content

    All limits are relative to max_disk_mb. We delete files immediately
    after probing, so the peak disk usage is just one file at a time.

    Returns (file_info, reason) or None.
    """
    candidates = []

    for f in files:
        name = f.filename.lower()
        size_mb = (f.size_bytes or 0) / 1e6

        if _should_skip_file(f.filename):
            continue
        if size_mb <= 0:
            continue
        if size_mb > max_disk_mb:
            continue

        if name.endswith(".h5ad"):
            candidates.append((0, size_mb, f, "h5ad_hdf5_inspect"))
        elif name.endswith(".h5") and not name.endswith(".h5ad"):
            candidates.append((1, size_mb, f, "h5_hdf5_inspect"))
        elif name.endswith(".h5seurat"):
            candidates.append((2, size_mb, f, "h5seurat_hdf5_inspect"))
        elif any(name.endswith(ext) for ext in [".mtx.gz", ".mtx"]):
            candidates.append((3, size_mb, f, "mtx_gzip_inspect"))
        elif any(name.endswith(ext) for ext in [".csv.gz", ".tsv.gz", ".txt.gz",
                                                  ".csv", ".tsv", ".txt",
                                                  ".tab.gz", ".tab"]):
            candidates.append((4, size_mb, f, "delimited_text_inspect"))
        elif any(name.endswith(ext) for ext in [".loom", ".loom.gz"]):
            candidates.append((5, size_mb, f, "loom_hdf5_inspect"))
        elif any(name.endswith(ext) for ext in [".tar.gz", ".tar", ".tgz", ".zip"]):
            candidates.append((6, size_mb, f, "archive_peek"))
        elif any(name.endswith(ext) for ext in [".rds", ".rdata", ".rda", ".RDS"]):
            candidates.append((7, size_mb, f, "rds_magic_inspect"))
        elif name.endswith(".gz"):
            # Unknown gzip — probe it to see what's inside
            candidates.append((8, size_mb, f, "unknown_gzip_inspect"))
        elif name.endswith(".mat"):
            candidates.append((9, size_mb, f, "matlab_inspect"))

    if not candidates:
        return None

    # Sort by priority, then by size (prefer smaller within same priority)
    candidates.sort(key=lambda x: (x[0], x[1]))
    priority, size_mb, file_info, reason = candidates[0]
    return (file_info, reason)


def _classify_filename_only(filename: str) -> str:
    """What would the filename alone suggest? (Pre-content-first behavior.)"""
    name = filename.lower()
    # h5ad
    if name.endswith(".h5ad") or name.endswith(".h5ad.gz"):
        return "H5AD"
    # 10x H5
    if "filtered_feature_bc_matrix" in name or "raw_feature_bc_matrix" in name:
        return "H5_10X"
    if name.endswith(".h5seurat"):
        return "H5SEURAT"
    # loom
    if ".loom" in name:
        return "LOOM"
    # MEX components
    if "matrix.mtx" in name:
        return "MEX_MATRIX"
    if "barcodes.tsv" in name:
        return "MEX_BARCODES"
    if "features.tsv" in name or "genes.tsv" in name:
        return "MEX_FEATURES"
    if name.endswith(".mtx.gz") or name.endswith(".mtx"):
        return "MTX_GENERIC"
    # RDS
    if name.endswith(".rds") or name.endswith(".rdata") or name.endswith(".rda"):
        return "RDS"
    # Generic H5
    if name.endswith(".h5") or name.endswith(".h5.gz"):
        return "H5_GENERIC"
    # Archive
    if any(name.endswith(ext) for ext in [".tar.gz", ".tar", ".tgz", ".zip"]):
        if any(pat in name for pat in ["_bc_matrix", ".matrix.", "_gene_bc_matrices",
                                         "_feature_bc_matrix", "cellranger", "outs"]):
            return "ARCHIVE_10X_FILENAME"
        return "ARCHIVE"
    # Delimited
    if any(name.endswith(ext) for ext in [".csv.gz", ".csv", ".tsv.gz", ".tsv",
                                            ".txt.gz", ".txt", ".tab.gz", ".tab"]):
        if "tpm" in name or "fpkm" in name or "rpkm" in name or "normalized" in name:
            return "NORMALIZED_MATRIX"
        if "count" in name or "raw" in name:
            return "RAW_COUNTS_MATRIX"
        return "DELIMITED_MATRIX"
    # Immune receptor
    if any(pat in name for pat in ["contig_annotation", "clonotype", "vdj"]):
        return "IMMUNE_RECEPTOR"
    # Peak
    if any(pat in name for pat in ["peaks", "narrowpeak", "broadpeak"]):
        return "PEAK_DATA"
    # Diff
    if any(pat in name for pat in ["deseq2", "diffbind", "differential_expression", "deg_"]):
        return "DIFF_RESULTS"
    # GFF/GTF
    if any(pat in name for pat in [".gff", ".gtf", ".gff3"]):
        return "GENE_ANNOTATION"
    # Spatial
    if any(pat in name for pat in ["tissue_positions", "scalefactors", "tissue_hires",
                                     "tissue_lowres", "spatial_barcodes"]):
        return "SPATIAL"
    # Gzip unknown
    if name.endswith(".gz"):
        return "GZIP_UNKNOWN"
    return "UNKNOWN"


async def probe_one_dataset(gse_id: str, probe: FileProbe) -> dict:
    """Download one file from a dataset and run content-first probing on it."""
    result = {
        "gse_id": gse_id,
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "status": "skipped",
        "downloaded_file": "",
        "downloaded_size_mb": 0,
        "content_detected_type": "n/a",
        "content_confidence": 0.0,
        "content_evidence": [],
        "filename_detected_type": "n/a",
        "content_matches_filename": None,
        "organism": [],
        "technology": "",
        "supplementary_file_count": 0,
        "error": "",
    }

    # --- Phase 1: Metadata ---
    try:
        meta = MetadataProbe()
        ctx = await meta.probe_all(gse_id)
        result["organism"] = ctx.organism
        result["technology"] = ctx.technology
    except Exception as e:
        result["error"] = f"Metadata failed: {e}"
        result["status"] = "metadata_failed"
        return result

    # --- Phase 2: List supplementary files ---
    try:
        lister = GeoSupplementaryLister()
        files = await lister.list_supplementary(gse_id)
        await lister.close()
        result["supplementary_file_count"] = len(files)
    except Exception as e:
        result["error"] = f"List failed: {e}"
        result["status"] = "list_failed"
        return result

    if not files:
        result["status"] = "no_files"
        return result

    # --- Phase 3: Pick a file and download ---
    pick = _pick_best_file(files)
    if pick is None:
        # All files too large or non-probable
        total_files = [f.filename for f in files if not _should_skip_file(f.filename)]
        if total_files:
            # List what files are available for logging
            result["status"] = "all_too_large"
            result["error"] = f"Files too large: {[f'{f.filename}({(f.size_bytes or 0)/1e6:.1f}MB)' for f in files if not _should_skip_file(f.filename)][:5]}"
        else:
            result["status"] = "no_probable_files"
        return result

    file_info, reason = pick
    size_mb = (file_info.size_bytes or 0) / 1e6
    result["downloaded_file"] = file_info.filename
    result["downloaded_size_mb"] = round(size_mb, 2)
    result["pick_reason"] = reason
    result["filename_detected_type"] = _classify_filename_only(file_info.filename)

    # --- Phase 4: Download ---
    dest_dir = DATA_DIR / gse_id
    dest_dir.mkdir(parents=True, exist_ok=True)
    dest_path = dest_dir / file_info.filename

    try:
        downloader = GeoFileDownloader(data_dir=DATA_DIR, max_concurrent=1)
        await downloader.download_file(
            file_info.url,
            dest_path,
            expected_size=file_info.size_bytes,
        )
        await downloader.close()
    except Exception as e:
        result["error"] = f"Download failed: {e}"
        result["status"] = "download_failed"
        log_error(gse_id, "download", traceback.format_exc())
        return result

    if not dest_path.exists():
        result["error"] = "File not found after download"
        result["status"] = "download_failed"
        return result

    # --- Phase 5: Content-first probing ---
    try:
        sidecar = probe.probe_file(dest_path, relative_root=dest_dir)
        if sidecar:
            result["content_detected_type"] = sidecar.detected_type.value
            result["content_confidence"] = sidecar.confidence
            result["content_evidence"] = sidecar.evidence[:10]  # cap for log size
            result["content_stage"] = sidecar.processing_stage.value

            # Compare with filename-only
            fn_type = result["filename_detected_type"]
            ct_type = sidecar.detected_type.value
            # Broad comparison: do they point to same general category?
            result["content_matches_filename"] = _types_compatible(fn_type, ct_type)
        else:
            result["content_detected_type"] = "NONE"
            result["content_confidence"] = 0.0
            result["content_matches_filename"] = False
            result["error"] = "probe_file returned None"
    except Exception as e:
        result["error"] = f"Probe failed: {e}"
        result["status"] = "probe_failed"
        log_error(gse_id, "probe_file", traceback.format_exc())
        return result

    # --- Determine status ---
    if result["content_detected_type"] in ("UNKNOWN", "NONE"):
        result["status"] = "probe_unknown"
    elif result["content_confidence"] < 0.5:
        result["status"] = "probe_low_confidence"
    elif result["content_matches_filename"]:
        result["status"] = "content_confirms_filename"
    else:
        result["status"] = "content_differs_from_filename"

    return result


def _types_compatible(fn_type: str, ct_type: str) -> bool:
    """Check if filename-based type and content-based type are broadly compatible."""
    # Map both to broad categories
    H5AD_LIKE = {"H5AD", "h5ad", "h5ad/h5seurat"}
    H5_LIKE = {"H5_10X", "H5_GENERIC", "h5_generic", "H5SEURAT", "h5"}
    MEX_LIKE = {"MEX_MATRIX", "MEX_BARCODES", "MEX_FEATURES", "10x_mex", "MTX_GENERIC"}
    ARCHIVE_LIKE = {"ARCHIVE", "ARCHIVE_10X_FILENAME", "archive", "archive_10x"}
    MATRIX_LIKE = {"DELIMITED_MATRIX", "RAW_COUNTS_MATRIX", "NORMALIZED_MATRIX",
                   "delimited_matrix", "CSV_MATRIX", "TSV_MATRIX",
                   "RAW_UMI_MATRIX", "NORMALIZED_MATRIX"}
    RDS_LIKE = {"RDS", "rds"}
    LOOM_LIKE = {"LOOM", "loom"}
    IMMUNE_LIKE = {"IMMUNE_RECEPTOR", "IMMUNE_RECEPTOR_DATA"}
    PEAK_LIKE = {"PEAK_DATA", "peak_data"}
    DIFF_LIKE = {"DIFF_RESULTS", "diff_results"}
    ANNOT_LIKE = {"GENE_ANNOTATION", "gene_annotation"}
    SPATIAL_LIKE = {"SPATIAL", "spatial"}

    def _category(t: str) -> str:
        t_upper = t.upper()
        if t in H5AD_LIKE or t_upper in H5AD_LIKE or "H5AD" in t_upper:
            return "h5ad"
        if t in H5_LIKE or t_upper in {"H5_10X", "H5_GENERIC", "H5SEURAT"}:
            return "h5"
        if t in MEX_LIKE or t_upper in {"MEX_MATRIX", "MEX_BARCODES", "MEX_FEATURES", "MTX_GENERIC"}:
            return "mex"
        if t in ARCHIVE_LIKE or "ARCHIVE" in t_upper:
            return "archive"
        if t in MATRIX_LIKE or "MATRIX" in t_upper:
            return "matrix"
        if t in RDS_LIKE or "RDS" in t_upper:
            return "rds"
        if t in LOOM_LIKE or "LOOM" in t_upper:
            return "loom"
        if t in IMMUNE_LIKE or "IMMUNE" in t_upper:
            return "immune"
        if t in PEAK_LIKE or "PEAK" in t_upper:
            return "peak"
        if t in DIFF_LIKE or "DIFF" in t_upper:
            return "diff"
        if t in ANNOT_LIKE or "ANNOTATION" in t_upper or "GFF" in t_upper or "GTF" in t_upper:
            return "annotation"
        if t in SPATIAL_LIKE or "SPATIAL" in t_upper:
            return "spatial"
        if t in ("GZIP_UNKNOWN", "UNKNOWN", "NONE", "n/a"):
            return "unknown"
        return "other"

    return _category(fn_type) == _category(ct_type)


async def main():
    probe = FileProbe()

    # Skip already-processed
    already_done = set()
    if RESULTS_FILE.exists():
        with open(RESULTS_FILE) as f:
            for line in f:
                try:
                    entry = json.loads(line.strip())
                    already_done.add(entry.get("gse_id", ""))
                except Exception:
                    pass

    pending = [g for g in GSE_LIST if g not in already_done]
    print(f"=== Content-First Probe Batch Test ===")
    print(f"{len(pending)} remaining of {len(GSE_LIST)} total\n")
    print(f"{'#':>4s}  {'GSE':>11s}  {'Status':>20s}  "
          f"{'Content Type':>28s}  {'Conf':>5s}  {'Filename Type':>28s}  Match")
    print("-" * 140)

    for i, gse_id in enumerate(pending, 1):
        print(f"\n[{i:>3d}/{len(pending)}] {gse_id} — probing...", end="", flush=True)
        try:
            result = await probe_one_dataset(gse_id, probe)
        except Exception as e:
            result = {
                "gse_id": gse_id,
                "timestamp": datetime.now(timezone.utc).isoformat(),
                "status": "crash",
                "error": f"{type(e).__name__}: {e}",
            }
            log_error(gse_id, "main_loop", traceback.format_exc())

        log_result(result)
        cleanup(gse_id)

    # Summary
    print(f"\n=== Done. Results: {RESULTS_FILE} ===")
    _print_summary()


def _print_summary():
    if not RESULTS_FILE.exists():
        return
    results = []
    with open(RESULTS_FILE) as f:
        for line in f:
            try:
                results.append(json.loads(line.strip()))
            except Exception:
                pass

    total = len(results)
    content_diff = sum(1 for r in results if r.get("status") == "content_differs_from_filename")
    content_confirm = sum(1 for r in results if r.get("status") == "content_confirms_filename")
    probe_unknown = sum(1 for r in results if r.get("status") == "probe_unknown")
    low_conf = sum(1 for r in results if r.get("status") == "probe_low_confidence")
    skipped = sum(1 for r in results if r.get("status") in ("skipped", "all_too_large", "no_probable_files"))
    failed = sum(1 for r in results if r.get("status") in ("download_failed", "probe_failed", "metadata_failed", "list_failed", "crash"))

    print(f"\n=== Summary ({total} datasets) ===")
    print(f"  Content confirms filename:    {content_confirm:>3d} ({content_confirm/total*100:.0f}%)")
    print(f"  Content DIFFERS from filename: {content_diff:>3d} ({content_diff/total*100:.0f}%)")
    print(f"  Probe returned UNKNOWN:        {probe_unknown:>3d} ({probe_unknown/total*100:.0f}%)")
    print(f"  Probe low confidence (<0.5):   {low_conf:>3d} ({low_conf/total*100:.0f}%)")
    print(f"  Skipped (no probable file):    {skipped:>3d} ({skipped/total*100:.0f}%)")
    print(f"  Failed (download/probe err):   {failed:>3d} ({failed/total*100:.0f}%)")
    print(f"\nTotal content-first improvements: {content_diff} datasets where content "
          f"detection would route differently than filename.")


if __name__ == "__main__":
    asyncio.run(main())
