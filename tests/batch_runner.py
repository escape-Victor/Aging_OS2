"""Batch test 50 GEO datasets — one at a time, log results, delete data after each."""
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

from gse_router.pipeline.orchestrator import Orchestrator
from gse_router.probes.metadata import MetadataProbe

DATA_DIR = Path("E:/project/Aging_OS2/data")
OUTPUT_DIR = Path("E:/project/Aging_OS2/output")
RESULTS_FILE = Path("E:/project/Aging_OS2/batch_results.jsonl")
ERROR_LOG = Path("E:/project/Aging_OS2/batch_errors.log")

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


def log_result(entry: dict) -> None:
    with open(RESULTS_FILE, "a") as f:
        f.write(json.dumps(entry) + "\n")
    # Print summary line
    status = entry.get("status", "???")
    gse = entry.get("gse_id", "???")
    fmt = entry.get("best_format", "???")
    err = entry.get("error", "")
    ts = entry.get("files_downloaded", "?")
    print(f"  [{status:>12s}] {gse}  format={fmt:>20s}  downloaded={ts}  {err}")


def log_error(gse_id: str, stage: str, error: str) -> None:
    ts = datetime.now(timezone.utc).isoformat()
    with open(ERROR_LOG, "a") as f:
        f.write(f"[{ts}] {gse_id} @ {stage}\n{error}\n---\n")


def cleanup(gse_id: str) -> None:
    """Remove all downloaded/processed data for this GSE."""
    paths = [
        DATA_DIR / gse_id,
        DATA_DIR / "supplementary" / gse_id,
        OUTPUT_DIR / "h5ad",
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
    # Also clean any _extracted_ dirs in data dir
    if DATA_DIR.exists():
        for d in DATA_DIR.iterdir():
            if d.is_dir() and d.name.startswith("_extracted_"):
                try:
                    shutil.rmtree(d)
                except Exception:
                    pass


async def probe_only(gse_id: str) -> dict:
    """Lightweight probe: check metadata + supplementary file listing only.
    Returns summary without downloading actual data files."""
    result = {
        "gse_id": gse_id,
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "status": "unknown",
        "organism": [],
        "platform": [],
        "technology": "",
        "supplementary_file_count": 0,
        "supplementary_files": [],
        "supplementary_formats": [],
        "best_format": "none",
        "error": "",
    }

    probe = MetadataProbe()

    try:
        ctx = await probe.probe_all(gse_id)
        result["organism"] = ctx.organism
        result["platform"] = ctx.platform
        result["technology"] = ctx.technology
        result["sample_count"] = len(ctx.samples)
    except Exception as e:
        result["error"] = f"Metadata probe failed: {e}"
        log_error(gse_id, "metadata_probe", traceback.format_exc())
        result["status"] = "metadata_failed"
        return result

    # List supplementary files (no download)
    try:
        from gse_router.tools.geo_downloader import GeoSupplementaryLister

        lister = GeoSupplementaryLister()
        files = await lister.list_supplementary(gse_id)
        if files:
            result["supplementary_file_count"] = len(files)
            result["supplementary_files"] = [
                {"name": f.filename, "size_mb": round(f.size_bytes / 1e6, 2) if f.size_bytes else None}
                for f in files[:30]  # cap at 30 for log readability
            ]
            # Classify formats from filenames (skip filelist.txt)
            formats = set()
            for f in files:
                name = f.filename.lower()
                if name == "filelist.txt":
                    continue  # not a data file
                if any(p in name for p in [".h5ad", ".h5seurat"]):
                    formats.add("h5ad/h5seurat")
                elif any(p in name for p in ["filtered_feature_bc_matrix.h5", "raw_feature_bc_matrix.h5"]):
                    formats.add("10x_h5")
                elif any(p in name for p in ["matrix.mtx", "barcodes.tsv", "features.tsv", "genes.tsv"]):
                    formats.add("10x_mex")
                elif any(p in name for p in [".loom"]):
                    formats.add("loom")
                elif any(p in name for p in [".rds", ".rdata"]):
                    formats.add("rds")
                elif any(p in name for p in ["_bc_matrix.tar.gz", ".matrix.tar.gz"]):
                    formats.add("archive_10x")
                elif any(p in name for p in [".csv.gz", ".csv", ".tsv.gz", ".tsv", ".txt.gz", ".txt"]):
                    formats.add("delimited_matrix")
                elif any(p in name for p in [".tar.gz", ".tar", ".zip"]):
                    formats.add("archive")
                elif any(p in name for p in [".fastq", ".fq", ".sra"]):
                    formats.add("fastq/sra")
                elif any(p in name for p in [".h5"]):
                    formats.add("h5_generic")
                elif any(p in name for p in [".mat"]):
                    formats.add("matlab")
            result["supplementary_formats"] = sorted(formats)
    except Exception as e:
        result["supplementary_file_count"] = -1
        result["error"] = (result["error"] + "; " if result["error"] else "") + f"File listing failed: {e}"
        log_error(gse_id, "list_supplementary", traceback.format_exc())

    # Determine best format
    fmts = result["supplementary_formats"]
    if "h5ad/h5seurat" in fmts:
        result["best_format"] = "h5ad"
    elif "10x_h5" in fmts:
        result["best_format"] = "10x_h5"
    elif "loom" in fmts:
        result["best_format"] = "loom"
    elif "10x_mex" in fmts or "archive_10x" in fmts:
        result["best_format"] = "10x_mex_or_archive"
    elif "delimited_matrix" in fmts:
        result["best_format"] = "delimited_matrix"
    elif "rds" in fmts:
        result["best_format"] = "rds"
    elif "h5_generic" in fmts:
        result["best_format"] = "h5_generic"
    elif "matlab" in fmts:
        result["best_format"] = "matlab_matrix"
    elif "archive" in fmts:
        result["best_format"] = "archive_unknown"
    elif "fastq/sra" in fmts:
        result["best_format"] = "fastq_only"
    else:
        result["best_format"] = "none"

    # Quick status
    if result["best_format"] in ("h5ad", "10x_h5", "loom", "10x_mex_or_archive"):
        result["status"] = "high_confidence"
    elif result["best_format"] in ("delimited_matrix", "rds", "h5_generic", "matlab_matrix"):
        result["status"] = "medium_confidence"
    elif result["best_format"] == "archive_unknown":
        result["status"] = "needs_extraction"
    elif result["best_format"] == "fastq_only":
        result["status"] = "fastq_only"
    elif result["best_format"] == "none" and result["supplementary_file_count"] == 0:
        result["status"] = "no_supplementary"
    else:
        result["status"] = "unknown"

    return result


async def process_one(gse_id: str, attempt_download: bool = True) -> dict:
    """Full pipeline for one GSE ID — download, probe, convert, validate."""
    result = {
        "gse_id": gse_id,
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "status": "unknown",
        "strategy": "",
        "best_format": "",
        "files_downloaded": 0,
        "groups_found": 0,
        "conversion_success": False,
        "output_path": "",
        "n_cells": None,
        "n_genes": None,
        "error": "",
        "warnings": [],
    }

    orch = Orchestrator(
        output_dir=OUTPUT_DIR,
        data_dir=DATA_DIR,
    )

    try:
        report = await orch.process(
            gse_id=gse_id,
            auto_download=attempt_download,
            skip_raw_fastq=True,
        )
        result["strategy"] = report.strategy
        result["status"] = report.status.value
        result["warnings"] = report.warnings

        # Count downloaded files
        data_path = DATA_DIR / gse_id
        if data_path.exists():
            result["files_downloaded"] = sum(1 for _ in data_path.rglob("*") if _.is_file())

        # Check if any step succeeded
        for step in report.steps:
            if step.status.value == "completed":
                result["conversion_success"] = True
                if step.output_artifacts:
                    result["output_path"] = str(step.output_artifacts[0])
                break

        # Check final h5ad if it exists
        final_h5ad = OUTPUT_DIR / "h5ad" / f"{gse_id}_final.h5ad"
        if final_h5ad.exists():
            try:
                import anndata
                adata = anndata.read_h5ad(final_h5ad)
                result["n_cells"] = adata.n_obs
                result["n_genes"] = adata.n_vars
                result["has_raw_counts"] = "counts" in adata.layers
            except Exception:
                pass

        if report.steps:
            # Determine best format from groups
            for step in report.steps:
                if step.method and "fastq" not in step.method.lower():
                    result["best_format"] = step.method
                    break

    except Exception as e:
        result["status"] = "error"
        result["error"] = f"{type(e).__name__}: {e}"
        log_error(gse_id, "process", traceback.format_exc())

    return result


async def main():
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
    print(f"=== Batch GSE Test: {len(pending)} remaining of {len(GSE_LIST)} total ===\n")
    print(f"{'#':>4s}  {'GSE':>11s}  {'Status':>14s}  {'Best Format':>22s}  {'DL':>4s}  Details")
    print("-" * 100)

    for i, gse_id in enumerate(pending, 1):
        print(f"\n[{i:>3d}/{len(pending)}] {gse_id} — probing metadata + file listing...")
        result = await probe_only(gse_id)
        log_result(result)
        cleanup(gse_id)

    print(f"\n=== Done. Results written to {RESULTS_FILE} ===")


if __name__ == "__main__":
    asyncio.run(main())
