"""Process GEO datasets one at a time — download, convert, validate, delete data.
Keeps results in batch_convert_results.jsonl."""
from __future__ import annotations

import asyncio
import json
import shutil
import sys
import traceback
from datetime import datetime, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from gse_router.pipeline.orchestrator import Orchestrator

DATA_DIR = Path("E:/project/Aging_OS2/data")
OUTPUT_DIR = Path("E:/project/Aging_OS2/output")
RESULTS_FILE = Path("E:/project/Aging_OS2/batch_convert_results.jsonl")

# Priority order: small + high-confidence first
BATCH = [
    # (gse_id, description, expected_format)
    ("GSE336973", "raw_counts_matrix_txt.gz", "delimited_raw"),
    ("GSE319592", "TPM_normalized_csv.gz", "delimited_normalized"),
    ("GSE337384", "salmon_counts_tpm", "delimited_mixed"),
    ("GSE318657", "small_RAW.tar_10x_v3", "archive_10x"),
    ("GSE328525", "small_RAW.tar_10x_v3", "archive_10x"),
    ("GSE328097", "small_RAW.tar", "archive_10x"),
    ("GSE248385", "small_RAW.tar", "archive_10x"),
    ("GSE327186", "small_RAW.tar", "archive_10x"),
    ("GSE276559", "medium_RAW.tar_10x_v3", "archive_10x"),
    ("GSE337368", "medium_RAW.tar", "archive_10x"),
    ("GSE309120", "10x_MEX_exposed", "10x_mex"),
    ("GSE327907", "seurat_v5_rds", "rds"),
    ("GSE314267", "seurat_rds", "rds"),
    ("GSE337162", "seurat_rds_large", "rds"),
    ("GSE301594", "RAW.tar_human", "archive_10x"),
    ("GSE301966", "RAW.tar_mouse", "archive_10x"),
    ("GSE335941", "RAW.tar_10x_v3", "archive_10x"),
    ("GSE335811", "RAW.tar_10x_v3", "archive_10x"),
    ("GSE327732", "RAW.tar_10x_v3", "archive_10x"),
    ("GSE317299", "RAW.tar_10x_v3_single_sample", "archive_10x"),
    ("GSE271348", "RAW.tar_10x_v3_large", "archive_10x"),
    ("GSE322527", "RAW.tar_10x_v2", "archive_10x"),
    ("GSE303713", "RAW.tar_novaseq", "archive_10x"),
    ("GSE311869", "RAW.tar", "archive_10x"),
    ("GSE308872", "RAW.tar", "archive_10x"),
    ("GSE332958", "RAW.tar", "archive_10x"),
    ("GSE329924", "RAW.tar", "archive_10x"),
    ("GSE336668", "RAW.tar", "archive_10x"),
    ("GSE336400", "RAW.tar", "archive_10x"),
    ("GSE314969", "RAW.tar", "archive_10x"),
    # Larger ones — still try
    ("GSE301717", "RAW.tar_10x_medium", "archive_10x"),
    ("GSE248386", "RAW.tar_large", "archive_10x"),
    ("GSE319066", "RAW.tar_large", "archive_10x"),
    ("GSE318800", "RAW.tar_large", "archive_10x"),
    ("GSE300504", "RAW.tar_zebrafish", "archive_10x"),
    ("GSE327808", "RAW.tar_rat", "archive_10x"),
    ("GSE329941", "RAW.tar_large", "archive_10x"),
    ("GSE337230", "RAW.tar_large", "archive_10x"),
    ("GSE337208", "RAW.tar_large", "archive_10x"),
    ("GSE337289", "RAW.tar_chicken", "archive_10x"),
    ("GSE327389", "RAW.tar_mouse_single", "archive_10x"),
    ("GSE325198", "rds_large", "rds"),
    ("GSE320454", "RAW.tar_plus_RDS", "mixed"),
    ("GSE336843", "RAW.tar_large", "archive_10x"),
    ("GSE299068", "RAW.tar_large", "archive_10x"),
    ("GSE326370", "RAW.tar_pig_large", "archive_10x"),
    ("GSE319693", "RAW.tar_huge", "archive_10x"),
    ("GSE337655", "RAW.tar_fly", "archive_10x"),
    ("GSE336861", "RAW.tar_large", "archive_10x"),
    ("GSE322748", "RAW.tar_very_large", "archive_10x"),
]


def cleanup(gse_id: str):
    """Remove all downloaded/processed data for this GSE."""
    for p in [
        DATA_DIR / gse_id,
        OUTPUT_DIR / "h5ad",
    ]:
        if p.exists():
            try:
                shutil.rmtree(p)
            except Exception:
                pass
    # Clean extracted dirs
    if DATA_DIR.exists():
        for d in DATA_DIR.iterdir():
            if d.is_dir() and d.name.startswith("_extracted_"):
                try:
                    shutil.rmtree(d)
                except Exception:
                    pass


async def process_one(gse_id: str) -> dict:
    result = {
        "gse_id": gse_id,
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "status": "unknown",
        "strategy": "",
        "route_method": "",
        "groups_found": 0,
        "group_types": [],
        "conversion_success": False,
        "output_n_cells": None,
        "output_n_genes": None,
        "has_raw_counts": False,
        "error": "",
        "warnings": [],
    }

    orch = Orchestrator(output_dir=OUTPUT_DIR, data_dir=DATA_DIR)

    try:
        report = await orch.process(
            gse_id=gse_id,
            auto_download=True,
            skip_raw_fastq=True,
        )

        result["strategy"] = report.strategy
        result["groups_found"] = len(report.steps)
        result["warnings"] = report.warnings

        for step in report.steps:
            if step.method:
                result["route_method"] = step.method
                break

        for step in report.steps:
            if step.status.value == "completed":
                result["conversion_success"] = True
                result["group_types"].append(step.method)
            elif step.status.value == "failed":
                result["error"] = step.error_message or "step failed"

        overall = report.status.value
        if result["conversion_success"]:
            result["status"] = "success"
        elif "limited" in str(result.get("warnings", "")).lower():
            result["status"] = "partial"
        elif overall == "failed":
            result["status"] = "failed"
        else:
            result["status"] = overall

        # Check final output
        final_h5ad = OUTPUT_DIR / "h5ad" / f"{gse_id}_final.h5ad"
        if final_h5ad.exists():
            try:
                import anndata
                adata = anndata.read_h5ad(final_h5ad)
                result["output_n_cells"] = adata.n_obs
                result["output_n_genes"] = adata.n_vars
                result["has_raw_counts"] = "counts" in adata.layers
            except Exception as e:
                result["error"] = (result["error"] + "; " if result["error"] else "") + f"Final read failed: {e}"

    except Exception as e:
        result["status"] = "exception"
        result["error"] = f"{type(e).__name__}: {e}"

    return result


async def main():
    already = set()
    if RESULTS_FILE.exists():
        with open(RESULTS_FILE) as f:
            for line in f:
                try:
                    already.add(json.loads(line.strip()).get("gse_id", ""))
                except Exception:
                    pass

    pending = [(g, d, e) for g, d, e in BATCH if g not in already]
    print(f"=== Batch Convert: {len(pending)} remaining ===\n")

    success = partial = failed = 0
    for i, (gse_id, desc, expected) in enumerate(pending, 1):
        n = len(pending)
        print(f"\n{'='*80}")
        print(f"[{i:>2d}/{n}] {gse_id} ({desc}) — expected: {expected}")
        print(f"{'='*80}")

        result = await process_one(gse_id)

        # Write result immediately
        with open(RESULTS_FILE, "a") as f:
            f.write(json.dumps(result) + "\n")

        # Summary line
        cells = result.get("output_n_cells") or 0
        genes = result.get("output_n_genes") or 0
        raw = "raw" if result.get("has_raw_counts") else "no_raw"
        err = result.get("error", "")[:100]
        print(f"  -> {result['status']:>10s} | strategy={result['strategy']:>25s} | "
              f"cells={cells:>6d} genes={genes:>6d} {raw} | {err}")

        if result["status"] == "success":
            success += 1
        elif result["status"] == "partial":
            partial += 1
        else:
            failed += 1

        print(f"  Running tally: {success} success | {partial} partial | {failed} failed")

        # Cleanup
        cleanup(gse_id)

    print(f"\n=== Final: {success} success | {partial} partial | {failed} failed ===")


if __name__ == "__main__":
    asyncio.run(main())
