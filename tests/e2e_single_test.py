"""End-to-end test: single GSE dataset -> .h5ad file using the real Orchestrator.

Usage: python -m tests.e2e_single_test [GSE_ID]
Default: GSE327186 (mouse, 38MB RAW.tar with 10x MEX triplet)
"""
from __future__ import annotations

import asyncio
import shutil
import sys
import time
from pathlib import Path

# Load .env before anything else
from gse_router.config import get_config
get_config()

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "src"))

DATA_DIR = PROJECT_ROOT / "data"
OUTPUT_DIR = PROJECT_ROOT / "output"


async def run_e2e(gse_id: str):
    """Full pipeline using the real Orchestrator."""
    from gse_router.pipeline.orchestrator import Orchestrator

    t0 = time.time()
    print(f"{'='*70}")
    print(f"E2E Test: {gse_id}")
    print(f"{'='*70}")

    # Clean any stale data
    dataset_dir = DATA_DIR / gse_id
    if dataset_dir.exists():
        shutil.rmtree(dataset_dir)

    orch = Orchestrator(output_dir=OUTPUT_DIR, data_dir=DATA_DIR)

    print("\nRunning full pipeline (download -> probe -> plan -> execute -> validate)...")
    report = await orch.process(
        gse_id=gse_id,
        data_path=None,
        skip_metadata_probe=False,
        dry_run=False,
        auto_download=True,
    )

    elapsed = time.time() - t0

    print(f"\n{'='*70}")
    print(f"RESULTS: {gse_id}")
    print(f"  Status: {report.status.value}")
    print(f"  Strategy: {report.strategy}")
    print(f"  Steps: {len(report.steps)}")
    print(f"  Total time: {elapsed:.1f}s")

    for step in report.steps:
        status_icon = "OK" if step.status.value == "completed" else "FAIL"
        print(f"    [{status_icon}] {step.method:30s} -> {step.output_artifacts}")
        if step.error_message:
            print(f"          error: {step.error_message}")

    if report.output and report.output.file_path:
        h5ad_path = Path(report.output.file_path)
        if h5ad_path.exists():
            print(f"\n  Output h5ad: {h5ad_path}")
            print(f"  Size: {h5ad_path.stat().st_size/1e6:.1f} MB")
            print(f"  Cells: {report.output.n_obs}")
            print(f"  Genes: {report.output.n_vars}")
            print(f"  Has raw counts: {report.output.has_raw_counts}")
            if report.qc_metrics:
                for m in report.qc_metrics:
                    print(f"  QC: {m.name} = {m.value}")
            print(f"\n[PASS] E2E TEST PASSED - .h5ad file produced successfully!")
        else:
            print(f"\n[FAIL] Output path recorded but file missing: {report.output.file_path}")
    else:
        print(f"\n[FAIL] E2E TEST FAILED - No .h5ad file produced")
        if report.errors:
            for err in report.errors:
                print(f"  Error: {err}")
        for step in report.steps:
            if step.error_message:
                print(f"  Step error [{step.step_id}]: {step.error_message}")

    # Cleanup
    if dataset_dir.exists():
        shutil.rmtree(dataset_dir)
        print(f"\nCleaned up {dataset_dir}")

    return report


if __name__ == "__main__":
    gse_id = sys.argv[1] if len(sys.argv) > 1 else "GSE327186"
    asyncio.run(run_e2e(gse_id))
