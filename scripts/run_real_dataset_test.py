#!/usr/bin/env python3
"""Real dataset test — runs the full pipeline on a GSE dataset from GSE_single_cell_RNAseq_50.txt.

Usage:
    python scripts/run_real_dataset_test.py [GSE_ID]

If no GSE_ID is provided, uses the first dataset in the list that has local data.
"""

from __future__ import annotations

import asyncio
import sys
import time
from pathlib import Path

# Add project root to path
PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from gse_router.pipeline.orchestrator import Orchestrator


async def main(gse_id: str, auto_download: bool = False) -> None:
    print(f"{'='*70}")
    print(f"  GSE Router — Real Dataset Test")
    print(f"  Dataset: {gse_id}")
    print(f"  Auto-download: {auto_download}")
    print(f"{'='*70}\n")

    data_dir = PROJECT_ROOT / "data"
    output_dir = PROJECT_ROOT / "output"

    orchestrator = Orchestrator(
        output_dir=output_dir,
        data_dir=data_dir,
        download_backend="auto",
    )

    t0 = time.monotonic()

    report = await orchestrator.process(
        gse_id=gse_id,
        skip_metadata_probe=False,
        skip_raw_fastq=True,
        auto_download=auto_download,
    )

    elapsed = time.monotonic() - t0

    # --- Report ---
    print(f"\n{'='*70}")
    print(f"  RESULTS")
    print(f"{'='*70}")
    print(f"  Status:       {report.status.value}")
    print(f"  Strategy:     {report.strategy}")
    print(f"  Reference:    {report.reference_genome}")
    print(f"  Duration:     {elapsed:.1f}s")
    if report.output:
        print(f"  Output:       {report.output.file_path} ({report.output.file_size_mb:.1f} MB)")
        print(f"  Cells:        {report.output.n_obs}")
        print(f"  Genes:        {report.output.n_vars}")
    if report.warnings:
        print(f"  Warnings:     {len(report.warnings)}")
        for w in report.warnings[:5]:
            print(f"    WARNING: {w}")
        if len(report.warnings) > 5:
            print(f"    ... and {len(report.warnings) - 5} more")

    print(f"\n  Steps ({len(report.steps)}):")
    for step in report.steps:
        icon = {"completed": "[OK]", "failed": "[FAIL]", "skipped": "[SKIP]", "pending": "[PEND]"}.get(
            step.status.value, "[???]"
        )
        print(f"    {icon} [{step.status.value.upper():11s}] {step.method}")
        if step.error_message:
            print(f"         Error: {step.error_message[:120]}")
        if step.warnings:
            for w in step.warnings[:2]:
                print(f"         WARNING: {w[:120]}")

    if report.qc_metrics:
        print(f"\n  QC Metrics ({len(report.qc_metrics)}):")
        for m in report.qc_metrics[:10]:
            status = "PASS" if m.passed else "FAIL"
            print(f"    {status} {m.name}: {m.value}")

    # Summary
    n_completed = sum(1 for s in report.steps if s.status.value == "completed")
    n_failed = sum(1 for s in report.steps if s.status.value == "failed")
    n_skipped = sum(1 for s in report.steps if s.status.value == "skipped")
    print(f"\n  Summary: {n_completed} completed, {n_failed} failed, {n_skipped} skipped")
    print(f"{'='*70}")

    # Write report to file
    report_path = output_dir / f"{gse_id}_report.json"
    report_path.write_text(report.model_dump_json(indent=2))
    print(f"\n  Full report: {report_path}")


if __name__ == "__main__":
    # Pick dataset
    if len(sys.argv) > 1:
        gse_id = sys.argv[1].strip().upper()
    else:
        # Read first ID from the list, preferring one with existing data
        ids_file = PROJECT_ROOT / "GSE_single_cell_RNAseq_50.txt"
        ids = [line.strip() for line in ids_file.read_text().splitlines() if line.strip()]
        # Prefer datasets with existing data
        existing = [gid for gid in ids if (PROJECT_ROOT / "data" / gid).exists()]
        gse_id = existing[0] if existing else ids[0]

    auto_download = "--download" in sys.argv

    asyncio.run(main(gse_id, auto_download=auto_download))
