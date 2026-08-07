#!/usr/bin/env python3
"""Batch evaluation: probe 50 GEO datasets, predict routing, estimate success.

Usage:
    python scripts/batch_evaluate.py [--process N] [--download]

Phases:
    1. Metadata probe (GEO API) — fast, all 50
    2. Supplementary listing (FTP) — fast, all 50
    3. Route prediction — classify each dataset's available formats
    4. Optional: process top-N most diverse datasets end-to-end
"""

from __future__ import annotations

import asyncio
import csv
import json
import os
import sys
import time
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path

# Add project root to path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import httpx


@dataclass
class DatasetAssessment:
    gse_id: str
    # Metadata
    organism: str = ""
    platform: str = ""
    technology: str = ""
    title: str = ""
    n_samples: int = 0
    metadata_error: str | None = None

    # Supplementary files
    n_supplementary: int = 0
    total_size_mb: float = 0.0
    file_types: set[str] = field(default_factory=set)
    has_h5ad: bool = False
    has_10x_h5: bool = False
    has_loom: bool = False
    has_mex: bool = False  # matrix.mtx + barcodes + features
    has_csv_matrix: bool = False
    has_rds: bool = False
    has_h5seurat: bool = False
    has_normalized: bool = False
    has_fastq: bool = False
    has_archive: bool = False
    has_archive_10x: bool = False
    supplementary_error: str | None = None

    # Route prediction
    best_priority: int = 99  # 1-8 from routing rules, 99 = no route
    best_route: str = ""
    predicted_outcome: str = "unknown"  # "clean_success", "partial", "hard_failure", "no_data"
    confidence: float = 0.0
    notes: list[str] = field(default_factory=list)

    # Processing result (filled after actual processing)
    actual_outcome: str = ""
    actual_error: str = ""
    processing_time_s: float = 0.0


async def probe_metadata(gse_id: str, client: httpx.AsyncClient) -> dict:
    """Fetch GEO metadata via the MCP tool's underlying API."""
    # Use GEOparse-style API: esummary + esearch
    import xml.etree.ElementTree as ET

    result = {
        "organism": "",
        "platform": "",
        "technology": "",
        "title": "",
        "n_samples": 0,
        "error": None,
    }

    try:
        # E-summary for the GSE record
        url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
        params = {
            "db": "gds",
            "id": gse_id,
            "retmode": "json",
        }
        resp = await client.get(url, params=params, timeout=30)
        resp.raise_for_status()
        data = resp.json()

        if "result" not in data or gse_id not in data.get("result", {}):
            result["error"] = f"GSE {gse_id} not found in GEO"
            return result

        record = data["result"][gse_id]
        result["title"] = record.get("title", "")
        result["organism"] = record.get("taxon", "")
        result["platform"] = record.get("gpl", "")
        result["n_samples"] = int(record.get("samples", 0))

        # Detect technology from title and summary
        summary = record.get("summary", "")
        combined = f"{result['title']} {summary}".lower()
        if "scrna" in combined or "single cell" in combined or "single-cell" in combined:
            result["technology"] = "scRNA-seq"
        elif "snrna" in combined or "single nucleus" in combined or "single-nucleus" in combined:
            result["technology"] = "snRNA-seq"
        elif "spatial" in combined:
            result["technology"] = "spatial transcriptomics"
        elif "atac" in combined:
            result["technology"] = "scATAC-seq"
        else:
            result["technology"] = "unknown"

    except Exception as exc:
        result["error"] = str(exc)

    return result


async def list_supplementary(gse_id: str, client: httpx.AsyncClient) -> dict:
    """List GEO supplementary files via FTP directory scrape."""
    result = {
        "files": [],
        "n_supplementary": 0,
        "total_size_mb": 0.0,
        "file_types": set(),
        "error": None,
    }

    try:
        # GEO FTP directory listing
        url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{gse_id[:-3]}nnn/{gse_id}/suppl/"
        resp = await client.get(url, timeout=30)

        if resp.status_code == 404:
            # Try without /suppl/
            url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{gse_id[:-3]}nnn/{gse_id}/"
            resp = await client.get(url, timeout=30)

        if resp.status_code != 200:
            result["error"] = f"FTP listing returned {resp.status_code}"
            return result

        text = resp.text

        # Parse Apache/nginx directory listing
        import re
        # Pattern: <a href="filename">filename</a> followed by size digits
        # Apache format: <tr><td><a href="...">...</a></td><td>...</td><td>size</td></tr>
        # Simpler: find all hrefs and then look for sizes nearby
        file_entries = re.findall(
            r'<a\s+href="([^"]+)"[^>]*>([^<]+)</a>[^<]*(?:</td><td[^>]*>)?\s*(\d+(?:\.\d+)?[KMGT]?)\s*(?:</td>)?',
            text,
            re.IGNORECASE,
        )

        seen = set()
        for href, name, size_str in file_entries:
            if href.startswith("/") or href.startswith("?") or href.endswith("/"):
                continue
            if name in seen:
                continue
            seen.add(name)

            # Parse size
            size_bytes = 0
            size_str = size_str.strip().upper()
            if size_str.endswith("G"):
                size_bytes = int(float(size_str[:-1]) * 1_000_000_000)
            elif size_str.endswith("M"):
                size_bytes = int(float(size_str[:-1]) * 1_000_000)
            elif size_str.endswith("K"):
                size_bytes = int(float(size_str[:-1]) * 1_000)
            elif size_str.isdigit():
                size_bytes = int(size_str)

            result["files"].append({
                "name": name,
                "size_bytes": size_bytes,
            })
            result["total_size_mb"] += size_bytes / 1_000_000

            # Classify file type
            name_lower = name.lower()
            if name_lower.endswith(".h5ad") or name_lower.endswith(".h5ad.gz"):
                result["file_types"].add("h5ad")
            elif name_lower.endswith(".h5") and not name_lower.endswith(".h5ad"):
                result["file_types"].add("h5")
            elif name_lower.endswith(".loom") or name_lower.endswith(".loom.gz"):
                result["file_types"].add("loom")
            elif name_lower.endswith(".mtx") or name_lower.endswith(".mtx.gz"):
                result["file_types"].add("mtx")
            elif "barcodes" in name_lower or "features" in name_lower or "genes.tsv" in name_lower:
                result["file_types"].add("mex_support")
            elif name_lower.endswith(".csv") or name_lower.endswith(".csv.gz") or name_lower.endswith(".tsv") or name_lower.endswith(".tsv.gz") or name_lower.endswith(".txt") or name_lower.endswith(".txt.gz"):
                result["file_types"].add("csv")
            elif name_lower.endswith(".rds") or name_lower.endswith(".rds.gz"):
                result["file_types"].add("rds")
            elif name_lower.endswith(".h5seurat"):
                result["file_types"].add("h5seurat")
            elif name_lower.endswith(".fastq") or name_lower.endswith(".fastq.gz") or name_lower.endswith(".fq") or name_lower.endswith(".fq.gz"):
                result["file_types"].add("fastq")
            elif name_lower.endswith(".bam") or name_lower.endswith(".bam.gz"):
                result["file_types"].add("bam")
            elif name_lower.endswith(".tar") or name_lower.endswith(".tar.gz") or name_lower.endswith(".tgz") or name_lower.endswith(".zip"):
                result["file_types"].add("archive")
            elif "raw" in name_lower:
                result["file_types"].add("archive")

        result["n_supplementary"] = len(result["files"])

    except Exception as exc:
        result["error"] = str(exc)

    return result


def predict_route(assessment: DatasetAssessment) -> None:
    """Apply 8-level routing rules to predict best processing route."""
    ft = assessment.file_types

    # Detect MEX triplets (need mtx + barcodes/features support)
    if "mtx" in ft and "mex_support" in ft:
        assessment.has_mex = True
    # Sometimes mtx alone + csv files that could be barcodes/features
    if "mtx" in ft and ("csv" in ft or len(assessment.file_types) >= 3):
        assessment.has_mex = True

    assessment.has_h5ad = "h5ad" in ft
    assessment.has_10x_h5 = "h5" in ft
    assessment.has_loom = "loom" in ft
    assessment.has_csv_matrix = "csv" in ft
    assessment.has_rds = "rds" in ft
    assessment.has_h5seurat = "h5seurat" in ft
    assessment.has_fastq = "fastq" in ft
    assessment.has_archive = "archive" in ft

    # Detect archive_10x pattern: tar.gz files containing MEX triplets
    if assessment.has_archive and not assessment.has_mex:
        assessment.has_archive_10x = True

    # Apply priority rules
    if assessment.has_h5ad:
        assessment.best_priority = 1
        assessment.best_route = "READ_H5AD (validate/standardize)"
        assessment.confidence = 0.95
        assessment.predicted_outcome = "clean_success"
    elif assessment.has_10x_h5:
        assessment.best_priority = 2
        assessment.best_route = "READ_10X_H5"
        assessment.confidence = 0.90
        assessment.predicted_outcome = "clean_success"
    elif assessment.has_loom:
        assessment.best_priority = 3
        assessment.best_route = "READ_LOOM"
        assessment.confidence = 0.85
        assessment.predicted_outcome = "clean_success"
    elif assessment.has_mex:
        assessment.best_priority = 4
        assessment.best_route = "READ_10X_MTX"
        assessment.confidence = 0.85
        assessment.predicted_outcome = "clean_success"
    elif assessment.has_csv_matrix and not assessment.has_archive:
        assessment.best_priority = 5
        assessment.best_route = "READ_RAW_UMI_MATRIX / READ_CSV_MATRIX"
        assessment.confidence = 0.50
        assessment.predicted_outcome = "partial"
        assessment.notes.append("CSV matrix — may need type coercion, may be normalized")
    elif assessment.has_rds or assessment.has_h5seurat:
        assessment.best_priority = 6
        assessment.best_route = "READ_RDS_SEURAT / READ_RDS_SCE / READ_H5SEURAT"
        assessment.confidence = 0.40
        assessment.predicted_outcome = "partial"
        assessment.notes.append("RDS/Seurat — requires R bridge (WSL2)")
    elif assessment.has_archive_10x:
        assessment.best_priority = 4  # Extract then MEX
        assessment.best_route = "ARCHIVE_10X → extract → READ_10X_MTX"
        assessment.confidence = 0.75
        assessment.predicted_outcome = "clean_success"
        assessment.notes.append("Archived 10x — requires extraction before MEX conversion")
    elif assessment.has_archive:
        assessment.best_priority = 7
        assessment.best_route = "ARCHIVE (unknown contents)"
        assessment.confidence = 0.20
        assessment.predicted_outcome = "partial"
        assessment.notes.append("Archive with unknown contents — needs extraction to assess")
    elif assessment.has_fastq and not assessment.has_mex:
        assessment.best_priority = 8
        assessment.best_route = "FASTQ/SRA (raw processing — disabled by default)"
        assessment.confidence = 0.10
        assessment.predicted_outcome = "hard_failure"
        assessment.notes.append("Only FASTQ available — raw processing disabled")
    elif assessment.n_supplementary == 0:
        assessment.best_priority = 99
        assessment.best_route = "NONE — no supplementary files"
        assessment.confidence = 0.0
        assessment.predicted_outcome = "no_data"
        assessment.notes.append("No supplementary files found on GEO FTP")
    else:
        assessment.best_priority = 99
        assessment.best_route = "UNKNOWN"
        assessment.confidence = 0.0
        assessment.predicted_outcome = "hard_failure"

    # Adjust confidence down for multi-sample complex datasets
    if assessment.n_samples > 4:
        assessment.confidence = max(0.0, assessment.confidence - 0.05)
        assessment.notes.append(f"Multi-sample ({assessment.n_samples}) — may need multi-triplet handling")

    # Adjust for non-standard organisms
    if assessment.organism.lower() not in ("homo sapiens", "mus musculus", ""):
        assessment.confidence = max(0.0, assessment.confidence - 0.05)
        assessment.notes.append(f"Non-standard organism: {assessment.organism}")


async def probe_one(gse_id: str, client: httpx.AsyncClient) -> DatasetAssessment:
    """Probe a single dataset: metadata + supplementary + route prediction."""
    a = DatasetAssessment(gse_id=gse_id)

    # Phase 1: Metadata
    meta = await probe_metadata(gse_id, client)
    a.organism = meta.get("organism", "")
    a.platform = meta.get("platform", "")
    a.technology = meta.get("technology", "")
    a.title = meta.get("title", "")[:100]
    a.n_samples = meta.get("n_samples", 0)
    a.metadata_error = meta.get("error")

    # Phase 2: Supplementary
    supp = await list_supplementary(gse_id, client)
    a.n_supplementary = supp["n_supplementary"]
    a.total_size_mb = supp["total_size_mb"]
    a.file_types = supp["file_types"]
    a.supplementary_error = supp.get("error")

    # Phase 3: Predict route
    predict_route(a)

    return a


async def probe_all(gse_ids: list[str], concurrency: int = 5) -> list[DatasetAssessment]:
    """Probe all datasets with controlled concurrency."""
    sem = asyncio.Semaphore(concurrency)

    async def probe_with_limit(gse_id: str) -> DatasetAssessment:
        async with sem:
            print(f"  Probing {gse_id}...")
            try:
                return await probe_one(gse_id, client)
            except Exception as exc:
                a = DatasetAssessment(gse_id=gse_id)
                a.metadata_error = str(exc)
                return a

    async with httpx.AsyncClient(timeout=30, follow_redirects=True) as client:
        tasks = [probe_with_limit(gid) for gid in gse_ids]
        results = await asyncio.gather(*tasks)

    return list(results)


def print_summary(assessments: list[DatasetAssessment]) -> None:
    """Print a summary table of results."""
    clean = [a for a in assessments if a.predicted_outcome == "clean_success"]
    partial = [a for a in assessments if a.predicted_outcome == "partial"]
    failed = [a for a in assessments if a.predicted_outcome == "hard_failure"]
    no_data = [a for a in assessments if a.predicted_outcome == "no_data"]

    print("\n" + "=" * 100)
    print(f"BATCH EVALUATION SUMMARY — {len(assessments)} datasets")
    print("=" * 100)

    print(f"\n{'Outcome':<20} {'Count':<8} {'%':<8}")
    print("-" * 40)
    print(f"{'Clean success':<20} {len(clean):<8} {len(clean)/len(assessments)*100:.1f}%")
    print(f"{'Partial':<20} {len(partial):<8} {len(partial)/len(assessments)*100:.1f}%")
    print(f"{'Hard failure':<20} {len(failed):<8} {len(failed)/len(assessments)*100:.1f}%")
    print(f"{'No data':<20} {len(no_data):<8} {len(no_data)/len(assessments)*100:.1f}%")

    print(f"\n--- Route Distribution ---")
    route_counts: dict[str, int] = {}
    for a in assessments:
        route_counts[a.best_route] = route_counts.get(a.best_route, 0) + 1
    for route, count in sorted(route_counts.items(), key=lambda x: -x[1]):
        print(f"  {route:<50} {count:>3}")

    print(f"\n--- Organism Distribution ---")
    org_counts: dict[str, int] = {}
    for a in assessments:
        org = a.organism or "unknown"
        org_counts[org] = org_counts.get(org, 0) + 1
    for org, count in sorted(org_counts.items(), key=lambda x: -x[1]):
        print(f"  {org:<40} {count:>3}")

    print(f"\n--- File Format Availability ---")
    format_keys = ["has_h5ad", "has_10x_h5", "has_loom", "has_mex", "has_csv_matrix",
                   "has_rds", "has_h5seurat", "has_archive_10x", "has_archive", "has_fastq"]
    for key in format_keys:
        count = sum(1 for a in assessments if getattr(a, key))
        if count > 0:
            print(f"  {key[4:]:<25} {count:>3} datasets ({count/len(assessments)*100:.0f}%)")

    # Print details for failures and partials
    if failed:
        print(f"\n--- HARD FAILURES ({len(failed)}) ---")
        for a in failed:
            print(f"  {a.gse_id}: {a.best_route}")
            for note in a.notes:
                print(f"    └─ {note}")

    if no_data:
        print(f"\n--- NO DATA ({len(no_data)}) ---")
        for a in no_data:
            print(f"  {a.gse_id}: {a.supplementary_error or 'No files found'}")

    print(f"\n--- DETAILED TABLE ---")
    print(f"{'GSE':<12} {'Organism':<25} {'Samples':<8} {'Files':<6} {'Size(MB)':<10} {'Best Route':<45} {'Outcome':<15} {'Conf':<6}")
    print("-" * 130)
    for a in sorted(assessments, key=lambda x: (x.predicted_outcome, -x.confidence)):
        org_short = a.organism[:24] if a.organism else "?"
        print(f"{a.gse_id:<12} {org_short:<25} {a.n_samples:<8} {a.n_supplementary:<6} "
              f"{a.total_size_mb:<10.1f} {a.best_route[:44]:<45} {a.predicted_outcome:<15} {a.confidence:.2f}")


def main():
    import argparse
    parser = argparse.ArgumentParser(description="Batch evaluate GEO datasets")
    parser.add_argument("--dataset-file", default="GSE_single_cell_RNAseq_50.txt",
                        help="File with GSE IDs (one per line)")
    parser.add_argument("--limit", type=int, default=0,
                        help="Limit to first N datasets (0 = all)")
    parser.add_argument("--output", default="data/batch_evaluation.json",
                        help="Output JSON path")
    parser.add_argument("--concurrency", type=int, default=5,
                        help="Concurrent probes")
    args = parser.parse_args()

    # Read GSE IDs
    dataset_file = Path(args.dataset_file)
    if not dataset_file.exists():
        dataset_file = Path(__file__).resolve().parent.parent / args.dataset_file

    gse_ids = []
    with open(dataset_file) as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith("#"):
                gse_ids.append(line)

    if args.limit:
        gse_ids = gse_ids[:args.limit]

    print(f"Loaded {len(gse_ids)} GSE IDs from {dataset_file}")
    print(f"Concurrency: {args.concurrency}")
    print()

    # Run probes
    t0 = time.monotonic()
    assessments = asyncio.run(probe_all(gse_ids, concurrency=args.concurrency))
    elapsed = time.monotonic() - t0

    print(f"\nProbing completed in {elapsed:.1f}s ({elapsed/len(gse_ids):.1f}s per dataset)")

    # Print summary
    print_summary(assessments)

    # Save to JSON
    output_path = Path(args.output)
    if not output_path.is_absolute():
        output_path = Path(__file__).resolve().parent.parent / args.output
    output_path.parent.mkdir(parents=True, exist_ok=True)

    data = []
    for a in assessments:
        data.append({
            "gse_id": a.gse_id,
            "organism": a.organism,
            "platform": a.platform,
            "technology": a.technology,
            "n_samples": a.n_samples,
            "n_supplementary": a.n_supplementary,
            "total_size_mb": round(a.total_size_mb, 1),
            "file_types": sorted(a.file_types),
            "best_priority": a.best_priority,
            "best_route": a.best_route,
            "predicted_outcome": a.predicted_outcome,
            "confidence": round(a.confidence, 2),
            "notes": a.notes,
            "metadata_error": a.metadata_error,
            "supplementary_error": a.supplementary_error,
        })

    with open(output_path, "w") as f:
        json.dump({
            "evaluated_at": datetime.now(timezone.utc).isoformat(),
            "n_datasets": len(assessments),
            "probe_time_s": round(elapsed, 1),
            "datasets": data,
        }, f, indent=2)

    print(f"\nResults saved to: {output_path}")


if __name__ == "__main__":
    main()
