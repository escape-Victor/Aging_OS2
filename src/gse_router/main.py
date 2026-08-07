"""CLI entry point for GSE Router.

Provides commands for probing, planning, and processing GEO datasets.

Configuration:
    Set NCBI_API_KEY and NCBI_EMAIL in a .env file or as environment variables.
    See .env.example for all available settings.
"""

from __future__ import annotations

import asyncio
import json
import sys
from pathlib import Path

import click

from .config import get_config
from .pipeline.orchestrator import Orchestrator


@click.group()
@click.version_option(version="0.1.0", prog_name="gse-router")
def cli() -> None:
    """GSE Router — Smart single-cell RNA-seq routing proxy.

    Converts arbitrary GEO single-cell datasets into standardized AnnData objects
    by probing files, reasoning about biology, and executing the right tools.

    \b
    Quick start:
      gse-router process GSE123456 --auto-download
      gse-router probe GSE123456
      gse-router list-supplementary GSE123456
    """
    pass


@cli.command()
@click.argument("gse_id")
@click.option("--data-dir", "-d", type=click.Path(exists=True, path_type=Path),
              help="Path to downloaded dataset files.")
@click.option("--output-dir", "-o", type=click.Path(path_type=Path),
              default=None, help="Output directory (default from config/env).")
@click.option("--skip-metadata", is_flag=True, help="Skip GEO/SRA metadata probe.")
@click.option("--dry-run", is_flag=True, help="Plan only, do not execute processing.")
@click.option("--email", "-e", help="Email for NCBI E-utilities (reads NCBI_EMAIL env var if unset).")
@click.option("--auto-download/--no-auto-download", default=None,
              help="Auto-download GEO supplementary files before processing (default from config).")
@click.option("--api-key", help="NCBI API key (reads NCBI_API_KEY env var if unset).")
def process(
    gse_id: str,
    data_dir: Path | None,
    output_dir: Path | None,
    skip_metadata: bool,
    dry_run: bool,
    email: str | None,
    auto_download: bool | None,
    api_key: str | None,
) -> None:
    """Process a GEO dataset end-to-end.

    GSE_ID is the GEO Series accession (e.g., GSE123456).

    Runs the full pipeline: (auto-download) -> metadata probe -> file probe ->
    manifest -> route plan -> execute -> validate -> report.

    \b
    API keys: set NCBI_API_KEY and NCBI_EMAIL in .env or as environment
    variables. Without an API key, NCBI rate-limits to 3 requests/second.
    Register at: https://www.ncbi.nlm.nih.gov/account/
    """
    cfg = get_config()

    # Resolve settings: CLI flags > env vars > defaults
    resolved_email = email or cfg.ncbi_email or None
    resolved_api_key = api_key or cfg.ncbi_api_key or None
    resolved_output = output_dir or cfg.output_dir
    resolved_auto_download = auto_download if auto_download is not None else cfg.auto_download

    orch = Orchestrator(
        output_dir=resolved_output,
        email=resolved_email,
        api_key=resolved_api_key,
    )
    report = asyncio.run(orch.process(
        gse_id=gse_id,
        data_path=data_dir,
        skip_metadata_probe=skip_metadata,
        dry_run=dry_run,
        auto_download=resolved_auto_download,
    ))
    click.echo(f"Status: {report.status.value}")
    click.echo(f"Steps: {len(report.steps)}")
    click.echo(f"Report: {resolved_output / f'{gse_id}_report.json'}")
    if report.errors:
        click.echo(f"Errors: {len(report.errors)}")
        for err in report.errors:
            click.echo(f"  - {err}")


@cli.command()
@click.argument("gse_id")
@click.option("--data-dir", "-d", type=click.Path(exists=True, path_type=Path),
              help="Path to downloaded dataset files.")
@click.option("--output-dir", "-o", type=click.Path(path_type=Path),
              default=Path.cwd() / "output", help="Output directory.")
def probe(gse_id: str, data_dir: Path | None, output_dir: Path) -> None:
    """Probe a GEO dataset and write sidecars + manifest (no processing)."""
    from .probes.file_probe import FileProbe
    from .router.classifier import FileClassifier

    data_path = data_dir or Path.cwd() / "data" / gse_id
    output_dir.mkdir(parents=True, exist_ok=True)

    click.echo(f"Probing {gse_id}...")

    # File probe
    fp = FileProbe()
    sidecars = fp.probe_directory(data_path)
    sidecar_dict = {sc.file_id: sc for sc in sidecars}

    # Write sidecars
    sidecar_dir = output_dir / "sidecars"
    sidecar_dir.mkdir(exist_ok=True)
    for sc in sidecars:
        sidecar_path = sidecar_dir / f"{sc.file_id.replace(':', '_')}.json"
        sidecar_path.write_text(sc.model_dump_json(indent=2))
    click.echo(f"  File sidecars: {len(sidecars)} -> {sidecar_dir}")

    # Classify groups
    classifier = FileClassifier()
    groups = classifier.group_files(sidecar_dict)
    click.echo(f"  File groups: {len(groups)}")

    # Build and write manifest
    from .schemas.manifest import DatasetManifest
    manifest = DatasetManifest(
        dataset_id=gse_id,
        files=sidecar_dict,
        file_groups=groups,
        selected_strategy=classifier.select_strategy(groups),
    )
    manifest_path = output_dir / f"{gse_id}_manifest.json"
    manifest_path.write_text(manifest.model_dump_json(indent=2))
    click.echo(f"  Manifest: {manifest_path}")
    click.echo(f"  Strategy: {manifest.selected_strategy.value}")


@cli.command()
@click.argument("gse_id")
@click.option("--data-dir", "-d", type=click.Path(exists=True, path_type=Path),
              help="Path to downloaded dataset files.")
@click.option("--manifest", "-m", "manifest_path", type=click.Path(exists=True, path_type=Path),
              help="Path to an existing manifest JSON (skips file probe).")
def plan(gse_id: str, data_dir: Path | None, manifest_path: Path | None) -> None:
    """Plan processing route for a GEO dataset (no execution)."""
    from .router.route_planner import RoutePlanner

    if manifest_path:
        from .schemas.manifest import DatasetManifest
        manifest = DatasetManifest.model_validate_json(manifest_path.read_text())
    else:
        click.echo("Run 'gse-router probe' first, or provide --manifest.")
        sys.exit(1)

    planner = RoutePlanner()
    route_plan = planner.plan(manifest)

    click.echo(json.dumps(json.loads(route_plan.model_dump_json()), indent=2))


@cli.command()
@click.argument("gse_id")
@click.option("--data-dir", type=click.Path(path_type=Path), default=None,
              help="Root data directory (default: ./data).")
@click.option("--file-types", multiple=True, default=None,
              help="File types to download (e.g. h5ad h5 loom mtx).")
@click.option("--max-size-mb", type=int, default=500,
              help="Maximum file size in MB (default: 500).")
def download(
    gse_id: str,
    data_dir: Path | None,
    file_types: tuple[str, ...],
    max_size_mb: int,
) -> None:
    """Download GEO supplementary files for a GSE accession.

    GSE_ID is the GEO Series accession (e.g., GSE123456).

    Downloads supplementary files from GEO FTP. By default, only processed
    data files are downloaded (archives and raw FASTQ/BAM are skipped).
    Use --file-types to restrict to specific file categories.
    """
    from .tools.geo_downloader import GeoFileDownloader

    data_dir = data_dir or Path.cwd() / "data"
    downloader = GeoFileDownloader(data_dir=data_dir)

    file_types_list = list(file_types) if file_types else None
    if file_types_list:
        click.echo(f"Downloading {gse_id} (types: {', '.join(file_types_list)})...")
    else:
        click.echo(f"Downloading {gse_id} (all processed files)...")

    try:
        dataset_dir, all_files, downloaded = asyncio.run(
            downloader.download_dataset(
                gse_id=gse_id,
                file_types=file_types_list,
                max_size_mb=max_size_mb,
            )
        )
    finally:
        asyncio.run(downloader.close())

    if not downloaded:
        click.echo(f"No files downloaded for {gse_id}.")
        if all_files:
            click.echo(f"  ({len(all_files)} files available but none matched filters)")
        return

    total_size = sum(p.stat().st_size for p in downloaded if p.exists())
    click.echo(f"Downloaded {len(downloaded)} files ({total_size / 1e6:.1f} MB) to:")
    click.echo(f"  {dataset_dir}")


@cli.command("list-supplementary")
@click.argument("gse_id")
def list_supplementary_cmd(gse_id: str) -> None:
    """List supplementary files available for a GSE accession.

    GSE_ID is the GEO Series accession (e.g., GSE123456).

    Queries the GEO FTP directory listing and displays available
    supplementary files with their sizes.
    """
    from .tools.geo_downloader import GeoSupplementaryLister

    lister = GeoSupplementaryLister()
    try:
        files = asyncio.run(lister.list_supplementary(gse_id))
    finally:
        asyncio.run(lister.close())

    if not files:
        click.echo(f"No supplementary files found for {gse_id}.")
        return

    click.echo(f"Supplementary files for {gse_id} ({len(files)} total):")
    click.echo("-" * 80)
    for f in files:
        size_str = f"{f.size_bytes / 1e6:.1f} MB" if f.size_bytes else "unknown"
        cat_str = f"[{f.category}]"
        click.echo(f"  {f.filename:60s} {size_str:>10s}  {cat_str:>8s}")
    click.echo("-" * 80)
    total_mb = sum(f.size_bytes or 0 for f in files) / 1e6
    click.echo(f"  Total: {len(files)} files, {total_mb:.1f} MB")


@cli.command()
@click.argument("h5ad_path", type=click.Path(exists=True, path_type=Path))
def validate(h5ad_path: Path) -> None:
    """Validate a single h5ad file against QC standards."""
    from .tools.validators import AnnDataValidator

    validator = AnnDataValidator()
    result = validator.validate(h5ad_path)

    click.echo(f"Passed: {result.passed}")
    click.echo(f"Metrics: {json.dumps(result.metrics, indent=2, default=str)}")
    if result.failures:
        click.echo(f"Failures ({len(result.failures)}):")
        for f in result.failures:
            click.echo(f"  - {f}")
    if result.warnings:
        click.echo(f"Warnings ({len(result.warnings)}):")
        for w in result.warnings:
            click.echo(f"  - {w}")


def main() -> None:
    """Entry point for console_scripts."""
    cli()


if __name__ == "__main__":
    main()
