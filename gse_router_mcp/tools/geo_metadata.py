"""GEO metadata MCP tools -- fetch_geo_metadata, list_geo_supplementary,
read_geo_filelist, download_geo_supplementary.

Each tool is a standalone async function callable from the MCP server.
Also provides a CLI entry point via ``python -m gse_router_mcp.tools.geo_metadata``.

References:
    config/tools.yaml -- tool-to-environment dispatch configuration
    src/gse_router/probes/metadata.py -- MetadataProbe (wrapped by fetch_geo_metadata)
"""

from __future__ import annotations

import asyncio
import hashlib
import json
import os
import re
import sys
from pathlib import Path
from typing import Any

import httpx
import structlog

# Ensure src/ is on sys.path so we can import gse_router
_src = Path(__file__).resolve().parent.parent.parent / "src"
if str(_src) not in sys.path:
    sys.path.insert(0, str(_src))

logger = structlog.get_logger(__name__)

# ---------------------------------------------------------------------------
# GEO FTP URL helpers
# ---------------------------------------------------------------------------

_GEO_FTP_BASE = "https://ftp.ncbi.nlm.nih.gov/geo/series"


def _geo_suppl_url(gse_id: str) -> str:
    """Build the supplementary directory URL for a GSE accession."""
    # GEO organises series into <prefix>nnn/ directories
    prefix = gse_id[:-3]
    return f"{_GEO_FTP_BASE}/{prefix}nnn/{gse_id}/suppl/"


def _geo_file_url(gse_id: str, filename: str) -> str:
    """Build the download URL for a single supplementary file."""
    return f"{_geo_suppl_url(gse_id)}{filename}"


def _geo_filelist_url(gse_id: str) -> str:
    """Build the URL for the filelist.txt of a GSE accession."""
    return f"{_geo_suppl_url(gse_id)}filelist.txt"


# ---------------------------------------------------------------------------
# GeoSupplementaryLister -- list and inspect GEO supplementary files
# ---------------------------------------------------------------------------


class GeoSupplementaryLister:
    """Lists and inspects supplementary files on GEO FTP without downloading.

    Uses the NCBI GEO FTP directory index (HTML) and the per-dataset
    ``filelist.txt`` manifest to discover what files are available.
    """

    def __init__(self, timeout: float = 30.0) -> None:
        self.timeout = timeout
        self._log = logger.bind(component="GeoSupplementaryLister")

    # -- public API ------------------------------------------------------

    async def list_supplementary(self, gse_id: str) -> dict[str, Any]:
        """Fetch the supplementary file listing for a GSE dataset.

        Scrapes the NCBI FTP directory index page.  Because the FTP server
        returns a plain HTML directory listing, we extract ``href``
        attributes from ``<a>`` tags.

        Args:
            gse_id: GEO Series accession (e.g. ``"GSE123456"``).

        Returns:
            A dict with keys ``gse_id``, ``files`` (list of
            ``{filename, size_bytes, url}``), and ``count``.
        """
        url = _geo_suppl_url(gse_id)
        self._log.info("Listing supplementary files", gse_id=gse_id, url=url)

        try:
            async with httpx.AsyncClient(timeout=self.timeout) as client:
                resp = await client.get(url)
        except httpx.TimeoutException:
            return _error(gse_id, f"Timeout listing supplementary files (>{self.timeout}s)")
        except httpx.RequestError as exc:
            return _error(gse_id, f"HTTP request failed: {exc}")

        if resp.status_code != 200:
            return _error(
                gse_id,
                f"GEO FTP returned HTTP {resp.status_code} for supplementary listing",
            )

        files: list[dict[str, Any]] = []
        # The FTP index page lists files in <a> tags with an href and optional
        # size in the following text node.  We use a simple regex to extract
        # both filename and size hint.
        text = resp.text
        # Pattern: <a href="filename">filename</a> ... optional size digits
        # Some FTP listings show size right after the link, e.g.:
        #   <a href="foo.gz">foo.gz</a>  2024-03-15 13:22  1.2M
        # Others just show the bare filename link.
        entries = re.findall(
            r'<a\s+href="([^"]+)"[^>]*>.*?</a>\s*[\d-]+\s+[\d:]+\s+([\d.]+[KMGT]?)?',
            text,
        )
        seen: set[str] = set()
        for href, size_str in entries:
            filename = href.rstrip("/")
            if not filename or filename in ("..", "Parent Directory", "/"):
                continue
            if filename.startswith("?"):
                continue
            if filename in seen:
                continue
            seen.add(filename)

            size_bytes = _parse_ftp_size(size_str.strip()) if size_str and size_str.strip() else None
            files.append({
                "filename": filename,
                "size_bytes": size_bytes,
                "url": _geo_file_url(gse_id, filename),
            })

        self._log.info("Supplementary listing complete", gse_id=gse_id, count=len(files))
        return {"gse_id": gse_id, "files": files, "count": len(files)}

    async def read_filelist(self, gse_id: str) -> dict[str, Any]:
        """Read and parse the ``filelist.txt`` manifest for a GSE dataset.

        GEO's ``filelist.txt`` uses the format ``<md5>  <filename>`` (one
        entry per line).  Reading this before downloading a ``RAW.tar`` is a
        cardinal rule of the router -- many ``RAW.tar`` archives contain
        processed matrices, not FASTQ.

        Args:
            gse_id: GEO Series accession (e.g. ``"GSE123456"``).

        Returns:
            A dict with keys ``gse_id``, ``entries`` (dict mapping
            ``filename`` -> ``md5``), and ``count``.
        """
        url = _geo_filelist_url(gse_id)
        self._log.info("Reading filelist.txt", gse_id=gse_id, url=url)

        try:
            async with httpx.AsyncClient(timeout=self.timeout) as client:
                resp = await client.get(url)
        except httpx.TimeoutException:
            return _error(gse_id, f"Timeout reading filelist.txt (>{self.timeout}s)")
        except httpx.RequestError as exc:
            return _error(gse_id, f"HTTP request failed: {exc}")

        if resp.status_code != 200:
            return _error(
                gse_id,
                f"GEO FTP returned HTTP {resp.status_code} for filelist.txt",
            )

        entries: dict[str, str] = {}
        for line in resp.text.strip().split("\n"):
            line = line.strip()
            if not line:
                continue
            # Format: "md5  filename"  (two space-separated tokens)
            parts = line.split(None, 1)  # split on any whitespace, max 2 parts
            if len(parts) == 2:
                md5_hash, filename = parts
                entries[filename.strip()] = md5_hash.strip()
            elif len(parts) == 1:
                # Edge case: line without md5 -- store with empty md5
                entries[parts[0].strip()] = ""

        self._log.info("Filelist read complete", gse_id=gse_id, count=len(entries))
        return {"gse_id": gse_id, "entries": entries, "count": len(entries)}


# ---------------------------------------------------------------------------
# GeoFileDownloader -- download GEO supplementary files
# ---------------------------------------------------------------------------


class GeoFileDownloader:
    """Downloads supplementary files from GEO FTP.

    Respects size limits and allows filtering by file type / extension.
    """

    def __init__(
        self,
        timeout: float = 120.0,
        chunk_size: int = 1024 * 1024,  # 1 MiB
    ) -> None:
        self.timeout = timeout
        self.chunk_size = chunk_size
        self._log = logger.bind(component="GeoFileDownloader")

    # -- public API ------------------------------------------------------

    async def download(
        self,
        gse_id: str,
        data_dir: str | Path,
        *,
        file_types: list[str] | None = None,
        max_size_mb: int = 500,
        download_method: str = "auto",
    ) -> dict[str, Any]:
        """Download supplementary files for a GSE dataset.

        Steps:
        1. List supplementary files via ``GeoSupplementaryLister``.
        2. Optionally filter by extension / type (``file_types``).
        3. Skip files whose reported size exceeds ``max_size_mb``.
        4. Download remaining files to ``data_dir`` via the selected backend.

        Args:
            gse_id: GEO Series accession.
            data_dir: Local directory to write downloaded files.
            file_types: If set, only download files whose extension (lowercase,
                without dot) is in this list.  Example: ``["h5ad", "h5", "loom",
                "mtx", "tsv", "gz", "txt"]``.
            max_size_mb: Skip any single file larger than this (megabytes).
            download_method: Backend selection — ``"auto"`` (smart selection),
                ``"aspera"``, ``"aria2"``, or ``"httpx"``.  Default: ``"auto"``.

        Returns:
            A dict with keys ``gse_id``, ``data_dir``, ``downloaded_files``
            (list of local paths), ``skipped_files`` (list of reasons),
            ``total_size_mb``, ``backend_used``, and ``error`` if something
            went wrong at the listing stage.
        """
        data_dir = Path(data_dir)
        data_dir.mkdir(parents=True, exist_ok=True)

        # 1. List
        lister = GeoSupplementaryLister(timeout=30.0)
        listing = await lister.list_supplementary(gse_id)
        if "error" in listing:
            return {**listing, "data_dir": str(data_dir), "downloaded_files": [],
                    "skipped_files": [], "total_size_mb": 0.0}

        files = listing.get("files", [])

        # 2. Filter by type
        if file_types:
            file_types_lower = {t.lower().lstrip(".") for t in file_types}
            files = [
                f for f in files
                if _ext(f["filename"]).lower() in file_types_lower
            ]
            self._log.info(
                "Filtered by type", gse_id=gse_id, remaining=len(files),
                types=sorted(file_types_lower),
            )

        # 3. Size filter
        max_size_bytes = max_size_mb * 1024 * 1024
        skipped: list[dict[str, str]] = []
        to_download: list[dict[str, Any]] = []
        for f in files:
            if f.get("size_bytes") and f["size_bytes"] > max_size_bytes:
                skipped.append({
                    "filename": f["filename"],
                    "reason": f"Exceeds max_size_mb ({max_size_mb} MB): {f['size_bytes'] / 1024 / 1024:.1f} MB",
                })
            else:
                to_download.append(f)

        self._log.info(
            "Download plan", gse_id=gse_id, to_download=len(to_download),
            skipped=len(skipped), download_method=download_method,
        )

        # 4. Download — use accelerator when available
        downloaded: list[str] = []
        total_bytes = 0
        backend_used = "httpx"

        # Try to set up the accelerator
        accelerator = None
        if download_method != "httpx":
            try:
                from gse_router.tools.download_accelerator import (
                    DownloadAccelerator, DownloadBackend,
                )
                from gse_router.tools.download_config import load_download_config
                dl_config = load_download_config()
                if download_method != "auto":
                    dl_config.backend = DownloadBackend(download_method)
                accelerator = DownloadAccelerator(
                    config=dl_config,
                    data_dir=data_dir,
                    max_concurrent=3,
                )
                await accelerator.detect_available_backends()
                backend_used = download_method
            except Exception as exc:
                self._log.info("Accelerator not available, using httpx", error=str(exc))

        if accelerator is not None:
            # Download via accelerator (Aspera/Aria2 with httpx fallback)
            for f_info in to_download:
                filename = f_info["filename"]
                url = f_info["url"]
                dest = data_dir / filename

                try:
                    self._log.info("Downloading (accelerator)", filename=filename)
                    result = await accelerator.download_file(
                        url=url, dest=dest,
                        expected_size=f_info.get("size_bytes"),
                    )
                    file_size = result.stat().st_size
                    total_bytes += file_size
                    downloaded.append(str(result))
                    self._log.info("Downloaded", filename=filename, size_mb=file_size / 1024 / 1024)
                except Exception as exc:
                    if dest.exists():
                        dest.unlink()
                    skipped.append({
                        "filename": filename,
                        "reason": f"Download error: {exc}",
                    })
                    self._log.warning("Download failed", filename=filename, error=str(exc))
        else:
            # Fallback: direct httpx streaming
            async with httpx.AsyncClient(timeout=self.timeout) as client:
                for f_info in to_download:
                    filename = f_info["filename"]
                    url = f_info["url"]
                    dest = data_dir / filename

                    try:
                        self._log.info("Downloading", filename=filename)
                        async with client.stream("GET", url) as resp:
                            if resp.status_code != 200:
                                skipped.append({
                                    "filename": filename,
                                    "reason": f"HTTP {resp.status_code}",
                                })
                                continue
                            with open(dest, "wb") as fh:
                                async for chunk in resp.aiter_bytes(chunk_size=self.chunk_size):
                                    fh.write(chunk)
                        file_size = dest.stat().st_size
                        total_bytes += file_size
                        downloaded.append(str(dest))
                        self._log.info("Downloaded", filename=filename, size_mb=file_size / 1024 / 1024)
                    except Exception as exc:
                        if dest.exists():
                            dest.unlink()
                        skipped.append({
                            "filename": filename,
                            "reason": f"Download error: {exc}",
                        })
                        self._log.warning("Download failed", filename=filename, error=str(exc))

        total_mb = total_bytes / 1024 / 1024
        self._log.info(
            "Download complete", gse_id=gse_id, downloaded=len(downloaded),
            skipped=len(skipped), total_size_mb=round(total_mb, 2),
            backend_used=backend_used,
        )

        return {
            "gse_id": gse_id,
            "data_dir": str(data_dir),
            "downloaded_files": downloaded,
            "skipped_files": skipped,
            "total_size_mb": round(total_mb, 2),
            "backend_used": backend_used,
        }


# ---------------------------------------------------------------------------
# MCP tool functions -- async, return structured dicts
# ---------------------------------------------------------------------------


async def fetch_geo_metadata(gse_id: str, email: str | None = None) -> dict[str, Any]:
    """Query GEO metadata and return structured study context.

    Wraps ``MetadataProbe.probe_all()`` from ``gse_router.probes.metadata``.

    Args:
        gse_id: GEO Series accession (e.g. ``"GSE123456"``).
        email: Optional NCBI contact email for E-utilities rate-limiting.

    Returns:
        A dict with keys ``gse_id``, ``title``, ``organism``, ``platform``,
        ``technology``, ``sample_count``, ``condition_groups``,
        ``sra_accessions``, ``confidence``, and ``error`` on failure.
    """
    log = logger.bind(tool="fetch_geo_metadata", gse_id=gse_id)
    log.info("Fetching GEO metadata")

    try:
        from gse_router.config import get_config
        from gse_router.probes.metadata import MetadataProbe

        cfg = get_config()
        resolved_email = email or cfg.ncbi_email or None
        resolved_api_key = cfg.ncbi_api_key or None

        probe = MetadataProbe(email=resolved_email, api_key=resolved_api_key)
        ctx = await probe.probe_all(gse_id)

        # Collect SRA accessions from samples
        sra_accessions: list[str] = []
        for sample in ctx.samples:
            sra_accessions.extend(sample.sra_runs)

        result: dict[str, Any] = {
            "gse_id": ctx.dataset_id,
            "title": ctx.title,
            "organism": ctx.organism,
            "platform": ctx.platform,
            "technology": ctx.technology,
            "sample_count": len(ctx.samples),
            "condition_groups": ctx.condition_groups,
            "sra_accessions": sorted(set(sra_accessions)),
            "confidence": ctx.confidence,
            "tissue_or_organ": ctx.tissue_or_organ,
            "disease": ctx.disease,
            "publication": ctx.publication,
            "probed_at": ctx.probed_at,
        }
        log.info("Metadata fetch complete", sample_count=len(ctx.samples))
        return result

    except Exception as exc:
        log.error("Metadata fetch failed", error=str(exc))
        return {"gse_id": gse_id, "error": str(exc)}


async def list_geo_supplementary(gse_id: str) -> dict[str, Any]:
    """List supplementary files available on GEO FTP for a dataset.

    Uses ``GeoSupplementaryLister.list_supplementary()``.

    Args:
        gse_id: GEO Series accession.

    Returns:
        A dict with ``gse_id``, ``files`` (list of
        ``{filename, size_bytes, url}``), ``count``, and ``error`` on failure.
    """
    log = logger.bind(tool="list_geo_supplementary", gse_id=gse_id)
    log.info("Listing GEO supplementary files")

    try:
        lister = GeoSupplementaryLister()
        return await lister.list_supplementary(gse_id)
    except Exception as exc:
        log.error("Supplementary listing failed", error=str(exc))
        return {"gse_id": gse_id, "error": str(exc), "files": [], "count": 0}


async def read_geo_filelist(gse_id: str) -> dict[str, Any]:
    """Read and parse the GEO supplementary ``filelist.txt`` manifest.

    Uses ``GeoSupplementaryLister.read_filelist()``.

    Args:
        gse_id: GEO Series accession.

    Returns:
        A dict with ``gse_id``, ``entries`` (dict mapping filename -> md5),
        ``count``, and ``error`` on failure.
    """
    log = logger.bind(tool="read_geo_filelist", gse_id=gse_id)
    log.info("Reading GEO filelist.txt")

    try:
        lister = GeoSupplementaryLister()
        return await lister.read_filelist(gse_id)
    except Exception as exc:
        log.error("Filelist read failed", error=str(exc))
        return {"gse_id": gse_id, "error": str(exc), "entries": {}, "count": 0}


async def download_geo_supplementary(
    gse_id: str,
    data_dir: str,
    file_types: list[str] | None = None,
    max_size_mb: int = 500,
    download_method: str = "auto",
) -> dict[str, Any]:
    """Download supplementary files for a GSE dataset.

    Uses ``GeoFileDownloader.download()``.

    Args:
        gse_id: GEO Series accession.
        data_dir: Local directory to write downloaded files.
        file_types: Optional list of file extensions to filter by (e.g.
            ``["h5ad", "h5", "loom", "mtx", "tsv", "gz"]``).
        max_size_mb: Maximum per-file size in megabytes (default 500).
        download_method: ``"auto"`` | ``"aspera"`` | ``"aria2"`` | ``"httpx"``.

    Returns:
        A dict with ``gse_id``, ``data_dir``, ``downloaded_files`` (list of
        local paths), ``skipped_files``, ``total_size_mb``, ``backend_used``,
        and ``error`` on failure.
    """
    log = logger.bind(tool="download_geo_supplementary", gse_id=gse_id)
    log.info("Downloading GEO supplementary files", data_dir=data_dir,
             download_method=download_method)

    try:
        downloader = GeoFileDownloader()
        return await downloader.download(
            gse_id=gse_id,
            data_dir=data_dir,
            file_types=file_types,
            max_size_mb=max_size_mb,
            download_method=download_method,
        )
    except Exception as exc:
        log.error("Download failed", error=str(exc))
        return {
            "gse_id": gse_id,
            "data_dir": data_dir,
            "error": str(exc),
            "downloaded_files": [],
            "skipped_files": [],
            "total_size_mb": 0.0,
        }


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _error(gse_id: str, message: str) -> dict[str, Any]:
    """Build a standardised error dict for a GeoSupplementaryLister call."""
    return {"gse_id": gse_id, "error": message, "files": [], "count": 0}


def _ext(filename: str) -> str:
    """Return the file extension (without dot), lowercased."""
    _, dot, ext = filename.rpartition(".")
    return ext if dot else ""


def _parse_ftp_size(raw: str) -> int | None:
    """Parse a human-readable size string from the NCBI FTP index.

    Examples:
        "1.2M" -> 1_258_291
        "340K" -> 348_160
        "5.3G" -> 5_690_687_488
        "1024" -> 1024
    """
    if not raw:
        return None
    raw = raw.strip().upper()
    multipliers = {"K": 1024, "M": 1024**2, "G": 1024**3, "T": 1024**4}
    try:
        if raw[-1] in multipliers:
            number = float(raw[:-1])
            return int(number * multipliers[raw[-1]])
        return int(float(raw))
    except (ValueError, IndexError):
        return None


# ---------------------------------------------------------------------------
# CLI entry point -- callable as ``python -m gse_router_mcp.tools.geo_metadata``
# ---------------------------------------------------------------------------


def _build_parser() -> argparse.ArgumentParser:
    import argparse

    parser = argparse.ArgumentParser(
        prog="python -m gse_router_mcp.tools.geo_metadata",
        description="GEO metadata tools -- fetch, list, and download.",
    )
    sub = parser.add_subparsers(dest="command", required=True)

    # fetch
    p_fetch = sub.add_parser("fetch", help="Fetch GEO metadata for a GSE accession")
    p_fetch.add_argument("gse_id", help="GEO Series accession (e.g. GSE123456)")
    p_fetch.add_argument("--email", default=None, help="NCBI contact email for E-utilities")
    p_fetch.add_argument("--pretty", action="store_true", default=True,
                         help="Pretty-print JSON output (default)")

    # list-supplementary
    p_list = sub.add_parser("list-supplementary", help="List GEO supplementary files")
    p_list.add_argument("gse_id", help="GEO Series accession")
    p_list.add_argument("--pretty", action="store_true", default=True)

    # read-filelist
    p_fl = sub.add_parser("read-filelist", help="Read GEO filelist.txt")
    p_fl.add_argument("gse_id", help="GEO Series accession")
    p_fl.add_argument("--pretty", action="store_true", default=True)

    # download
    p_dl = sub.add_parser("download", help="Download GEO supplementary files")
    p_dl.add_argument("gse_id", help="GEO Series accession")
    p_dl.add_argument("data_dir", help="Local directory for downloaded files")
    p_dl.add_argument("--file-types", nargs="*", default=None,
                      help="Only download files with these extensions (e.g. h5ad loom mtx)")
    p_dl.add_argument("--max-size-mb", type=int, default=500,
                      help="Skip files larger than this (default: 500)")
    p_dl.add_argument("--pretty", action="store_true", default=True)

    return parser


def main() -> None:
    """CLI entry point -- dispatched by micromamba commands in config/tools.yaml."""
    import argparse

    parser = _build_parser()
    args = parser.parse_args()

    async def _run() -> dict[str, Any]:
        if args.command == "fetch":
            return await fetch_geo_metadata(args.gse_id, email=args.email)
        elif args.command == "list-supplementary":
            return await list_geo_supplementary(args.gse_id)
        elif args.command == "read-filelist":
            return await read_geo_filelist(args.gse_id)
        elif args.command == "download":
            return await download_geo_supplementary(
                gse_id=args.gse_id,
                data_dir=args.data_dir,
                file_types=args.file_types,
                max_size_mb=args.max_size_mb,
            )
        else:
            return {"error": f"Unknown command: {args.command}"}

    result = asyncio.run(_run())

    indent = 2 if getattr(args, "pretty", True) else None
    print(json.dumps(result, indent=indent, default=str))
    sys.exit(0 if "error" not in result else 1)


if __name__ == "__main__":
    main()
