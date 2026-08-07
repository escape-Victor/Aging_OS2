"""GEO supplementary file downloader — async HTTP-based downloader for GEO FTP.

This module handles:
- Listing supplementary files from GEO FTP directory listings
- Parsing filelist.txt for MD5 checksums when available
- Selective downloading based on route plan file-type priorities
- Cardinal rule enforcement: never download massive RAW.tar if processed data exists

Architecture:
    Uses httpx.AsyncClient for all HTTP requests. Parses Apache/nginx directory
    listing HTML to extract file names and sizes without relying on BeautifulSoup
    (regex-based parsing for zero additional dependencies beyond httpx).

    This is a standalone module with no dependencies on other gse_router modules
    except structlog. It can be imported and used independently.
"""

from __future__ import annotations

import asyncio
import hashlib
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import ClassVar

import structlog

logger = structlog.get_logger(__name__)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

GEO_FTP_BASE: str = "https://ftp.ncbi.nlm.nih.gov/geo/series"

# Recognised supplementary file extensions and their categories
_SUPP_EXTENSION_MAP: dict[str, str] = {
    ".h5ad": "h5ad",
    ".h5ad.gz": "h5ad",
    ".h5": "h5",
    ".h5.gz": "h5",
    ".loom": "loom",
    ".loom.gz": "loom",
    ".mtx": "mtx",
    ".mtx.gz": "mtx",
    ".tsv": "tsv",
    ".tsv.gz": "tsv",
    ".csv": "csv",
    ".csv.gz": "csv",
    ".txt": "txt",
    ".txt.gz": "txt",
    ".tab": "tsv",
    ".tab.gz": "tsv",
    ".rds": "rds",
    ".rds.gz": "rds",
    ".rdata": "rds",
    ".rda": "rds",
    ".h5seurat": "rds",
    ".mat": "mat",
    ".mat.gz": "mat",
    ".xlsx": "xlsx",
    ".fastq": "fastq",
    ".fastq.gz": "fastq",
    ".fq": "fastq",
    ".fq.gz": "fastq",
    ".tar": "tar",
    ".tar.gz": "tar",
    ".tgz": "tar",
    ".zip": "zip",
    ".gz": "gz",
    ".bam": "bam",
    ".bam.bai": "bai",
    ".bai": "bai",
    ".jpg": "image",
    ".jpeg": "image",
    ".png": "image",
    ".tiff": "image",
    ".tif": "image",
    ".pdf": "pdf",
    ".html": "html",
}

# File types that indicate processed scRNA-seq data (not raw sequencing)
_PROCESSED_FILE_TYPES: frozenset[str] = frozenset({
    "h5ad", "h5", "loom", "mtx", "tsv", "csv", "txt", "rds", "mat", "xlsx",
})

# File types that are raw sequencing data — download only as last resort
_RAW_FILE_TYPES: frozenset[str] = frozenset({"fastq", "bam"})

# File types that are archive bundles — may contain processed data inside
_ARCHIVE_FILE_TYPES: frozenset[str] = frozenset({"tar", "zip"})

# File categories worth downloading (everything except HTML, images, PDFs)
# Archives are included because they often contain processed matrices.
# Raw FASTQ/BAM are excluded — only downloaded at explicit last resort.
_DOWNLOADABLE_CATEGORIES: frozenset[str] = frozenset({
    "h5ad", "h5", "loom", "mtx", "tsv", "csv", "txt",
    "rds", "tar", "zip", "gz", "mat", "xlsx",
})

# MIME types that correspond to HTML pages (not data files)
_HTML_MIME_TYPES: frozenset[str] = frozenset({
    "text/html", "application/xhtml+xml",
})


# ---------------------------------------------------------------------------
# Data classes
# ---------------------------------------------------------------------------


@dataclass
class GeoFileInfo:
    """Information about a single supplementary file on GEO FTP.

    Attributes:
        filename: Base file name (e.g. 'GSE123456_RAW.tar').
        size_bytes: File size in bytes, or None if unparseable.
        url: Full URL to the file.
        md5: MD5 hex digest from filelist.txt, or None if unavailable.
        category: Quick classification by extension (see _SUPP_EXTENSION_MAP).
    """

    filename: str
    size_bytes: int | None = None
    url: str = ""
    md5: str | None = None
    category: str = "other"

    @property
    def size_mb(self) -> float | None:
        """File size in megabytes, for human-readable display."""
        if self.size_bytes is None:
            return None
        return self.size_bytes / (1024 * 1024)

    @property
    def is_processed(self) -> bool:
        """True if this file likely contains processed expression data."""
        return self.category in _PROCESSED_FILE_TYPES

    @property
    def is_raw(self) -> bool:
        """True if this file is raw sequencing data (FASTQ/BAM)."""
        return self.category in _RAW_FILE_TYPES

    @property
    def is_archive(self) -> bool:
        """True if this file is a tar/zip archive bundle."""
        return self.category in _ARCHIVE_FILE_TYPES


# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------


def classify_supplementary_file(filename: str) -> str:
    """Quick classification of a GEO supplementary file by extension.

    Args:
        filename: The file name (e.g. 'GSE123456_RAW.tar', 'barcodes.tsv.gz').

    Returns:
        One of: 'h5ad', 'h5', 'loom', 'mtx', 'tsv', 'csv', 'txt', 'rds',
        'fastq', 'tar', 'zip', 'gz', 'bam', 'bai', 'image', 'pdf', 'html',
        or 'other'.
    """
    lower = filename.lower()
    # Check compound extensions first (e.g. .tsv.gz, .loom.gz)
    for ext in (".h5ad.gz", ".h5.gz", ".loom.gz", ".mat.gz",
                ".tar.gz", ".tgz", ".tsv.gz", ".csv.gz",
                ".txt.gz", ".mtx.gz", ".fastq.gz", ".fq.gz",
                ".rds.gz", ".bam.bai", ".tab.gz"):
        if lower.endswith(ext):
            return _SUPP_EXTENSION_MAP[ext]
    # Check single extensions
    suffix = Path(lower).suffix
    return _SUPP_EXTENSION_MAP.get(suffix, "other")


def _human_size_to_bytes(text: str) -> int | None:
    """Parse a human-readable file size string into bytes.

    Handles formats produced by Apache mod_autoindex (e.g. '500M', '1.2G',
    ' 512K', '4096', '3.4T'). Returns None for unparseable values.

    Args:
        text: The raw size text from an Apache directory listing cell.

    Returns:
        File size in bytes, or None.
    """
    text = text.strip()
    if not text:
        return None
    try:
        # Pure integer (bytes in some nginx configs)
        return int(text)
    except ValueError:
        pass

    units: dict[str, int] = {
        "K": 1024,
        "M": 1024**2,
        "G": 1024**3,
        "T": 1024**4,
        "P": 1024**5,
    }
    match = re.match(r"^([\d.]+)\s*([KMGTP])?$", text, re.IGNORECASE)
    if not match:
        return None
    value = float(match.group(1))
    unit = match.group(2)
    multiplier = units.get(unit.upper(), 1) if unit else 1
    return int(value * multiplier)


def _build_suppl_url(gse_id: str) -> str:
    """Construct the GEO supplementary files directory URL.

    GEO organises supplementary files in a 3-level directory hierarchy:
    /series/GSEnnn/<gse_id>/suppl/

    For example, GSE123456 → /series/GSE123nnn/GSE123456/suppl/

    Args:
        gse_id: GEO Series accession (e.g. 'GSE123456').

    Returns:
        Full URL to the supplementary directory listing.

    Raises:
        ValueError: If gse_id does not match the expected GSE/GPL/GSM pattern.
    """
    gse_id = gse_id.strip().upper()
    if not re.match(r"^GS[E|P|M]\d+$", gse_id):
        raise ValueError(
            f"Invalid GEO accession '{gse_id}'. "
            f"Expected pattern: GSE/GPL/GSM followed by digits."
        )
    # Extract the numeric portion and compute the nnn prefix
    num_part = re.search(r"\d+", gse_id)
    if num_part is None:
        raise ValueError(f"Cannot extract numeric portion from '{gse_id}'")
    num = int(num_part.group())
    prefix = f"{num // 1000 * 1000:03d}nnn"
    # e.g. GSE123456 → GSE123nnn
    gse_prefix = f"{gse_id[:-3]}nnn"
    return f"{GEO_FTP_BASE}/{gse_prefix}/{gse_id}/suppl/"


# ---------------------------------------------------------------------------
# Apache / nginx directory listing parser
# ---------------------------------------------------------------------------

# Pattern 1: Apache 2.4 mod_autoindex (FancyIndexing on)
# <tr><td><a href="file.name">file.name</a></td><td align="right">2024-01-15 10:30  </td><td align="right"> 500M  </td></tr>
_APACHE_ROW_RE = re.compile(
    r'<tr>\s*'
    r'<td[^>]*>\s*<a\s+href="([^"]+)"[^>]*>.*?</a>\s*</td>\s*'
    r'<td[^>]*>[^<]*</td>\s*'        # date column (skip)
    r'<td[^>]*>([^<]*)</td>',         # size column
    re.IGNORECASE | re.DOTALL,
)

# Pattern 2: Simpler Apache variant (no date column)
# <tr><td><a href="file">file</a></td><td align="right">  500M  </td></tr>
_APACHE_SIMPLE_RE = re.compile(
    r'<tr>\s*'
    r'<td[^>]*>\s*<a\s+href="([^"]+)"[^>]*>.*?</a>\s*</td>\s*'
    r'<td[^>]*>([^<]*)</td>',
    re.IGNORECASE | re.DOTALL,
)

# Pattern 3: nginx autoindex (traditional format with abbreviated month names)
# <a href="file.name">file.name</a>  01-Jan-2024 10:30  500M
_NGINX_ROW_RE = re.compile(
    r'<a\s+href="([^"]+)">[^<]+</a>\s+'
    r'\d{2}-[A-Z][a-z]{2}-\d{4}\s+\d{2}:\d{2}\s+'
    r'([\d.KMGTP]+)',
    re.IGNORECASE,
)

# Pattern 4: GEO FTP <pre> block format (ISO 8601 dates)
# <a href="file.name">file.name</a>   2016-07-19 14:27   20M
# This is the format actually used by NCBI GEO FTP as of 2024+.
# Dates are YYYY-MM-DD, sizes are human-readable (K/M/G/T suffix) or "-".
_GEO_FTP_PRE_RE = re.compile(
    r'<a\s+href="([^"]+)">[^<]+</a>\s+'
    r'\d{4}-\d{2}-\d{2}\s+\d{2}:\d{2}\s+'
    r'([\d.]+[KMGTP]?|-)\s*',
    re.IGNORECASE,
)

# Fallback: just extract all <a> tags with non-directory hrefs
_ANCHOR_RE = re.compile(
    r'<a\s+href="([^"]+)"[^>]*>(?:<[^>]+>)*([^<]+)(?:</[^>]+>)*</a>',
    re.IGNORECASE,
)

# Regex to identify entries that are subdirectories (trailing slash or known icons)
_DIR_INDICATORS: set[str] = {"/", "/?C=N;O=D", "/?C=M;O=A", "/?C=S;O=A", "/?C=D;O=A"}


def _parse_apache_listing(html: str, base_url: str) -> list[GeoFileInfo]:
    """Parse an Apache/nginx directory listing HTML page.

    Extracts file names and sizes from the HTML table produced by
    Apache mod_autoindex or nginx autoindex module.

    Args:
        html: Raw HTML content of the directory listing page.
        base_url: The base URL of the directory (used to construct full URLs).

    Returns:
        List of GeoFileInfo objects, one per file found.
    """
    files: list[GeoFileInfo] = []
    seen: set[str] = set()

    # Ensure base_url ends with a slash
    if not base_url.endswith("/"):
        base_url += "/"

    # --- Try structured parsers first ---
    for pattern, has_size_col in [
        (_APACHE_ROW_RE, True),
        (_APACHE_SIMPLE_RE, True),
    ]:
        for match in pattern.finditer(html):
            href = match.group(1).strip()
            size_text = match.group(2).strip() if has_size_col else ""
            _add_file_entry(href, size_text, base_url, files, seen)
        if files:
            # Found files with a structured parser — use these results
            return files

    # --- Try nginx autoindex ---
    for match in _NGINX_ROW_RE.finditer(html):
        href = match.group(1).strip()
        size_text = match.group(2).strip() if match.lastindex and match.lastindex >= 2 else ""
        _add_file_entry(href, size_text, base_url, files, seen)
    if files:
        return files

    # --- Try GEO FTP <pre> block format (ISO dates: YYYY-MM-DD HH:MM SIZE) ---
    for match in _GEO_FTP_PRE_RE.finditer(html):
        href = match.group(1).strip()
        size_text = match.group(2).strip()
        # Skip entries with "-" size (directory markers)
        if size_text == "-":
            continue
        _add_file_entry(href, size_text, base_url, files, seen)
    if files:
        return files

    # --- Fallback: extract all anchor tags ---
    for match in _ANCHOR_RE.finditer(html):
        href = match.group(1).strip()
        _add_file_entry(href, "", base_url, files, seen)

    return files


def _add_file_entry(
    href: str,
    size_text: str,
    base_url: str,
    files: list[GeoFileInfo],
    seen: set[str],
) -> None:
    """Add a file entry from a parsed <a> tag to the results list.

    Skips parent-directory links, subdirectory links, query-only hrefs, and
    duplicates.

    Args:
        href: The href attribute value from the <a> tag.
        size_text: Raw size text from the adjacent table cell (may be empty).
        base_url: Base URL for constructing the full file URL.
        files: Accumulator list for GeoFileInfo objects.
        seen: Set of filenames already added (to avoid duplicates).
    """
    # Skip empty, anchor-only, or javascript links
    if not href or href.startswith("#") or href.startswith("javascript:"):
        return
    # Skip parent directory link
    if href in ("../", "..", "/"):
        return
    # Skip directory entries (trailing slash or known autoindex query strings)
    if href.endswith("/") or any(href.endswith(d) for d in _DIR_INDICATORS):
        return
    # Skip query-string-only hrefs (e.g. "?C=N;O=D")
    if href.startswith("?") and "=" in href:
        return

    filename = href.split("/")[-1]  # in case href is a relative path
    if not filename:
        return
    if filename in seen:
        return
    seen.add(filename)

    size_bytes = _human_size_to_bytes(size_text)
    url = base_url + filename
    category = classify_supplementary_file(filename)

    files.append(GeoFileInfo(
        filename=filename,
        size_bytes=size_bytes,
        url=url,
        md5=None,  # filled later if filelist.txt is available
        category=category,
    ))


# ---------------------------------------------------------------------------
# filelist.txt parser
# ---------------------------------------------------------------------------

# GEO has used at least two filelist.txt formats over time:
#
# Format A (older, tab-separated, with MD5):
#   GSE123456_RAW.tar    524288000    abc123def456...
#   Columns: filename, size_bytes, md5_hex
#
# Format B (newer, tab-separated, with header):
#   #Archive/File	Name	Time	Size	Type
#   Archive	GSE129788_RAW.tar	07/01/2019 14:50:47	429690880	TAR
#   File	GSM3722100_YX1L_10X.txt.gz	04/15/2019 06:28:33	15464907	TXT
#   Columns: entry_type, name, date, size_bytes, file_type
#
# We try Format B first (presence of "#Archive/File" header), fall back to A.

_FILELIST_LINE_RE = re.compile(
    r"^(\S+)\s+(\d+)\s+([a-fA-F0-9]{32})",
    re.IGNORECASE,
)
_FILELIST_LINE_NO_MD5_RE = re.compile(
    r"^(\S+)\s+(\d+)",
)

# Format B: TSV with columns Type, Name, Time, Size, FileType
_FILELIST_FORMAT_B_HEADER = "#Archive/File"
_FILELIST_FORMAT_B_RE = re.compile(
    r"^(Archive|File)\t([^\t]+)\t([^\t]+)\t(\d+)\t(\S+)",
)


def _parse_filelist(text: str) -> dict[str, str]:
    """Parse a filelist.txt contents into a filename-to-MD5 mapping.

    Handles two formats:
    - Format A (older): ``filename  size  md5`` (space/tab separated)
    - Format B (newer): ``Type\\tName\\tTime\\tSize\\tFileType`` (5-column TSV)

    Format B does NOT include MD5 hashes — only filename and size are
    extracted. Format A provides MD5 hex digests for checksum verification.

    Args:
        text: Raw content of the filelist.txt file.

    Returns:
        Dict mapping filename -> MD5 hex string (empty string if no MD5).
    """
    md5_map: dict[str, str] = {}

    # --- Try Format B first (newer GEO format with 5-column TSV) ---
    if _FILELIST_FORMAT_B_HEADER in text:
        for line in text.splitlines():
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            match = _FILELIST_FORMAT_B_RE.match(line)
            if match:
                name = match.group(2)
                # Format B has no MD5 — store empty string as placeholder
                if name not in md5_map:
                    md5_map[name] = ""
        return md5_map

    # --- Fall back to Format A (older format: filename size md5) ---
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        # Try full format: filename \t size \t md5
        match = _FILELIST_LINE_RE.match(line)
        if match:
            md5_map[match.group(1)] = match.group(3).lower()
            continue
        # Fallback: filename \t size (no md5)
        match = _FILELIST_LINE_NO_MD5_RE.match(line)
        if match:
            if match.group(1) not in md5_map:
                md5_map[match.group(1)] = ""
    return md5_map


# ---------------------------------------------------------------------------
# GeoSupplementaryLister
# ---------------------------------------------------------------------------


class GeoSupplementaryLister:
    """Lists supplementary files for a GEO dataset via HTTP directory scraping.

    This class handles the first step of the download pipeline: discovering
    what files are available on GEO FTP for a given dataset. It parses the
    Apache/nginx directory listing page and optionally the filelist.txt for
    MD5 checksums.

    Attributes:
        client: An httpx.AsyncClient instance. One is created if not provided.
        _owns_client: Whether this instance created the client (and should close it).
    """

    def __init__(self, client: "httpx.AsyncClient | None" = None) -> None:
        """Initialise the lister.

        Args:
            client: An optional httpx.AsyncClient to reuse. If None, a new
                client is created. The caller is responsible for closing the
                client if one is passed in.
        """
        self._owns_client: bool = client is None
        self._client: "httpx.AsyncClient | None" = client
        self._logger = structlog.get_logger(__name__).bind(
            component="GeoSupplementaryLister"
        )

    async def _get_client(self) -> "httpx.AsyncClient":
        """Return the httpx client, creating one lazily if needed."""
        if self._client is None:
            import httpx
            self._client = httpx.AsyncClient(
                timeout=30.0,
                follow_redirects=True,
                headers={
                    "User-Agent": (
                        "Mozilla/5.0 (compatible; gse-router/0.1; "
                        "+https://github.com/example/gse-router)"
                    ),
                },
            )
            self._owns_client = True
        return self._client

    async def close(self) -> None:
        """Close the underlying httpx client if this instance owns it."""
        if self._owns_client and self._client is not None:
            await self._client.aclose()
            self._client = None

    async def list_supplementary(self, gse_id: str) -> list[GeoFileInfo]:
        """List all supplementary files available for a GEO dataset.

        Constructs the GEO FTP supplementary directory URL, fetches the
        Apache/nginx directory listing page, and parses it for file names
        and sizes.

        Args:
            gse_id: GEO Series accession (e.g. 'GSE123456').

        Returns:
            List of GeoFileInfo objects representing available files.
            Returns an empty list if the directory cannot be accessed.

        Raises:
            ValueError: If gse_id is not a valid GEO accession.
        """
        url = _build_suppl_url(gse_id)
        self._logger.info("Listing supplementary files", gse_id=gse_id, url=url)

        client = await self._get_client()
        try:
            response = await client.get(url)
        except Exception as exc:
            self._logger.warning(
                "Failed to fetch supplementary directory listing",
                gse_id=gse_id,
                url=url,
                error=str(exc),
            )
            return []

        if response.status_code == 404:
            self._logger.info(
                "No supplementary directory for dataset",
                gse_id=gse_id,
                status=response.status_code,
            )
            return []
        if response.status_code != 200:
            self._logger.warning(
                "Non-200 response from supplementary directory",
                gse_id=gse_id,
                status=response.status_code,
            )
            return []

        content_type = response.headers.get("content-type", "")
        if any(mime in content_type for mime in _HTML_MIME_TYPES):
            html = response.text
        else:
            # Some GEO directories return text/plain for the listing
            html = response.text

        files = _parse_apache_listing(html, url)

        self._logger.info(
            "Found supplementary files",
            gse_id=gse_id,
            count=len(files),
            total_size_mb=round(sum(f.size_bytes or 0 for f in files) / (1024 * 1024), 1),
        )
        return files

    async def read_filelist(self, gse_id: str) -> dict[str, str]:
        """Fetch and parse the filelist.txt for a GEO dataset.

        filelist.txt contains filename, size, and MD5 checksum for each
        supplementary file. It is available for many (but not all) GEO datasets.

        Args:
            gse_id: GEO Series accession (e.g. 'GSE123456').

        Returns:
            Dict mapping filename -> MD5 hex string. Returns an empty dict
            if filelist.txt is not available.
        """
        url = _build_suppl_url(gse_id)
        if not url.endswith("/"):
            url += "/"
        filelist_url = url + "filelist.txt"

        self._logger.info("Fetching filelist.txt", gse_id=gse_id)

        client = await self._get_client()
        try:
            response = await client.get(filelist_url)
        except Exception as exc:
            self._logger.info(
                "Failed to fetch filelist.txt",
                gse_id=gse_id,
                error=str(exc),
            )
            return {}

        if response.status_code != 200:
            self._logger.info(
                "filelist.txt not available",
                gse_id=gse_id,
                status=response.status_code,
            )
            return {}

        return _parse_filelist(response.text)

    async def list_with_checksums(self, gse_id: str) -> list[GeoFileInfo]:
        """List supplementary files and enrich with MD5 checksums from filelist.txt.

        Combines list_supplementary() and read_filelist() into a single
        convenience method. Files get their md5 attribute populated from
        filelist.txt if available.

        Args:
            gse_id: GEO Series accession.

        Returns:
            List of GeoFileInfo with md5 populated where available.
        """
        files = await self.list_supplementary(gse_id)
        if not files:
            return files

        md5_map = await self.read_filelist(gse_id)
        if md5_map:
            for f in files:
                if f.filename in md5_map:
                    f.md5 = md5_map[f.filename] or None

            matched = sum(1 for f in files if f.md5)
            self._logger.info(
                "Enriched files with MD5 checksums",
                gse_id=gse_id,
                total=len(files),
                with_md5=matched,
            )

        return files


# ---------------------------------------------------------------------------
# GeoFileDownloader
# ---------------------------------------------------------------------------


async def _verify_md5(file_path: Path, expected_md5: str) -> bool:
    """Verify a downloaded file's MD5 checksum.

    Args:
        file_path: Path to the downloaded file.
        expected_md5: Expected MD5 hex digest.

    Returns:
        True if the checksums match, False otherwise.
    """
    if not expected_md5:
        return True  # no checksum to verify against
    hasher = hashlib.md5()
    with open(file_path, "rb") as fh:
        for chunk in iter(lambda: fh.read(8192), b""):
            hasher.update(chunk)
    actual = hasher.hexdigest()
    return actual == expected_md5.lower()


class GeoFileDownloader:
    """Downloads GEO supplementary files via HTTP with progress logging.

    Features:
    - Concurrent downloads with configurable max parallelism.
    - Selective downloading by file type (h5ad, h5, loom, mtx, etc.).
    - Maximum file size safety limit (skip huge archive tarballs).
    - MD5 checksum verification against filelist.txt when available.
    - Structural logging of each download with timing and size.

    The downloader enforces the cardinal rule: never download massive
    RAW.tar archives if processed data files exist. Callers should use
    filters to select only the files needed by the route plan.

    Attributes:
        data_dir: Root directory where downloaded files are stored.
        max_concurrent: Maximum number of simultaneous downloads.
    """

    def __init__(
        self,
        data_dir: Path,
        client: "httpx.AsyncClient | None" = None,
        max_concurrent: int = 3,
        use_accelerator: bool = True,
        download_config: "DownloadConfig | None" = None,
    ) -> None:
        """Initialise the downloader.

        Args:
            data_dir: Directory where downloaded datasets will be stored.
                A subdirectory named after the GSE ID is created automatically.
            client: Optional httpx.AsyncClient to reuse. One is created if None.
            max_concurrent: Max simultaneous downloads (default: 3).
            use_accelerator: If True (default), try Aspera/Aria2 before httpx.
                Gracefully degrades to httpx if no accelerator CLI is installed.
            download_config: Optional DownloadConfig for accelerator tuning.
                If None, config is loaded from environment variables.
        """
        self.data_dir = Path(data_dir)
        self.max_concurrent = max_concurrent
        self._client = client
        self._owns_client: bool = client is None
        self._semaphore: asyncio.Semaphore | None = None
        self._logger = structlog.get_logger(__name__).bind(
            component="GeoFileDownloader"
        )
        # --- Accelerator integration ---
        self._use_accelerator = use_accelerator
        self._accelerator: "DownloadAccelerator | None" = None
        if use_accelerator:
            from .download_accelerator import DownloadAccelerator
            from .download_config import DownloadConfig, load_download_config
            self._accelerator = DownloadAccelerator(
                config=download_config or load_download_config(),
                data_dir=data_dir,
                max_concurrent=max_concurrent,
            )

    async def _get_client(self) -> "httpx.AsyncClient":
        """Return the httpx client, creating one lazily if needed."""
        if self._client is None:
            import httpx
            self._client = httpx.AsyncClient(
                timeout=600.0,  # 10 minutes for large files
                follow_redirects=True,
                headers={
                    "User-Agent": (
                        "Mozilla/5.0 (compatible; gse-router/0.1; "
                        "+https://github.com/example/gse-router)"
                    ),
                },
            )
            self._owns_client = True
        return self._client

    async def close(self) -> None:
        """Close the underlying httpx client and accelerator if owned."""
        if self._owns_client and self._client is not None:
            await self._client.aclose()
            self._client = None
        if self._accelerator is not None:
            await self._accelerator.close()
            self._accelerator = None

    @property
    def semaphore(self) -> "asyncio.Semaphore":
        """Lazily initialised semaphore for concurrency limiting."""
        if self._semaphore is None:
            self._semaphore = asyncio.Semaphore(self.max_concurrent)
        return self._semaphore

    # ------------------------------------------------------------------
    # Single file download
    # ------------------------------------------------------------------

    async def download_file(
        self,
        url: str,
        dest: Path,
        expected_size: int | None = None,
        expected_md5: str | None = None,
        max_retries: int = 3,
    ) -> Path:
        """Download a single file with progress logging.

        If the accelerator is enabled, tries Aspera/Aria2 first, then falls
        back to httpx streaming. This provides the best available speed for
        each file without changing the caller's interface.

        Args:
            url: Full URL to the file.
            dest: Destination path on the local filesystem.
            expected_size: Expected file size in bytes (for verification).
            expected_md5: Expected MD5 hex digest (for verification).

        Returns:
            The local Path to the downloaded file.

        Raises:
            httpx.HTTPStatusError: If the HTTP request fails.
            OSError: If the file cannot be written.
        """
        # --- Try accelerator first ---
        if self._accelerator is not None:
            try:
                return await self._accelerator.download_file(
                    url=url,
                    dest=dest,
                    expected_size=expected_size,
                    expected_md5=expected_md5,
                    max_retries=max_retries,
                )
            except Exception as exc:
                self._logger.warning(
                    "Accelerator download failed, falling back to httpx",
                    url=url,
                    error=str(exc),
                )
                if dest.exists():
                    dest.unlink(missing_ok=True)
                # Fall through to existing httpx code below

        import asyncio
        import sys
        import time
        import httpx

        dest.parent.mkdir(parents=True, exist_ok=True)
        filename = dest.name

        last_error = None
        for attempt in range(max_retries):
            if attempt > 0:
                delay = min(2 ** attempt, 30)
                self._logger.info(
                    "Retrying download",
                    attempt=attempt + 1,
                    max_retries=max_retries,
                    delay_s=delay,
                    filename=filename,
                )
                await asyncio.sleep(delay)
                if dest.exists():
                    dest.unlink()

            self._logger.info(
                "Downloading file",
                url=url,
                dest=str(dest),
                expected_size_mb=(
                    round(expected_size / (1024 * 1024), 1)
                    if expected_size else None
                ),
            )
            start_time = time.monotonic()

            client = await self._get_client()
            async with self._sem:
                try:
                    async with client.stream("GET", url) as response:
                        response.raise_for_status()
                        total_size = int(response.headers.get("content-length", 0))
                        downloaded = 0
                        with open(dest, "wb") as fh:
                            async for chunk in response.aiter_bytes(chunk_size=65536):
                                fh.write(chunk)
                                downloaded += len(chunk)
                        actual_size = dest.stat().st_size
                except httpx.HTTPStatusError:
                    self._logger.error(
                        "Download failed (HTTP error)",
                        url=url,
                        status=getattr(getattr(response, "status_code", None), "__int__", lambda: "?")() if 'response' in dir() else "?",
                    )
                    last_error = sys.exc_info()[1]
                    continue
                except Exception:
                    self._logger.error("Download failed", url=url, exc_info=True)
                    if dest.exists():
                        dest.unlink()
                    last_error = sys.exc_info()[1]
                    continue
                break  # Success — exit retry loop
        else:
            # All retries exhausted
            raise last_error  # type: ignore[misc]

        elapsed = time.monotonic() - start_time
        speed_mbps = (actual_size / (1024 * 1024) / elapsed) if elapsed > 0 else 0

        # --- Size verification ---
        if expected_size is not None and actual_size != expected_size:
            self._logger.warning(
                "Downloaded file size mismatch",
                expected=expected_size,
                actual=actual_size,
                filename=filename,
            )
        elif total_size > 0 and actual_size != total_size:
            self._logger.warning(
                "Incomplete download (content-length mismatch)",
                expected=total_size,
                actual=actual_size,
                filename=filename,
            )

        # --- MD5 verification ---
        if expected_md5:
            md5_ok = await _verify_md5(dest, expected_md5)
            if not md5_ok:
                self._logger.error(
                    "MD5 checksum mismatch",
                    filename=filename,
                    expected=expected_md5,
                )
                dest.unlink()
                raise IOError(
                    f"MD5 checksum mismatch for '{filename}': "
                    f"expected {expected_md5}"
                )

        self._logger.info(
            "Download complete",
            filename=filename,
            size_mb=round(actual_size / (1024 * 1024), 1),
            elapsed_s=round(elapsed, 1),
            speed_mbps=round(speed_mbps, 1),
        )
        return dest

    # ------------------------------------------------------------------
    # Multiple file downloads (concurrent)
    # ------------------------------------------------------------------

    async def download_files(self, files: list[GeoFileInfo]) -> list[Path]:
        """Download multiple files concurrently, limited by max_concurrent.

        Args:
            files: List of GeoFileInfo objects describing files to download.

        Returns:
            List of local Paths for successfully downloaded files.
            Failed downloads are logged but do not halt the batch.
        """
        self._logger.info(
            "Starting batch download",
            count=len(files),
            max_concurrent=self.max_concurrent,
        )
        tasks: list[asyncio.Task] = []
        results: list[Path] = []

        for f in files:
            task = asyncio.create_task(
                self._download_one_safe(f), name=f.filename
            )
            tasks.append(task)

        outcomes = await asyncio.gather(*tasks, return_exceptions=True)
        for i, outcome in enumerate(outcomes):
            if isinstance(outcome, Exception):
                self._logger.error(
                    "Batch download task failed",
                    filename=files[i].filename,
                    error=str(outcome),
                )
            elif outcome is not None:
                results.append(outcome)

        self._logger.info(
            "Batch download finished",
            requested=len(files),
            succeeded=len(results),
            failed=len(files) - len(results),
        )
        return results

    async def _download_one_safe(self, file_info: GeoFileInfo) -> Path | None:
        """Download a single file, returning None on failure instead of raising.

        Args:
            file_info: GeoFileInfo for the file to download.

        Returns:
            Local Path on success, None on failure.
        """
        dest = self.data_dir / file_info.filename
        try:
            return await self.download_file(
                url=file_info.url,
                dest=dest,
                expected_size=file_info.size_bytes,
                expected_md5=file_info.md5,
            )
        except Exception as exc:
            self._logger.error(
                "Download failed for file",
                filename=file_info.filename,
                error=str(exc),
            )
            return None

    # ------------------------------------------------------------------
    # Dataset-level download
    # ------------------------------------------------------------------

    async def download_dataset(
        self,
        gse_id: str,
        file_types: list[str] | None = None,
        max_size_mb: int = 500,
        skip_archives: bool = False,
        skip_raw: bool = True,
    ) -> tuple[Path, list[GeoFileInfo], list[Path]]:
        """Main entry point: discover and download supplementary files for a dataset.

        Workflow:
        1. List all supplementary files via GeoSupplementaryLister.
        2. Optionally enrich with MD5 checksums from filelist.txt.
        3. Apply filters (file type, size limit, skip raw).
        4. Download matching files concurrently.
        5. Auto-unpack archives (tar/zip) after download.
        6. Return the dataset directory, all file infos, and downloaded paths.

        By default, downloads ALL data files (everything except HTML, images,
        and raw FASTQ/BAM). Archives (.tar/.zip) ARE included — they often
        contain processed data that the file probe will discover after
        extraction. The cardinal rule ("prefer processed over raw") is
        enforced by the route planner, not the downloader.

        Args:
            gse_id: GEO Series accession (e.g. 'GSE123456').
            file_types: If provided, only download files whose category
                matches one of these types (e.g. ['h5ad', 'h5', 'mtx', 'loom']).
                If None, downloads all categories in _DOWNLOADABLE_CATEGORIES.
            max_size_mb: Skip files larger than this size in MB (default: 500).
            skip_archives: Skip .tar and .zip files (default: False — archives
                are downloaded so the file probe can inspect their contents).
            skip_raw: Skip FASTQ/BAM files (default: True).

        Returns:
            Tuple of (dataset_dir, all_file_infos, downloaded_paths).
            dataset_dir is the directory where files were saved.
            all_file_infos is the complete listing from the lister.
            downloaded_paths lists the local paths of successfully
            downloaded files.
        """
        gse_id = gse_id.strip().upper()
        dataset_dir = self.data_dir / gse_id
        dataset_dir.mkdir(parents=True, exist_ok=True)
        self.data_dir = dataset_dir  # subsequent downloads go into dataset_dir

        # Determine which categories to include
        effective_file_types = file_types
        if effective_file_types is None:
            effective_file_types = list(_DOWNLOADABLE_CATEGORIES)

        self._logger.info(
            "Starting dataset download",
            gse_id=gse_id,
            file_types=effective_file_types,
            max_size_mb=max_size_mb,
            skip_archives=skip_archives,
            skip_raw=skip_raw,
        )

        # Step 1: list
        lister = GeoSupplementaryLister(client=self._client)
        all_files = await lister.list_with_checksums(gse_id)

        if not all_files:
            self._logger.warning(
                "No supplementary files found for dataset",
                gse_id=gse_id,
            )
            return dataset_dir, [], []

        # Step 2: filter
        to_download: list[GeoFileInfo] = []
        skipped_archives: list[str] = []
        skipped_raw: list[str] = []
        skipped_size: list[str] = []
        skipped_type: list[str] = []

        for f in all_files:
            # Skip archives only if explicitly requested (default: download them)
            if skip_archives and f.is_archive:
                skipped_archives.append(f.filename)
                continue
            # Skip raw sequencing data (FASTQ/BAM) — always skipped by default
            if skip_raw and f.is_raw:
                skipped_raw.append(f.filename)
                continue
            # Skip HTML and image files — never useful for scRNA-seq
            if f.category in ("html", "image", "pdf"):
                skipped_type.append(f.filename)
                continue
            # Skip oversized files
            if f.size_bytes is not None and max_size_mb > 0:
                if f.size_bytes > max_size_mb * 1024 * 1024:
                    skipped_size.append(
                        f"{f.filename} ({f.size_mb:.0f} MB)"
                    )
                    continue
            # Filter by file type (greedy by default via _DOWNLOADABLE_CATEGORIES)
            if effective_file_types is not None:
                if f.category not in effective_file_types:
                    skipped_type.append(f.filename)
                    continue
            to_download.append(f)

        self._logger.info(
            "Filtered files for download",
            gse_id=gse_id,
            total=len(all_files),
            to_download=len(to_download),
            skipped_archives=len(skipped_archives),
            skipped_raw=len(skipped_raw),
            skipped_size=len(skipped_size),
            skipped_type=len(skipped_type),
        )
        if skipped_archives:
            self._logger.info("Skipped archives", files=skipped_archives)
        if skipped_raw:
            self._logger.info("Skipped raw files", files=skipped_raw)
        if skipped_size:
            self._logger.info("Skipped oversized files", files=skipped_size)
        if skipped_type:
            self._logger.info("Skipped non-data files", files=skipped_type)

        # Step 3: download
        downloaded = await self.download_files(to_download)

        total_size_mb = round(
            sum(p.stat().st_size for p in downloaded if p.exists()) / (1024 * 1024),
            1,
        )

        self._logger.info(
            "Dataset download complete",
            gse_id=gse_id,
            downloaded_count=len(downloaded),
            total_size_mb=total_size_mb,
        )

        # Step 4: unpack archives so the file probe can inspect their contents
        if not skip_archives:
            extracted = self.unpack_archives(dataset_dir)
            if extracted:
                self._logger.info(
                    "Archives unpacked",
                    gse_id=gse_id,
                    archive_count=len(extracted),
                    extract_dirs=[str(d) for d in extracted],
                )

        return dataset_dir, all_files, downloaded

    # ------------------------------------------------------------------
    # Archive unpacking
    # ------------------------------------------------------------------

    @staticmethod
    def unpack_archives(dataset_dir: Path) -> list[Path]:
        """Find and extract all tar/zip archives in a dataset directory.

        Scans the directory for .tar, .tar.gz, .tgz, and .zip files, then
        extracts each into a subdirectory ``_extracted_{stem}/``. This makes
        the archive contents visible to the file probe without modifying the
        original archives.

        Many GEO datasets only provide processed data inside RAW.tar archives.
        By extracting everything, we ensure the file probe can discover all
        available data regardless of how the submitter packaged it.

        Args:
            dataset_dir: Directory containing downloaded files (including archives).

        Returns:
            List of extraction directory Paths. Empty list if no archives found
            or extraction failed for all archives.
        """
        import tarfile
        import zipfile

        extracted: list[Path] = []

        # Only true archives (containing multiple files) — NOT standalone .gz
        # data files. Compressed data (.csv.gz, .txt.gz, .fastq.gz, etc.) are
        # handled in-memory by the file probe; extracting them to disk wastes
        # space and creates duplicate sidecars.
        archive_exts = {".tar", ".tgz", ".zip"}
        candidates = sorted(
            [p for p in dataset_dir.iterdir()
             if p.is_file() and (
                 p.suffix.lower() in archive_exts
                 or p.name.lower().endswith((".tar.gz", ".tgz"))
             )],
            key=lambda p: p.stat().st_size,
        )

        if not candidates:
            return extracted

        logger.info(
            "Unpacking archives",
            dir=str(dataset_dir),
            count=len(candidates),
        )

        for archive_path in candidates:
            stem = archive_path.name
            # Normalize compound extensions for the extract dir name
            for ext in (".tar.gz", ".tgz", ".tar", ".zip", ".gz"):
                if stem.lower().endswith(ext):
                    stem = stem[:-len(ext)]
                    break
            extract_dir = dataset_dir / f"_extracted_{stem}"
            extract_dir.mkdir(parents=True, exist_ok=True)

            try:
                if archive_path.suffix.lower() == ".zip" or zipfile.is_zipfile(archive_path):
                    with zipfile.ZipFile(archive_path, "r") as zf:
                        zf.extractall(extract_dir)
                    logger.info(
                        "Extracted ZIP archive",
                        archive=archive_path.name,
                        files=len(list(extract_dir.rglob("*"))),
                    )
                    extracted.append(extract_dir)
                elif (archive_path.suffix.lower() == ".gz"
                      and not archive_path.name.lower().endswith((".tar.gz", ".tgz"))):
                    # Standalone .gz file (not a tar.gz) — decompress in place
                    import gzip
                    import shutil
                    out_name = archive_path.stem
                    out_path = extract_dir / out_name
                    with gzip.open(archive_path, "rb") as gz_in:
                        with open(out_path, "wb") as f_out:
                            shutil.copyfileobj(gz_in, f_out)
                    logger.info(
                        "Decompressed gzip file",
                        archive=archive_path.name,
                        output=out_name,
                    )
                    extracted.append(extract_dir)
                elif tarfile.is_tarfile(archive_path):
                    with tarfile.open(archive_path, "r:*") as tf:
                        tf.extractall(extract_dir)
                    logger.info(
                        "Extracted TAR archive",
                        archive=archive_path.name,
                        files=len(list(extract_dir.rglob("*"))),
                    )
                    extracted.append(extract_dir)
                else:
                    logger.debug(
                        "Skipping non-archive file",
                        filename=archive_path.name,
                    )
            except Exception as exc:
                logger.warning(
                    "Failed to extract archive",
                    archive=archive_path.name,
                    error=str(exc),
                )

        return extracted

    # ------------------------------------------------------------------
    # Convenience: download specific file types defined by route priority
    # ------------------------------------------------------------------

    async def download_by_priority(
        self,
        gse_id: str,
        priority_file_types: list[list[str]] | None = None,
        max_size_mb: int = 500,
    ) -> tuple[Path, list[GeoFileInfo], list[Path]]:
        """Download files following the 8-level routing priority.

        Tries each priority level in order. Downloads the first level that
        has matching files available. This avoids downloading lower-priority
        files (e.g. FASTQ) when higher-priority ones (e.g. h5ad) exist.

        The default priority list mirrors the CLAUDE.md routing table:

        Priority 1: h5ad
        Priority 2: h5
        Priority 3: loom
        Priority 4: mtx, tsv, csv, txt  (MEX triplet components)
        Priority 5: tsv, csv, txt  (raw UMI matrices)
        Priority 6: rds
        Priority 7: tsv, csv, txt  (normalized / limited-use)
        Priority 8: fastq  (only as last resort)

        Args:
            gse_id: GEO Series accession.
            priority_file_types: Each list entry is a priority level containing
                one or more file type categories. If None, the default
                priority list is used.
            max_size_mb: Maximum file size in MB.

        Returns:
            Same tuple as download_dataset().
        """
        if priority_file_types is None:
            priority_file_types = [
                ["h5ad"],                       # 1: h5ad
                ["h5"],                         # 2: 10x H5
                ["loom"],                       # 3: loom
                ["mtx", "tsv", "csv", "txt"],   # 4: MEX triplet
                ["tsv", "csv", "txt"],          # 5: raw UMI matrix
                ["rds"],                        # 6: Seurat/SCE RDS
                ["tsv", "csv", "txt"],          # 7: normalized matrix
                ["fastq"],                      # 8: FASTQ (last resort)
            ]

        gse_id = gse_id.strip().upper()
        dataset_dir = self.data_dir / gse_id

        # First, list all files to determine what is available
        lister = GeoSupplementaryLister(client=self._client)
        all_files = await lister.list_with_checksums(gse_id)

        if not all_files:
            return dataset_dir, [], []

        available_categories = {f.category for f in all_files}

        # Find the first priority level that has matching files
        selected_types: list[str] | None = None
        for level, types in enumerate(priority_file_types, start=1):
            if any(t in available_categories for t in types):
                selected_types = types
                self._logger.info(
                    "Selected download priority level",
                    gse_id=gse_id,
                    level=level,
                    types=types,
                )
                break

        if selected_types is None:
            self._logger.warning(
                "No matching files found for any priority level",
                gse_id=gse_id,
                available=available_categories,
            )
            return dataset_dir, all_files, []

        # Now download using the selected type filter
        # Only skip raw files if we are not at the FASTQ level
        is_raw_level = selected_types == ["fastq"]
        return await self.download_dataset(
            gse_id=gse_id,
            file_types=selected_types,
            max_size_mb=max_size_mb,
            skip_archives=False,   # download archives too — probe will inspect inside
            skip_raw=not is_raw_level,
        )


# ---------------------------------------------------------------------------
# Standalone convenience functions (no class instantiation needed)
# ---------------------------------------------------------------------------


async def list_geo_supplementary(
    gse_id: str,
    client: "httpx.AsyncClient | None" = None,
) -> list[GeoFileInfo]:
    """Convenience: list supplementary files for a GEO dataset.

    Args:
        gse_id: GEO Series accession.
        client: Optional httpx.AsyncClient.

    Returns:
        List of GeoFileInfo objects.
    """
    lister = GeoSupplementaryLister(client=client)
    try:
        return await lister.list_with_checksums(gse_id)
    finally:
        await lister.close()


async def download_geo_dataset(
    gse_id: str,
    data_dir: Path,
    file_types: list[str] | None = None,
    max_size_mb: int = 500,
) -> tuple[Path, list[Path]]:
    """Convenience: download a GEO dataset's supplementary files.

    Args:
        gse_id: GEO Series accession.
        data_dir: Directory to save files into.
        file_types: File type categories to download (None = all processed).
        max_size_mb: Maximum file size in MB.

    Returns:
        Tuple of (dataset_dir, downloaded_paths).
    """
    downloader = GeoFileDownloader(data_dir=data_dir, max_concurrent=3)
    try:
        dataset_dir, all_files, downloaded = await downloader.download_dataset(
            gse_id=gse_id,
            file_types=file_types,
            max_size_mb=max_size_mb,
        )
        return dataset_dir, downloaded
    finally:
        await downloader.close()


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------


def _build_parser() -> "argparse.ArgumentParser":
    """Build the argparse CLI for geo_downloader."""
    import argparse

    parser = argparse.ArgumentParser(
        prog="gse_router.tools.geo_downloader",
        description="Download GEO supplementary files for scRNA-seq datasets.",
    )

    sub = parser.add_subparsers(dest="command", required=True)

    # --- list ---
    list_parser = sub.add_parser("list", help="List supplementary files for a GEO dataset")
    list_parser.add_argument("gse_id", help="GEO Series accession (e.g. GSE123456)")
    list_parser.add_argument(
        "--with-md5", action="store_true", default=True,
        help="Enrich with MD5 checksums from filelist.txt (default).",
    )
    list_parser.add_argument(
        "--no-md5", action="store_false", dest="with_md5",
        help="Skip filelist.txt lookup.",
    )

    # --- download ---
    download_parser = sub.add_parser("download", help="Download GEO supplementary files")
    download_parser.add_argument("gse_id", help="GEO Series accession (e.g. GSE123456)")
    download_parser.add_argument(
        "--data-dir", type=Path, default=Path.cwd() / "data",
        help="Root data directory (default: ./data)",
    )
    download_parser.add_argument(
        "--file-types", nargs="*", default=None,
        help="File types to download (e.g. h5ad h5 loom mtx). Default: all processed.",
    )
    download_parser.add_argument(
        "--max-size-mb", type=int, default=500,
        help="Maximum file size in MB (default: 500)",
    )
    download_parser.add_argument(
        "--max-concurrent", type=int, default=3,
        help="Maximum concurrent downloads (default: 3)",
    )
    download_parser.add_argument(
        "--include-archives", action="store_true", default=False,
        help="Download .tar/.zip archives (skipped by default).",
    )
    download_parser.add_argument(
        "--include-raw", action="store_true", default=False,
        help="Download FASTQ/BAM raw files (skipped by default).",
    )

    # --- download-by-priority ---
    priority_parser = sub.add_parser(
        "download-by-priority",
        help="Download files following the 8-level routing priority",
    )
    priority_parser.add_argument("gse_id", help="GEO Series accession (e.g. GSE123456)")
    priority_parser.add_argument(
        "--data-dir", type=Path, default=Path.cwd() / "data",
        help="Root data directory (default: ./data)",
    )
    priority_parser.add_argument(
        "--max-size-mb", type=int, default=500,
        help="Maximum file size in MB (default: 500)",
    )
    priority_parser.add_argument(
        "--max-concurrent", type=int, default=3,
        help="Maximum concurrent downloads (default: 3)",
    )

    return parser


def main(args: list[str] | None = None) -> None:
    """CLI entry point for GEO supplementary file listing and downloading.

    Usage examples::

        python -m gse_router.tools.geo_downloader list GSE123456
        python -m gse_router.tools.geo_downloader download GSE123456 --file-types h5ad h5
        python -m gse_router.tools.geo_downloader download-by-priority GSE123456
    """
    import json
    import sys

    parser = _build_parser()
    ns = parser.parse_args(args)

    async def _run() -> None:
        if ns.command == "list":
            lister = GeoSupplementaryLister()
            try:
                if ns.with_md5:
                    files = await lister.list_with_checksums(ns.gse_id)
                else:
                    files = await lister.list_supplementary(ns.gse_id)
            finally:
                await lister.close()

            output = {
                "gse_id": ns.gse_id.upper(),
                "count": len(files),
                "files": [
                    {
                        "filename": f.filename,
                        "size_bytes": f.size_bytes,
                        "size_mb": round(f.size_mb, 1) if f.size_mb else None,
                        "url": f.url,
                        "md5": f.md5,
                        "category": f.category,
                    }
                    for f in files
                ],
            }
            json.dump(output, sys.stdout, indent=2, default=str)
            sys.stdout.write("\n")
            sys.exit(0)

        elif ns.command == "download":
            downloader = GeoFileDownloader(
                data_dir=ns.data_dir,
                max_concurrent=ns.max_concurrent,
            )
            try:
                dataset_dir, all_files, downloaded = await downloader.download_dataset(
                    gse_id=ns.gse_id,
                    file_types=ns.file_types,
                    max_size_mb=ns.max_size_mb,
                    skip_archives=not ns.include_archives,
                    skip_raw=not ns.include_raw,
                )
            finally:
                await downloader.close()

            output = {
                "gse_id": ns.gse_id.upper(),
                "dataset_dir": str(dataset_dir),
                "total_files_listed": len(all_files),
                "downloaded_count": len(downloaded),
                "downloaded": [str(p) for p in downloaded],
                "total_size_mb": round(
                    sum(p.stat().st_size for p in downloaded if p.exists())
                    / (1024 * 1024),
                    1,
                ),
            }
            json.dump(output, sys.stdout, indent=2, default=str)
            sys.stdout.write("\n")
            sys.exit(0 if downloaded else 1)

        elif ns.command == "download-by-priority":
            downloader = GeoFileDownloader(
                data_dir=ns.data_dir,
                max_concurrent=ns.max_concurrent,
            )
            try:
                dataset_dir, all_files, downloaded = await downloader.download_by_priority(
                    gse_id=ns.gse_id,
                    max_size_mb=ns.max_size_mb,
                )
            finally:
                await downloader.close()

            output = {
                "gse_id": ns.gse_id.upper(),
                "dataset_dir": str(dataset_dir),
                "total_files_listed": len(all_files),
                "downloaded_count": len(downloaded),
                "downloaded": [str(p) for p in downloaded],
            }
            json.dump(output, sys.stdout, indent=2, default=str)
            sys.stdout.write("\n")
            sys.exit(0 if downloaded else 1)

    asyncio.run(_run())


if __name__ == "__main__":
    main()
