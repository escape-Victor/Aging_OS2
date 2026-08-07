"""Unified download backend — Aspera FASP, Aria2 multi-threaded, and httpx fallback.

This module provides a single ``DownloadAccelerator`` facade that selects the
best available download transport for each file based on URL host, file size,
and backend availability. All backends implement the same async interface:

1. **Aspera** (``ascp`` CLI) — UDP-based FASP protocol. Ideal for large files
   on NCBI FTP. Requires Aspera Connect to be installed.
2. **Aria2** (``aria2c`` CLI) — multi-threaded HTTP with chunked parallelism.
   Best for HTTP/HTTPS files > 50 MB. Requires aria2 to be installed.
3. **Httpx** (pure Python) — always-available fallback using
   ``httpx.AsyncClient`` streaming. Used for small files and when no
   accelerator CLI is installed.

Backend selection follows a deterministic priority chain:
   NCBI host + large file → Aspera → Aria2 → httpx
   Any HTTP host + large file → Aria2 → httpx
   Small file (<50 MB) → httpx (avoid subprocess overhead)

On failure, each backend falls through to the next available one. If all
backends fail for a file, the error is raised to the caller.
"""

from __future__ import annotations

import asyncio
import hashlib
import os
import re
import shutil
import sys
import time
from pathlib import Path
from typing import ClassVar

import structlog

from .download_config import DownloadBackend, DownloadConfig, load_download_config

logger = structlog.get_logger(__name__)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

# Files larger than this (bytes) are candidates for accelerator backends.
# Smaller files go straight to httpx — subprocess overhead isn't worth it.
_ACCELERATOR_MIN_SIZE: int = 50 * 1024 * 1024  # 50 MB

# NCBI hosts that support Aspera FASP transfers.
_NCBI_HOSTS: frozenset[str] = frozenset({
    "ftp.ncbi.nlm.nih.gov",
    "ftp.ncbi.nih.gov",
    "ftp-trace.ncbi.nlm.nih.gov",
})

# Default Aspera remote user (NCBI's public anonymous Aspera account).
_DEFAULT_ASPERA_USER: str = "anonftp"

# ---------------------------------------------------------------------------
# Aspera downloader
# ---------------------------------------------------------------------------


class AsperaDownloader:
    """Wraps the IBM Aspera ``ascp`` CLI for FASP-based high-speed transfers.

    FASP (Fast and Secure Protocol) uses UDP instead of TCP, making it
    immune to TCP's latency-based throughput limits. This is the fastest
    way to download large files from NCBI.

    Command pattern::

        ascp -v -i <key_path> -k 1 -T -l 300m
             anonftp@ftp.ncbi.nlm.nih.gov:/geo/series/GSE123nnn/GSE123456/suppl/file.tar
             /local/path/

    Requirements:
        - Aspera Connect installed (https://www.ibm.com/aspera/connect/)
        - ``ascp`` on PATH
        - Aspera public key (freely distributed by NCBI)
    """

    def __init__(self, config: DownloadConfig) -> None:
        self.config = config
        self._logger = logger.bind(component="AsperaDownloader")

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    async def download(
        self,
        remote_path: str,
        local_dir: Path,
        remote_host: str = "ftp.ncbi.nlm.nih.gov",
        remote_user: str = _DEFAULT_ASPERA_USER,
        timeout: float = 3600.0,
    ) -> Path:
        """Download a file via Aspera FASP.

        Args:
            remote_path: Absolute path on the remote Aspera server, e.g.
                ``"/geo/series/GSE123nnn/GSE123456/suppl/GSE123456_RAW.tar"``.
            local_dir: Directory to save the downloaded file. The filename
                is derived from the last component of ``remote_path``.
            remote_host: Aspera server hostname.
            remote_user: Aspera username (default: ``"anonftp"`` for NCBI).
            timeout: Maximum time in seconds.

        Returns:
            Path to the downloaded file.

        Raises:
            RuntimeError: If ascp exits with a non-zero code.
            FileNotFoundError: If ascp is not on PATH.
        """
        local_dir = Path(local_dir)
        local_dir.mkdir(parents=True, exist_ok=True)

        filename = remote_path.rstrip("/").split("/")[-1]
        dest = local_dir / filename

        cmd = self._build_command(remote_path, local_dir, remote_host, remote_user)
        self._logger.info(
            "Starting Aspera download",
            remote_path=remote_path,
            dest=str(local_dir),
            cmd=" ".join(cmd),
        )

        started = time.monotonic()
        await self._run_ascp(cmd, timeout=timeout)
        elapsed = time.monotonic() - started

        if not dest.exists():
            raise RuntimeError(
                f"Aspera download completed but file not found: {dest}"
            )

        size_mb = dest.stat().st_size / (1024 * 1024)
        speed_mbps = (dest.stat().st_size / (1024 * 1024) / elapsed) if elapsed > 0 else 0

        self._logger.info(
            "Aspera download complete",
            filename=filename,
            size_mb=round(size_mb, 1),
            elapsed_s=round(elapsed, 1),
            speed_mbps=round(speed_mbps, 1),
        )
        return dest

    # ------------------------------------------------------------------
    # Command building
    # ------------------------------------------------------------------

    def _build_command(
        self,
        remote_path: str,
        local_dir: Path,
        remote_host: str,
        remote_user: str,
    ) -> list[str]:
        """Build the ascp command line from configuration.

        Flags:
            -v     Verbose (progress goes to stderr)
            -i     Identity file (private key)
            -k 1   Resume partial transfers (checkpoint restart)
            -T     Disable in-transit encryption (faster on trusted networks)
            -l N   Bandwidth cap (e.g. "300m" = 300 Mbps)
        """
        cmd = ["ascp", "-v"]

        if self.config.aspera_key_path:
            cmd.extend(["-i", self.config.aspera_key_path])

        if self.config.aspera_resume:
            cmd.extend(["-k", "1"])

        if self.config.aspera_disable_encryption:
            cmd.append("-T")

        if self.config.aspera_bandwidth_limit:
            cmd.extend(["-l", self.config.aspera_bandwidth_limit])

        # Remote spec: user@host:path
        remote_spec = f"{remote_user}@{remote_host}:{remote_path}"
        cmd.append(remote_spec)
        cmd.append(str(local_dir))

        return cmd

    async def _run_ascp(self, cmd: list[str], timeout: float) -> None:
        """Execute ascp as a subprocess, streaming stderr for progress logging.

        ascp writes progress information to stderr. We capture and log it
        at DEBUG level for troubleshooting.

        Args:
            cmd: Full command line as a list of strings.
            timeout: Maximum runtime in seconds.

        Raises:
            RuntimeError: If ascp exits with a non-zero return code.
            asyncio.TimeoutError: If the process exceeds ``timeout``.
        """
        try:
            process = await asyncio.wait_for(
                asyncio.create_subprocess_exec(
                    *cmd,
                    stdout=asyncio.subprocess.PIPE,
                    stderr=asyncio.subprocess.PIPE,
                ),
                timeout=10.0,  # process start timeout
            )
        except FileNotFoundError:
            raise FileNotFoundError(
                "ascp not found. Install Aspera Connect: "
                "https://www.ibm.com/aspera/connect/"
            )

        # Stream stderr in the background while waiting for the process
        async def _stream_stderr() -> None:
            if process.stderr is None:
                return
            async for line in process.stderr:
                decoded = line.decode("utf-8", errors="replace").strip()
                if decoded:
                    self._logger.debug("ascp", line=decoded)

        stream_task = asyncio.create_task(_stream_stderr())

        try:
            await asyncio.wait_for(process.wait(), timeout=timeout)
        except asyncio.TimeoutError:
            process.kill()
            await process.wait()
            raise RuntimeError(
                f"Aspera download timed out after {timeout:.0f}s"
            )

        await stream_task

        if process.returncode != 0:
            raise RuntimeError(
                f"ascp exited with code {process.returncode}"
            )


# ---------------------------------------------------------------------------
# Aria2 downloader
# ---------------------------------------------------------------------------


class Aria2Downloader:
    """Wraps ``aria2c`` CLI for multi-threaded HTTP downloads.

    Aria2 splits a single file into multiple chunks and downloads them in
    parallel over multiple TCP connections, dramatically improving throughput
    compared to single-threaded HTTP clients like wget or httpx.

    Command pattern::

        aria2c -x 16 -s 16 -c --min-split-size=20M -d /local/dir "URL"

    Flags:
        -x N   Max connections per server (default 16)
        -s N   Number of file splits (default 16)
        -c     Continue/resume partial downloads
        -d DIR Destination directory

    Requirements:
        - aria2 installed (``aria2c`` on PATH)
        - https://aria2.github.io/
    """

    def __init__(self, config: DownloadConfig) -> None:
        self.config = config
        self._logger = logger.bind(component="Aria2Downloader")

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    async def download(
        self,
        url: str,
        local_dir: Path,
        output_filename: str | None = None,
        timeout: float = 3600.0,
    ) -> Path:
        """Download a file via aria2c.

        Args:
            url: Full HTTP/HTTPS URL to download.
            local_dir: Directory to save the file.
            output_filename: If provided, use this as the output filename.
                Otherwise, aria2c derives the name from the URL.
            timeout: Maximum time in seconds.

        Returns:
            Path to the downloaded file.

        Raises:
            RuntimeError: If aria2c exits with a non-zero code.
            FileNotFoundError: If aria2c is not on PATH.
        """
        local_dir = Path(local_dir)
        local_dir.mkdir(parents=True, exist_ok=True)

        cmd = self._build_command(url, local_dir, output_filename)
        self._logger.info(
            "Starting aria2 download",
            url=url,
            dest=str(local_dir),
            cmd=" ".join(cmd),
        )

        started = time.monotonic()
        await self._run_aria2(cmd, timeout=timeout)
        elapsed = time.monotonic() - started

        # Aria2 names the output file after the URL if -o is not given
        if output_filename:
            dest = local_dir / output_filename
        else:
            dest = local_dir / url.rstrip("/").split("/")[-1].split("?")[0]

        if not dest.exists():
            # Aria2 may have added a .aria2 control file suffix — try without it
            raise RuntimeError(
                f"aria2c download completed but file not found: {dest}"
            )

        size_mb = dest.stat().st_size / (1024 * 1024)
        speed_mbps = (dest.stat().st_size / (1024 * 1024) / elapsed) if elapsed > 0 else 0

        self._logger.info(
            "Aria2 download complete",
            filename=dest.name,
            size_mb=round(size_mb, 1),
            elapsed_s=round(elapsed, 1),
            speed_mbps=round(speed_mbps, 1),
        )
        return dest

    # ------------------------------------------------------------------
    # Command building
    # ------------------------------------------------------------------

    def _build_command(
        self,
        url: str,
        local_dir: Path,
        output_filename: str | None,
    ) -> list[str]:
        """Build the aria2c command line."""
        cmd = [
            "aria2c",
            "-x", str(self.config.aria2_max_connections),
            "-s", str(self.config.aria2_split),
            "--min-split-size", self.config.aria2_min_split_size,
            "-d", str(local_dir),
            "--summary-interval=30",  # Progress summary every 30s
        ]

        if self.config.aria2_continue:
            cmd.append("-c")

        if output_filename:
            cmd.extend(["-o", output_filename])

        # Suppress console output (we log via structlog)
        cmd.append("--console-log-level=warn")

        cmd.append(url)
        return cmd

    async def _run_aria2(self, cmd: list[str], timeout: float) -> None:
        """Execute aria2c as a subprocess.

        Args:
            cmd: Full command line.
            timeout: Maximum runtime in seconds.

        Raises:
            RuntimeError: On non-zero exit or timeout.
            FileNotFoundError: If aria2c is not on PATH.
        """
        try:
            process = await asyncio.wait_for(
                asyncio.create_subprocess_exec(
                    *cmd,
                    stdout=asyncio.subprocess.PIPE,
                    stderr=asyncio.subprocess.PIPE,
                ),
                timeout=10.0,
            )
        except FileNotFoundError:
            raise FileNotFoundError(
                "aria2c not found. Install aria2: https://aria2.github.io/"
            )

        async def _stream_output() -> None:
            # Aria2 writes progress to stderr by default
            if process.stderr is None:
                return
            async for line in process.stderr:
                decoded = line.decode("utf-8", errors="replace").strip()
                if decoded:
                    self._logger.debug("aria2c", line=decoded)

        stream_task = asyncio.create_task(_stream_output())

        try:
            await asyncio.wait_for(process.wait(), timeout=timeout)
        except asyncio.TimeoutError:
            process.kill()
            await process.wait()
            raise RuntimeError(
                f"Aria2 download timed out after {timeout:.0f}s"
            )

        await stream_task

        if process.returncode != 0:
            raise RuntimeError(
                f"aria2c exited with code {process.returncode}"
            )


# ---------------------------------------------------------------------------
# Httpx downloader (always-available fallback)
# ---------------------------------------------------------------------------


class HttpxDownloader:
    """Httpx-based fallback downloader. Always available — no CLI tools needed.

    Uses the same streaming + MD5 verification pattern as the existing
    ``GeoFileDownloader``. This is extracted into its own class so the
    accelerator can treat it as just another backend.
    """

    def __init__(self, config: DownloadConfig) -> None:
        self.config = config
        self._logger = logger.bind(component="HttpxDownloader")
        self._client: "httpx.AsyncClient | None" = None

    async def _get_client(self) -> "httpx.AsyncClient":
        if self._client is None:
            import httpx
            self._client = httpx.AsyncClient(
                timeout=self.config.httpx_timeout,
                follow_redirects=True,
                headers={
                    "User-Agent": (
                        "Mozilla/5.0 (compatible; gse-router/0.1; "
                        "+https://github.com/example/gse-router)"
                    ),
                },
            )
        return self._client

    async def close(self) -> None:
        if self._client is not None:
            await self._client.aclose()
            self._client = None

    async def download(
        self,
        url: str,
        dest: Path,
        expected_size: int | None = None,
        expected_md5: str | None = None,
        max_retries: int = 3,
    ) -> Path:
        """Download a file via httpx streaming.

        Args:
            url: Full URL.
            dest: Destination path.
            expected_size: Expected file size in bytes (for verification).
            expected_md5: Expected MD5 hex digest.
            max_retries: Maximum retry attempts.

        Returns:
            Path to the downloaded file.

        Raises:
            httpx.HTTPStatusError: On HTTP errors after retries exhausted.
            OSError: On filesystem errors.
            IOError: On MD5 mismatch.
        """
        import httpx

        dest.parent.mkdir(parents=True, exist_ok=True)
        filename = dest.name

        last_error: Exception | None = None
        for attempt in range(max_retries):
            if attempt > 0:
                delay = min(2 ** attempt, 30)
                self._logger.info(
                    "Retrying httpx download",
                    attempt=attempt + 1,
                    max_retries=max_retries,
                    delay_s=delay,
                    filename=filename,
                )
                await asyncio.sleep(delay)
                if dest.exists():
                    dest.unlink()

            self._logger.info("Downloading via httpx", url=url, dest=str(dest))
            start_time = time.monotonic()

            client = await self._get_client()
            try:
                async with client.stream("GET", url) as response:
                    response.raise_for_status()
                    total_size = int(response.headers.get("content-length", 0))
                    with open(dest, "wb") as fh:
                        async for chunk in response.aiter_bytes(chunk_size=65536):
                            fh.write(chunk)
                    actual_size = dest.stat().st_size
            except httpx.HTTPStatusError as exc:
                self._logger.error(
                    "Httpx download HTTP error",
                    url=url,
                    status=exc.response.status_code,
                )
                last_error = exc
                continue
            except Exception as exc:
                self._logger.error("Httpx download failed", url=url, error=str(exc))
                if dest.exists():
                    dest.unlink()
                last_error = exc
                continue

            elapsed = time.monotonic() - start_time
            speed_mbps = (actual_size / (1024 * 1024) / elapsed) if elapsed > 0 else 0

            # Size verification
            if expected_size is not None and actual_size != expected_size:
                self._logger.warning(
                    "Size mismatch", expected=expected_size, actual=actual_size,
                )
            elif total_size > 0 and actual_size != total_size:
                self._logger.warning(
                    "Incomplete download (content-length mismatch)",
                    expected=total_size, actual=actual_size,
                )

            # MD5 verification
            if expected_md5 and self.config.verify_checksums:
                md5_ok = await self._verify_md5(dest, expected_md5)
                if not md5_ok:
                    self._logger.error("MD5 mismatch", filename=filename)
                    dest.unlink()
                    raise IOError(
                        f"MD5 checksum mismatch for '{filename}': expected {expected_md5}"
                    )

            self._logger.info(
                "Httpx download complete",
                filename=filename,
                size_mb=round(actual_size / (1024 * 1024), 1),
                elapsed_s=round(elapsed, 1),
                speed_mbps=round(speed_mbps, 1),
            )
            return dest

        # All retries exhausted
        raise last_error  # type: ignore[misc]

    @staticmethod
    async def _verify_md5(file_path: Path, expected_md5: str) -> bool:
        """Verify a downloaded file's MD5 checksum."""
        if not expected_md5:
            return True
        hasher = hashlib.md5()
        with open(file_path, "rb") as fh:
            for chunk in iter(lambda: fh.read(8192), b""):
                hasher.update(chunk)
        return hasher.hexdigest() == expected_md5.lower()


# ---------------------------------------------------------------------------
# Download accelerator (main facade)
# ---------------------------------------------------------------------------


class DownloadAccelerator:
    """Unified download facade — selects the best backend per file.

    Typical usage::

        accel = DownloadAccelerator(data_dir=Path("./data"))
        dest = await accel.download_file(
            url="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE123nnn/GSE123456/suppl/matrix.mtx.gz",
            dest=Path("./data/GSE123456/matrix.mtx.gz"),
            expected_size=52428800,
        )

    For multiple files::

        files = [
            ("https://.../matrix.mtx.gz", Path("./data/GSE123456/matrix.mtx.gz")),
            ("https://.../barcodes.tsv.gz", Path("./data/GSE123456/barcodes.tsv.gz")),
        ]
        results = await accel.download_files(files)

    The accelerator auto-detects available backends on first use. Detection
    is cached so subsequent calls are instantaneous.
    """

    # URL pattern to detect NCBI GEO FTP URLs.
    _NCBI_URL_RE: ClassVar[re.Pattern] = re.compile(
        r"https?://[^/]*ncbi\.[^/]*nih\.gov/",
        re.IGNORECASE,
    )

    def __init__(
        self,
        config: DownloadConfig | None = None,
        data_dir: Path | None = None,
        max_concurrent: int = 3,
    ) -> None:
        """Initialise the accelerator.

        Args:
            config: Download configuration. If None, loaded from env vars.
            data_dir: Root directory for downloads.
            max_concurrent: Maximum simultaneous downloads (semaphore-limited).
        """
        self.config = config or load_download_config()
        self.data_dir = data_dir or Path(".")
        self.max_concurrent = max_concurrent
        self._semaphore: asyncio.Semaphore | None = None
        self._available: dict[DownloadBackend, bool] = {}
        self._detection_done: bool = False
        self._logger = logger.bind(component="DownloadAccelerator")

        # Backend instances (lazily initialised)
        self._aspera: AsperaDownloader | None = None
        self._aria2: Aria2Downloader | None = None
        self._httpx: HttpxDownloader | None = None

    # ------------------------------------------------------------------
    # Backend detection
    # ------------------------------------------------------------------

    async def detect_available_backends(self) -> dict[DownloadBackend, bool]:
        """Check which backends are installed and usable. Results are cached."""
        if self._detection_done:
            return dict(self._available)

        results: dict[DownloadBackend, bool] = {
            DownloadBackend.HTTPX: True,  # Always available
            DownloadBackend.ASPERA: await self._check_aspera_available(),
            DownloadBackend.ARIA2: await self._check_aria2_available(),
        }

        self._available = results
        self._detection_done = True

        self._logger.info(
            "Backend availability",
            **{k.value: v for k, v in results.items()},
        )
        return dict(results)

    async def _check_aspera_available(self) -> bool:
        """Check whether ascp is on PATH or at known locations."""
        # Standard PATH check
        if shutil.which("ascp"):
            # Also check that the key is available
            if self.config.aspera_key_path and Path(self.config.aspera_key_path).is_file():
                return True
            # If no key configured, still mark as available — the key check
            # will fail at download time with a clear error
            return True

        # Windows: check WSL
        if sys.platform == "win32":
            try:
                result = await asyncio.create_subprocess_exec(
                    "wsl", "which", "ascp",
                    stdout=asyncio.subprocess.PIPE,
                    stderr=asyncio.subprocess.PIPE,
                )
                stdout, _ = await asyncio.wait_for(result.communicate(), timeout=5.0)
                if result.returncode == 0 and stdout.strip():
                    return True
            except Exception:
                pass

        return False

    async def _check_aria2_available(self) -> bool:
        """Check whether aria2c is on PATH."""
        if shutil.which("aria2c"):
            return True

        # Windows: check WSL
        if sys.platform == "win32":
            try:
                result = await asyncio.create_subprocess_exec(
                    "wsl", "which", "aria2c",
                    stdout=asyncio.subprocess.PIPE,
                    stderr=asyncio.subprocess.PIPE,
                )
                stdout, _ = await asyncio.wait_for(result.communicate(), timeout=5.0)
                if result.returncode == 0 and stdout.strip():
                    return True
            except Exception:
                pass

        return False

    # ------------------------------------------------------------------
    # Smart backend selection
    # ------------------------------------------------------------------

    def select_backend(
        self,
        url: str,
        file_size: int | None = None,
    ) -> DownloadBackend:
        """Choose the best backend for a given URL and file characteristics.

        Decision matrix:

        ==================== ============ ======= ======= =======
        URL host              File size    Aspera  Aria2   Httpx
        ==================== ============ ======= ======= =======
        NCBI (ftp.ncbi.nlm)   > 50 MB      ✓✅      ✓       ✓
        NCBI (ftp.ncbi.nlm)   ≤ 50 MB      —       ✓✅      ✓
        Any HTTP              > 50 MB      —       ✓✅      ✓
        Any HTTP              ≤ 50 MB      —       —       ✓✅
        ==================== ============ ======= ======= =======

        ✓✅ = selected; ✓ = fallback; — = skipped

        Args:
            url: Full download URL.
            file_size: File size in bytes (from directory listing).

        Returns:
            The recommended backend enum value.
        """
        if self.config.backend != DownloadBackend.AUTO:
            return self.config.backend

        is_ncbi = self._is_ncbi_url(url)
        is_large = file_size is None or file_size >= _ACCELERATOR_MIN_SIZE

        available = self._available
        if not self._detection_done:
            # Synchronous fallback: check PATH directly (fast, no I/O)
            available = {
                DownloadBackend.HTTPX: True,
                DownloadBackend.ASPERA: shutil.which("ascp") is not None,
                DownloadBackend.ARIA2: shutil.which("aria2c") is not None,
            }

        # NCBI + large + Aspera available → Aspera
        if is_ncbi and is_large and available.get(DownloadBackend.ASPERA, False):
            return DownloadBackend.ASPERA

        # NCBI + small + Aria2 available → Aria2 (better than single-thread httpx)
        if is_ncbi and available.get(DownloadBackend.ARIA2, False):
            return DownloadBackend.ARIA2

        # Any HTTP + large + Aria2 available → Aria2
        if is_large and available.get(DownloadBackend.ARIA2, False):
            return DownloadBackend.ARIA2

        # Everything else → httpx
        return DownloadBackend.HTTPX

    @classmethod
    def _is_ncbi_url(cls, url: str) -> bool:
        """Check whether a URL points to an NCBI server."""
        return bool(cls._NCBI_URL_RE.search(url))

    # ------------------------------------------------------------------
    # Core download API
    # ------------------------------------------------------------------

    async def download_file(
        self,
        url: str,
        dest: Path,
        expected_size: int | None = None,
        expected_md5: str | None = None,
        backend: DownloadBackend | None = None,
        max_retries: int | None = None,
    ) -> Path:
        """Download a single file with automatic backend selection and fallback.

        Tries the preferred backend; if it fails, falls back through the
        chain (Aspera → Aria2 → httpx). On mirror-aware calls, also tries
        mirror URLs before giving up.

        Args:
            url: Full download URL.
            dest: Destination path on the local filesystem.
            expected_size: Expected file size in bytes.
            expected_md5: Expected MD5 hex digest.
            backend: Force a specific backend (overrides auto-selection).
            max_retries: Max attempts per backend (default from config).

        Returns:
            Path to the downloaded file.

        Raises:
            RuntimeError: If all backends fail.
        """
        # Ensure backends are detected
        if not self._detection_done:
            await self.detect_available_backends()

        max_retries = max_retries if max_retries is not None else self.config.max_retries
        selected = backend or self.select_backend(url, expected_size)

        # Build fallback chain
        backends_to_try: list[DownloadBackend] = [selected]
        if selected != DownloadBackend.ASPERA and self._available.get(DownloadBackend.ASPERA, False) and self._is_ncbi_url(url):
            backends_to_try.insert(0, DownloadBackend.ASPERA)
        if DownloadBackend.ARIA2 not in backends_to_try and self._available.get(DownloadBackend.ARIA2, False):
            backends_to_try.append(DownloadBackend.ARIA2)
        if DownloadBackend.HTTPX not in backends_to_try:
            backends_to_try.append(DownloadBackend.HTTPX)

        # Deduplicate while preserving order
        seen: set[DownloadBackend] = set()
        backends_to_try = [b for b in backends_to_try if not (b in seen or seen.add(b))]  # type: ignore[func-returns-value]

        self._logger.debug(
            "Download fallback chain",
            url=url,
            chain=[b.value for b in backends_to_try],
        )

        last_error: Exception | None = None
        for try_backend in backends_to_try:
            try:
                return await self._download_with_backend(
                    url=url,
                    dest=dest,
                    backend=try_backend,
                    expected_size=expected_size,
                    expected_md5=expected_md5,
                    max_retries=max_retries,
                )
            except Exception as exc:
                self._logger.warning(
                    "Backend failed, trying fallback",
                    backend=try_backend.value,
                    url=url,
                    error=str(exc),
                )
                last_error = exc
                if dest.exists():
                    dest.unlink(missing_ok=True)

        raise RuntimeError(
            f"All download backends failed for '{url}'. "
            f"Last error: {last_error}"
        )

    async def _download_with_backend(
        self,
        url: str,
        dest: Path,
        backend: DownloadBackend,
        expected_size: int | None,
        expected_md5: str | None,
        max_retries: int,
    ) -> Path:
        """Dispatch a download to the appropriate backend implementation."""
        if backend == DownloadBackend.ASPERA:
            return await self._download_via_aspera(url, dest)
        elif backend == DownloadBackend.ARIA2:
            return await self._download_via_aria2(url, dest)
        else:
            return await self._download_via_httpx(
                url, dest, expected_size, expected_md5, max_retries,
            )

    async def _download_via_aspera(self, url: str, dest: Path) -> Path:
        """Convert an HTTP URL to an Aspera remote path and download."""
        if self._aspera is None:
            self._aspera = AsperaDownloader(self.config)

        # Convert https://ftp.ncbi.nlm.nih.gov/geo/series/GSE123nnn/GSE123456/suppl/file.tar
        # → /geo/series/GSE123nnn/GSE123456/suppl/file.tar
        remote_path = self._url_to_aspera_path(url)
        return await self._aspera.download(
            remote_path=remote_path,
            local_dir=dest.parent,
        )

    async def _download_via_aria2(self, url: str, dest: Path) -> Path:
        """Download via aria2c."""
        if self._aria2 is None:
            self._aria2 = Aria2Downloader(self.config)

        return await self._aria2.download(
            url=url,
            local_dir=dest.parent,
            output_filename=dest.name,
        )

    async def _download_via_httpx(
        self,
        url: str,
        dest: Path,
        expected_size: int | None,
        expected_md5: str | None,
        max_retries: int,
    ) -> Path:
        """Download via httpx streaming."""
        if self._httpx is None:
            self._httpx = HttpxDownloader(self.config)

        return await self._httpx.download(
            url=url,
            dest=dest,
            expected_size=expected_size,
            expected_md5=expected_md5,
            max_retries=max_retries,
        )

    @staticmethod
    def _url_to_aspera_path(url: str) -> str:
        """Convert an HTTPS NCBI URL to an Aspera remote path.

        https://ftp.ncbi.nlm.nih.gov/geo/series/GSE123nnn/GSE123456/suppl/file.tar
        → /geo/series/GSE123nnn/GSE123456/suppl/file.tar
        """
        # Strip protocol and host
        path = re.sub(r"^https?://[^/]+", "", url)
        # Ensure leading slash
        if not path.startswith("/"):
            path = "/" + path
        return path

    # ------------------------------------------------------------------
    # Multiple file downloads
    # ------------------------------------------------------------------

    async def download_files(
        self,
        files: list[tuple[str, Path]],
        max_concurrent: int | None = None,
    ) -> list[Path]:
        """Download multiple files with semaphore-limited concurrency.

        Args:
            files: List of ``(url, dest)`` pairs.
            max_concurrent: Override the instance-level concurrency cap.

        Returns:
            List of local Paths for successfully downloaded files.
            Failed downloads are logged but do not halt the batch.
        """
        cap = max_concurrent or self.max_concurrent
        self._logger.info("Starting batch accelerator download", count=len(files), max_concurrent=cap)

        tasks: list[asyncio.Task] = []
        for url, dest in files:
            task = asyncio.create_task(
                self._download_one_safe(url, dest),
                name=dest.name,
            )
            tasks.append(task)

        outcomes = await asyncio.gather(*tasks, return_exceptions=True)
        results: list[Path] = []
        for i, outcome in enumerate(outcomes):
            if isinstance(outcome, Exception):
                self._logger.error(
                    "Batch download task failed",
                    url=files[i][0],
                    error=str(outcome),
                )
            elif outcome is not None:
                results.append(outcome)

        self._logger.info(
            "Batch accelerator download finished",
            requested=len(files),
            succeeded=len(results),
            failed=len(files) - len(results),
        )
        return results

    async def _download_one_safe(self, url: str, dest: Path) -> Path | None:
        """Download a single file, returning None on failure."""
        try:
            return await self.download_file(url=url, dest=dest)
        except Exception as exc:
            self._logger.error(
                "Download failed",
                url=url,
                filename=dest.name,
                error=str(exc),
            )
            return None

    # ------------------------------------------------------------------
    # Cleanup
    # ------------------------------------------------------------------

    async def close(self) -> None:
        """Close backend clients."""
        if self._httpx is not None:
            await self._httpx.close()
            self._httpx = None
        self._detection_done = False
        self._available.clear()
