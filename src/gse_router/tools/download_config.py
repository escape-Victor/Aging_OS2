"""Centralised download configuration — env vars, YAML overrides, and auto-detection.

This module has zero dependencies on other gse_router modules. It can be
imported standalone and is used by both the download accelerator and the
MCP tool wrappers.

Configuration precedence (highest to lowest):
1. Environment variables (e.g. ``ASPERA_KEY_PATH``)
2. YAML config passed via ``load_download_config(yaml_config=...)``
3. Sensible defaults suitable for most environments
"""

from __future__ import annotations

import os
import sys
from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path


# ---------------------------------------------------------------------------
# Backend enum (shared across all download modules)
# ---------------------------------------------------------------------------


class DownloadBackend(str, Enum):
    """Available download transport backends.

    ``AUTO`` (the default) auto-selects the best available backend per file
    based on URL host and file size.
    """

    ASPERA = "aspera"   # ascp CLI — FASP protocol (UDP-based, high-speed)
    ARIA2 = "aria2"     # aria2c CLI — multi-threaded HTTP with resume
    HTTPX = "httpx"     # httpx.AsyncClient — always available pure-Python fallback
    AUTO = "auto"       # smart selection: Aspera > Aria2 > httpx


# ---------------------------------------------------------------------------
# Configuration dataclass
# ---------------------------------------------------------------------------


@dataclass
class DownloadConfig:
    """Download accelerator configuration.

    All fields have sensible defaults. Use ``load_download_config()`` to
    populate from environment variables and optional YAML overrides.

    Attributes:
        backend: Default download backend (``"auto"`` for smart selection).
        aspera_key_path: Path to Aspera's ``asperaweb_id_dsa.openssh`` private key.
            Auto-detected from common install locations if left empty.
        aspera_bandwidth_limit: Bandwidth cap for ascp (``-l`` flag).
        aspera_disable_encryption: Disable in-transit encryption for speed (``-T``).
        aspera_resume: Enable checkpoint resume (``-k 1``).
        aria2_max_connections: Max connections per server (``-x``).
        aria2_split: Number of file splits for parallel download (``-s``).
        aria2_continue: Resume partial downloads (``-c``).
        aria2_min_split_size: Minimum chunk size for splitting.
        mirror_priority: Ordered list of mirror sources to try.
        mirror_health_check_timeout: Seconds for mirror HEAD requests.
        httpx_timeout: Timeout for httpx-based downloads.
        httpx_listing_timeout: Timeout for directory listing requests.
        http_proxy: HTTP proxy URL.
        https_proxy: HTTPS proxy URL.
        max_retries: Retry count per file before falling back.
        verify_checksums: Verify MD5 after download when available.
    """

    # --- Backend selection ---
    backend: DownloadBackend = DownloadBackend.AUTO

    # --- Aspera / FASP ---
    aspera_key_path: str = ""
    aspera_bandwidth_limit: str = "300m"
    aspera_disable_encryption: bool = True
    aspera_resume: bool = True

    # --- Aria2 ---
    aria2_max_connections: int = 16
    aria2_split: int = 16
    aria2_continue: bool = True
    aria2_min_split_size: str = "20M"

    # --- Mirrors ---
    mirror_priority: list[str] = field(default_factory=lambda: [
        "ncbi",
        "ebi_ae",
        "ena",
        "cncb",
        "ddbj",
    ])
    mirror_health_check_timeout: float = 10.0

    # --- HTTP / httpx ---
    httpx_timeout: float = 600.0
    httpx_listing_timeout: float = 30.0
    http_proxy: str = ""
    https_proxy: str = ""

    # --- General ---
    max_retries: int = 3
    verify_checksums: bool = True
    max_download_size_mb: float = 500.0


# ---------------------------------------------------------------------------
# Aspera key auto-detection
# ---------------------------------------------------------------------------

# Known install locations for ``asperaweb_id_dsa.openssh`` across platforms.
# NCBI distributes this public-private keypair freely — it is NOT a secret.
_ASPERA_KEY_CANDIDATES: list[tuple[str, ...]] = [
    # Linux / macOS — Aspera Connect from IBM website
    (".aspera", "connect", "etc", "asperaweb_id_dsa.openssh"),
    (".aspera", "connect", "etc", "asperaweb_id_dsa.putty"),
    # macOS — Homebrew
    ("opt", "aspera-connect", "etc", "asperaweb_id_dsa.openssh"),
    # Linux — system-wide install
    ("opt", "aspera", "etc", "asperaweb_id_dsa.openssh"),
    ("usr", "local", "aspera", "etc", "asperaweb_id_dsa.openssh"),
    # Conda / micromamba
    ("etc", "asperaweb_id_dsa.openssh"),
]


def _windows_aspera_paths() -> list[Path]:
    """Build candidate Aspera key paths for Windows."""
    candidates: list[Path] = []
    local_appdata = os.environ.get("LOCALAPPDATA", "")
    program_files = os.environ.get("PROGRAMFILES", "")
    program_files_x86 = os.environ.get("ProgramFiles(x86)", "")

    for base in (local_appdata, program_files, program_files_x86):
        if not base:
            continue
        base_path = Path(base)
        candidates.extend([
            base_path / "Aspera" / "Aspera Connect" / "etc" / "asperaweb_id_dsa.openssh",
            base_path / "Aspera" / "Aspera Connect" / "etc" / "asperaweb_id_dsa.putty",
            base_path / "IBM" / "Aspera Connect" / "etc" / "asperaweb_id_dsa.openssh",
        ])
    return candidates


def _default_aspera_key() -> str:
    """Auto-detect the Aspera private key across common install locations.

    Checks, in order:
    1. ``ASPERA_KEY_PATH`` environment variable (explicit override)
    2. Known platform-specific install directories
    3. ``shutil.which("ascp")`` sibling path
    4. Returns ``""`` if no key found (Aspera backend will be disabled)

    Returns:
        Absolute path to the key file as a string, or ``""`` if not found.
    """
    # 1. Explicit env override
    env_key = os.environ.get("ASPERA_KEY_PATH", "")
    if env_key:
        env_path = Path(env_key)
        if env_path.is_file():
            return str(env_path)

    # 2. Known install locations
    if sys.platform == "win32":
        for candidate in _windows_aspera_paths():
            if candidate.is_file():
                return str(candidate)
    else:
        home = Path.home()
        for parts in _ASPERA_KEY_CANDIDATES:
            candidate = home.joinpath(*parts)
            if candidate.is_file():
                return str(candidate)

    # 3. Check alongside ascp binary
    import shutil
    ascp_path = shutil.which("ascp")
    if ascp_path:
        ascp_dir = Path(ascp_path).parent.parent / "etc"
        for key_name in ("asperaweb_id_dsa.openssh", "asperaweb_id_dsa.putty"):
            key_path = ascp_dir / key_name
            if key_path.is_file():
                return str(key_path)

    return ""


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def load_download_config(
    yaml_config: dict | None = None,
) -> DownloadConfig:
    """Build ``DownloadConfig`` from environment variables and optional YAML.

    Environment variables (highest precedence):

    .. list-table::
       :header-rows: 1

       * - Variable
         - Type
         - Default
         - Description
       * - ``GSE_ROUTER_DOWNLOAD_BACKEND``
         - str
         - ``"auto"``
         - ``aspera`` | ``aria2`` | ``httpx`` | ``auto``
       * - ``ASPERA_KEY_PATH``
         - str
         - auto-detect
         - Path to asperaweb_id_dsa.openssh
       * - ``ASPERA_BANDWIDTH_LIMIT``
         - str
         - ``"300m"``
         - Bandwidth cap (``-l`` flag)
       * - ``ARIA2_MAX_CONNECTIONS``
         - int
         - ``16``
         - Max connections per server
       * - ``ARIA2_SPLIT``
         - int
         - ``16``
         - File split count
       * - ``GSE_ROUTER_MIRRORS``
         - str
         - ``"ncbi,ebi_ae,ena,cncb,ddbj"``
         - Comma-separated mirror priority
       * - ``GSE_ROUTER_HTTP_PROXY``
         - str
         - ``""``
         - HTTP proxy URL
       * - ``GSE_ROUTER_HTTPS_PROXY``
         - str
         - ``""``
         - HTTPS proxy URL
       * - ``GSE_ROUTER_DOWNLOAD_TIMEOUT``
         - int
         - ``600``
         - Download timeout (seconds)
       * - ``GSE_ROUTER_LISTING_TIMEOUT``
         - int
         - ``30``
         - Directory listing timeout (seconds)

    Args:
        yaml_config: Optional dict from ``config/tools.yaml`` or similar.
            Keys should match the field names in ``DownloadConfig``.

    Returns:
        A fully populated ``DownloadConfig`` instance.
    """
    # Start with defaults
    config = DownloadConfig()

    # --- YAML overrides (medium precedence) ---
    if yaml_config:
        _apply_yaml_overrides(config, yaml_config)

    # --- Environment variables (highest precedence) ---
    _apply_env_overrides(config)

    # --- Auto-detect Aspera key if not already set ---
    if not config.aspera_key_path:
        config.aspera_key_path = _default_aspera_key()

    return config


def _apply_yaml_overrides(config: DownloadConfig, yaml: dict) -> None:
    """Apply YAML config values to a DownloadConfig instance.

    Only sets fields that exist on DownloadConfig and are present in yaml.
    """
    # Top-level download section
    dl_section = yaml.get("download", yaml)

    str_fields = (
        "backend", "aspera_key_path", "aspera_bandwidth_limit",
        "http_proxy", "https_proxy",
    )
    for field_name in str_fields:
        if field_name in dl_section:
            setattr(config, field_name, dl_section[field_name])

    int_fields = (
        "aria2_max_connections", "aria2_split", "max_retries",
    )
    for field_name in int_fields:
        if field_name in dl_section:
            try:
                setattr(config, field_name, int(dl_section[field_name]))
            except (TypeError, ValueError):
                pass

    float_fields = (
        "httpx_timeout", "httpx_listing_timeout", "mirror_health_check_timeout",
        "max_download_size_mb",
    )
    for field_name in float_fields:
        if field_name in dl_section:
            try:
                setattr(config, field_name, float(dl_section[field_name]))
            except (TypeError, ValueError):
                pass

    bool_fields = ("aspera_disable_encryption", "aspera_resume", "aria2_continue", "verify_checksums")
    for field_name in bool_fields:
        if field_name in dl_section:
            try:
                setattr(config, field_name, bool(dl_section[field_name]))
            except (TypeError, ValueError):
                pass

    # Mirror priority
    if "mirrors" in dl_section:
        mirrors = dl_section["mirrors"]
        if isinstance(mirrors, list):
            config.mirror_priority = mirrors

    # Nested sections (aspera:, aria2:)
    for section, prefix, field_map in [
        ("aspera", "", {
            "key_path": "aspera_key_path",
            "bandwidth_limit": "aspera_bandwidth_limit",
            "disable_encryption": "aspera_disable_encryption",
            "resume": "aspera_resume",
        }),
        ("aria2", "", {
            "max_connections": "aria2_max_connections",
            "split": "aria2_split",
            "continue": "aria2_continue",
            "min_split_size": "aria2_min_split_size",
        }),
        ("timeouts", "", {
            "listing": "httpx_listing_timeout",
            "download": "httpx_timeout",
            "mirror_health_check": "mirror_health_check_timeout",
        }),
    ]:
        nested = dl_section.get(section, {})
        if isinstance(nested, dict):
            for yaml_key, config_key in field_map.items():
                if yaml_key in nested:
                    val = nested[yaml_key]
                    if isinstance(val, str) and config_key in int_fields:
                        try:
                            val = int(val)
                        except (TypeError, ValueError):
                            pass
                    setattr(config, config_key, val)


def _apply_env_overrides(config: DownloadConfig) -> None:
    """Apply environment variable overrides to a DownloadConfig instance."""
    _env_str(config, "GSE_ROUTER_DOWNLOAD_BACKEND", "backend")
    _env_str(config, "ASPERA_KEY_PATH", "aspera_key_path")
    _env_str(config, "ASPERA_BANDWIDTH_LIMIT", "aspera_bandwidth_limit")
    _env_str(config, "GSE_ROUTER_HTTP_PROXY", "http_proxy")
    _env_str(config, "GSE_ROUTER_HTTPS_PROXY", "https_proxy")
    _env_str(config, "ARIA2_MIN_SPLIT_SIZE", "aria2_min_split_size")

    _env_int(config, "ARIA2_MAX_CONNECTIONS", "aria2_max_connections")
    _env_int(config, "ARIA2_SPLIT", "aria2_split")
    _env_int(config, "GSE_ROUTER_DOWNLOAD_TIMEOUT", "httpx_timeout")
    _env_int(config, "GSE_ROUTER_LISTING_TIMEOUT", "httpx_listing_timeout")
    _env_int(config, "GSE_ROUTER_MAX_RETRIES", "max_retries")
    _env_float(config, "GSE_ROUTER_MAX_DOWNLOAD_SIZE_MB", "max_download_size_mb")

    _env_bool(config, "ASPERA_DISABLE_ENCRYPTION", "aspera_disable_encryption")
    _env_bool(config, "ASPERA_RESUME", "aspera_resume")
    _env_bool(config, "ARIA2_CONTINUE", "aria2_continue")
    _env_bool(config, "GSE_ROUTER_VERIFY_CHECKSUMS", "verify_checksums")

    # Mirror priority from env
    mirrors_env = os.environ.get("GSE_ROUTER_MIRRORS", "")
    if mirrors_env:
        config.mirror_priority = [m.strip() for m in mirrors_env.split(",") if m.strip()]


def _env_str(config: DownloadConfig, env_name: str, attr: str) -> None:
    val = os.environ.get(env_name, "")
    if val:
        setattr(config, attr, val)


def _env_int(config: DownloadConfig, env_name: str, attr: str) -> None:
    val = os.environ.get(env_name, "")
    if val:
        try:
            setattr(config, attr, int(val))
        except (TypeError, ValueError):
            pass


def _env_float(config: DownloadConfig, env_name: str, attr: str) -> None:
    val = os.environ.get(env_name, "")
    if val:
        try:
            setattr(config, attr, float(val))
        except (TypeError, ValueError):
            pass


def _env_bool(config: DownloadConfig, env_name: str, attr: str) -> None:
    val = os.environ.get(env_name, "")
    if val:
        setattr(config, attr, val.lower() in ("1", "true", "yes", "on"))
