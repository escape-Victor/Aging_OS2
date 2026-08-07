"""Configuration loader — reads settings from environment variables and .env file.

Provides a single source of truth for all configuration. Settings are loaded
from, in priority order:
1. Environment variables (highest)
2. .env file in the project root
3. Built-in defaults (lowest)

Usage:
    from gse_router.config import get_config
    cfg = get_config()
    print(cfg.ncbi_api_key)
"""

from __future__ import annotations

import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import ClassVar


def _find_project_root() -> Path:
    """Find the project root directory (where .env or pyproject.toml lives)."""
    # Start from this file's location
    current = Path(__file__).resolve().parent.parent.parent
    # Check for pyproject.toml
    markers = ["pyproject.toml", ".env", ".git"]
    for marker in markers:
        if (current / marker).exists():
            return current
    # Fall back to current working directory
    return Path.cwd()


def _load_dotenv(env_path: Path | None = None) -> None:
    """Load environment variables from a .env file.

    Does NOT overwrite existing environment variables (os.environ takes priority).
    Handles comments, blank lines, quoted values, and inline comments.

    Args:
        env_path: Path to .env file. If None, auto-discovers from project root.
    """
    if env_path is None:
        env_path = _find_project_root() / ".env"

    if not env_path.is_file():
        return

    with open(env_path, encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            # Skip comments and blank lines
            if not line or line.startswith("#"):
                continue

            # Skip if no equals sign
            if "=" not in line:
                continue

            key, _, value = line.partition("=")
            key = key.strip()
            value = value.strip()

            # Remove surrounding quotes
            if len(value) >= 2 and value[0] == value[-1] and value[0] in ('"', "'"):
                value = value[1:-1]

            # Remove inline comments (but not inside quotes)
            if "#" in value and not (value.startswith('"') and value.endswith('"')):
                value = value.split("#")[0].strip()

            # Only set if not already in environment
            if key and key not in os.environ:
                os.environ[key] = value


@dataclass
class GseRouterConfig:
    """Configuration for GSE Router.

    All values can be set via environment variables or .env file.
    Use get_config() to obtain the singleton instance.
    """

    # --- NCBI / GEO API ---
    ncbi_api_key: str = ""
    ncbi_email: str = ""

    # --- Directories ---
    data_dir: Path = field(default_factory=lambda: Path("./data"))
    output_dir: Path = field(default_factory=lambda: Path("./output"))

    # --- Pipeline ---
    auto_download: bool = False
    max_download_size_mb: int = 500
    skip_raw_fastq: bool = True

    # --- Logging ---
    log_level: str = "INFO"

    # --- Singleton ---
    _instance: ClassVar[GseRouterConfig | None] = None

    @classmethod
    def from_env(cls) -> GseRouterConfig:
        """Create a config instance from environment variables."""
        return cls(
            ncbi_api_key=os.environ.get("NCBI_API_KEY", ""),
            ncbi_email=os.environ.get("NCBI_EMAIL", ""),
            data_dir=Path(os.environ.get("GSE_ROUTER_DATA_DIR", "./data")),
            output_dir=Path(os.environ.get("GSE_ROUTER_OUTPUT_DIR", "./output")),
            auto_download=os.environ.get("GSE_ROUTER_AUTO_DOWNLOAD", "").lower() in ("true", "1", "yes"),
            max_download_size_mb=int(os.environ.get("GSE_ROUTER_MAX_DOWNLOAD_SIZE_MB", "500")),
            skip_raw_fastq=os.environ.get("GSE_ROUTER_SKIP_RAW_FASTQ", "true").lower() in ("true", "1", "yes"),
            log_level=os.environ.get("GSE_ROUTER_LOG_LEVEL", "INFO"),
        )


def get_config(reload: bool = False) -> GseRouterConfig:
    """Get the global configuration singleton.

    Loads .env on first call. Subsequent calls return the cached instance.

    Args:
        reload: If True, reload from environment (ignoring cache).

    Returns:
        GseRouterConfig with all settings populated.
    """
    if GseRouterConfig._instance is None or reload:
        _load_dotenv()
        GseRouterConfig._instance = GseRouterConfig.from_env()
    return GseRouterConfig._instance


# ---------------------------------------------------------------------------
# Convenience: re-exported for CLI and tools
# ---------------------------------------------------------------------------


def get_ncbi_api_key() -> str:
    """Get the NCBI API key, checking multiple sources.

    Priority:
    1. NCBI_API_KEY env var / .env file
    2. NCBI_API_KEY environment variable set by parent process
    """
    return get_config().ncbi_api_key


def get_ncbi_email() -> str:
    """Get the NCBI contact email, checking multiple sources.

    Priority:
    1. NCBI_EMAIL env var / .env file
    2. --email CLI flag (must be passed separately)
    """
    return get_config().ncbi_email
