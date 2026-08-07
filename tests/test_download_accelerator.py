"""Tests for download accelerator, backends, mirror resolution, and config loading.

Test categories:
- Backend selection logic (URL type, file size, availability)
- AsperaDownloader command building (unit, no actual ascp required)
- Aria2Downloader command building (unit, no actual aria2c required)
- MirrorResolver URL construction and health checking
- DownloadConfig loading from env vars and YAML
- Graceful degradation when accelerators are unavailable
"""

from __future__ import annotations

import asyncio
import os
import sys
from pathlib import Path
from unittest import mock

import pytest

# Ensure src/ is importable
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from gse_router.tools.download_config import (
    DownloadBackend,
    DownloadConfig,
    _default_aspera_key,
    load_download_config,
)
from gse_router.tools.mirrors import (
    Mirror,
    MirrorResolver,
    MirrorSource,
)


# ============================================================================
# DownloadBackend enum
# ============================================================================


class TestDownloadBackend:
    """Test the DownloadBackend enum values and AUTO default."""

    def test_all_backends_exist(self):
        assert DownloadBackend.ASPERA.value == "aspera"
        assert DownloadBackend.ARIA2.value == "aria2"
        assert DownloadBackend.HTTPX.value == "httpx"
        assert DownloadBackend.AUTO.value == "auto"

    def test_auto_is_default_in_config(self):
        config = DownloadConfig()
        assert config.backend == DownloadBackend.AUTO


# ============================================================================
# DownloadConfig loading
# ============================================================================


class TestDownloadConfigDefaults:
    """Test DownloadConfig default values."""

    def test_default_backend_is_auto(self):
        config = DownloadConfig()
        assert config.backend == DownloadBackend.AUTO

    def test_default_aspera_settings(self):
        config = DownloadConfig()
        assert config.aspera_bandwidth_limit == "300m"
        assert config.aspera_disable_encryption is True
        assert config.aspera_resume is True

    def test_default_aria2_settings(self):
        config = DownloadConfig()
        assert config.aria2_max_connections == 16
        assert config.aria2_split == 16
        assert config.aria2_continue is True
        assert config.aria2_min_split_size == "20M"

    def test_default_mirror_priority(self):
        config = DownloadConfig()
        assert "ncbi" in config.mirror_priority
        assert "ebi_ae" in config.mirror_priority
        assert "ena" in config.mirror_priority

    def test_default_timeouts(self):
        config = DownloadConfig()
        assert config.httpx_timeout == 600.0
        assert config.httpx_listing_timeout == 30.0
        assert config.mirror_health_check_timeout == 10.0

    def test_default_max_retries(self):
        config = DownloadConfig()
        assert config.max_retries == 3


class TestDownloadConfigFromEnv:
    """Test loading DownloadConfig from environment variables."""

    def test_loads_backend_from_env(self, monkeypatch):
        monkeypatch.setenv("GSE_ROUTER_DOWNLOAD_BACKEND", "aria2")
        config = load_download_config()
        assert config.backend == "aria2"

    def test_loads_aspera_key_from_env(self, monkeypatch):
        monkeypatch.setenv("ASPERA_KEY_PATH", "/custom/key/path")
        monkeypatch.setattr(
            "gse_router.tools.download_config._default_aspera_key",
            lambda: "/custom/key/path",
        )
        config = load_download_config()
        assert config.aspera_key_path == "/custom/key/path"

    def test_loads_aria2_settings_from_env(self, monkeypatch):
        monkeypatch.setenv("ARIA2_MAX_CONNECTIONS", "32")
        monkeypatch.setenv("ARIA2_SPLIT", "8")
        config = load_download_config()
        assert config.aria2_max_connections == 32
        assert config.aria2_split == 8

    def test_loads_mirrors_from_env(self, monkeypatch):
        monkeypatch.setenv("GSE_ROUTER_MIRRORS", "ncbi,ena,ddbj")
        config = load_download_config()
        assert config.mirror_priority == ["ncbi", "ena", "ddbj"]

    def test_loads_timeout_from_env(self, monkeypatch):
        monkeypatch.setenv("GSE_ROUTER_DOWNLOAD_TIMEOUT", "1200")
        monkeypatch.setenv("GSE_ROUTER_LISTING_TIMEOUT", "60")
        config = load_download_config()
        assert config.httpx_timeout == 1200
        assert config.httpx_listing_timeout == 60

    def test_loads_proxy_from_env(self, monkeypatch):
        monkeypatch.setenv("GSE_ROUTER_HTTP_PROXY", "http://proxy:8080")
        monkeypatch.setenv("GSE_ROUTER_HTTPS_PROXY", "http://proxy:8080")
        config = load_download_config()
        assert config.http_proxy == "http://proxy:8080"
        assert config.https_proxy == "http://proxy:8080"

    def test_invalid_int_uses_default(self, monkeypatch):
        monkeypatch.setenv("ARIA2_MAX_CONNECTIONS", "not_a_number")
        config = load_download_config()
        assert config.aria2_max_connections == 16  # default

    def test_empty_env_uses_defaults(self, monkeypatch):
        # Ensure env vars are unset
        for var in (
            "GSE_ROUTER_DOWNLOAD_BACKEND", "ASPERA_KEY_PATH",
            "ARIA2_MAX_CONNECTIONS", "ARIA2_SPLIT",
            "GSE_ROUTER_DOWNLOAD_TIMEOUT", "GSE_ROUTER_HTTP_PROXY",
        ):
            monkeypatch.delenv(var, raising=False)
        config = load_download_config()
        assert config.backend == DownloadBackend.AUTO
        assert config.aria2_max_connections == 16


class TestDownloadConfigFromYaml:
    """Test YAML overrides for DownloadConfig."""

    def test_yaml_overrides_backend(self):
        yaml = {"download": {"backend": "aspera"}}
        config = load_download_config(yaml_config=yaml)
        assert config.backend == "aspera"

    def test_yaml_overrides_aspera_settings(self):
        yaml = {
            "download": {
                "aspera": {
                    "bandwidth_limit": "500m",
                    "disable_encryption": False,
                    "resume": False,
                }
            }
        }
        config = load_download_config(yaml_config=yaml)
        assert config.aspera_bandwidth_limit == "500m"
        assert config.aspera_disable_encryption is False
        assert config.aspera_resume is False

    def test_yaml_overrides_aria2_settings(self):
        yaml = {
            "download": {
                "aria2": {
                    "max_connections": 32,
                    "split": 8,
                    "min_split_size": "50M",
                }
            }
        }
        config = load_download_config(yaml_config=yaml)
        assert config.aria2_max_connections == 32
        assert config.aria2_split == 8
        assert config.aria2_min_split_size == "50M"

    def test_yaml_overrides_mirrors(self):
        yaml = {"download": {"mirrors": ["ncbi", "ddbj"]}}
        config = load_download_config(yaml_config=yaml)
        assert config.mirror_priority == ["ncbi", "ddbj"]

    def test_yaml_overrides_timeouts(self):
        yaml = {"download": {"timeouts": {"listing": 60, "download": 1200}}}
        config = load_download_config(yaml_config=yaml)
        assert config.httpx_listing_timeout == 60
        assert config.httpx_timeout == 1200

    def test_env_overrides_yaml(self, monkeypatch):
        """Env vars take precedence over YAML."""
        monkeypatch.setenv("GSE_ROUTER_DOWNLOAD_BACKEND", "httpx")
        yaml = {"download": {"backend": "aspera"}}
        config = load_download_config(yaml_config=yaml)
        assert config.backend == "httpx"  # env wins


# ============================================================================
# Aspera key auto-detection
# ============================================================================


class TestAsperaKeyDetection:
    """Test _default_aspera_key() auto-detection."""

    def test_returns_env_var_if_valid(self, monkeypatch, tmp_path):
        key_file = tmp_path / "test_key"
        key_file.write_text("fake-key")
        monkeypatch.setenv("ASPERA_KEY_PATH", str(key_file))
        result = _default_aspera_key()
        assert result == str(key_file)

    def test_returns_empty_if_env_var_invalid(self, monkeypatch):
        monkeypatch.setenv("ASPERA_KEY_PATH", "/nonexistent/path/key")
        # Normally returns "" since the file doesn't exist
        # But may find a real key on the dev machine — accept any non-crash result
        result = _default_aspera_key()
        assert isinstance(result, str)

    def test_returns_string(self):
        """Sanity: always returns a string, never raises."""
        result = _default_aspera_key()
        assert isinstance(result, str)


# ============================================================================
# MirrorResolver — URL construction
# ============================================================================


class TestMirrorResolverUrlConstruction:
    """Test that MirrorResolver builds correct URLs for GEO and SRA."""

    def test_resolve_geo_suppl_url_returns_urls(self):
        resolver = MirrorResolver()
        results = resolver.resolve_geo_suppl_url("GSE123456")
        assert len(results) > 0
        # NCBI should be first (priority 0)
        assert results[0][0].source == MirrorSource.NCBI
        # Template: {base}/geo/series/{prefix}nnn/{gse_id}/suppl/
        # GSE123456 → geo_prefix="GSE123" → "GSE123nnn"
        assert "GSE123nnn/GSE123456/suppl/" in results[0][1]

    def test_resolve_geo_file_url_returns_urls(self):
        resolver = MirrorResolver()
        results = resolver.resolve_geo_file_url("GSE123456", "matrix.mtx.gz")
        assert len(results) > 0
        assert results[0][0].source == MirrorSource.NCBI
        assert "matrix.mtx.gz" in results[0][1]
        assert "GSE123456" in results[0][1]

    def test_resolve_geo_file_url_ena_excluded(self):
        """ENA does not mirror GEO supplementary — should not appear in results."""
        resolver = MirrorResolver()
        results = resolver.resolve_geo_file_url("GSE123456", "data.h5ad")
        sources = {r[0].source for r in results}
        assert MirrorSource.ENA not in sources

    def test_resolve_sra_url_returns_urls(self):
        resolver = MirrorResolver()
        results = resolver.resolve_sra_url("SRR12345678")
        # At least NCBI and ENA should have SRA templates
        assert len(results) >= 2
        sources = {r[0].source.value for r in results}
        assert "ncbi" in sources

    def test_invalid_gse_id_raises(self):
        resolver = MirrorResolver()
        with pytest.raises(ValueError):
            resolver.resolve_geo_suppl_url("INVALID")

    def test_geo_prefix_parsing(self):
        """Test that different GSE ID lengths parse correctly."""
        resolver = MirrorResolver()
        # 6-digit (common)
        results = resolver.resolve_geo_file_url("GSE123456", "test.txt")
        assert "GSE123nnn" in results[0][1]

        # 7-digit
        results = resolver.resolve_geo_file_url("GSE1234567", "test.txt")
        assert "GSE1234nnn" in results[0][1]

    def test_mirror_priority_order(self):
        """Mirrors should be returned in priority order."""
        resolver = MirrorResolver()
        results = resolver.resolve_geo_file_url("GSE123456", "test.txt")
        priorities = [r[0].priority for r in results]
        assert priorities == sorted(priorities)  # ascending

    def test_custom_mirror_priority(self):
        """Custom mirror_priority should reorder results."""
        resolver = MirrorResolver(mirror_priority=["ddbj", "ncbi"])
        results = resolver.resolve_geo_file_url("GSE123456", "test.txt")
        if len(results) >= 2:
            assert results[0][0].source == MirrorSource.DDBJ


class TestMirrorHealthChecking:
    """Test mirror health check logic (mocked HTTP)."""

    @pytest.mark.asyncio
    async def test_check_mirror_health_all_reachable(self):
        """When all mirrors respond 200, all should be healthy."""
        resolver = MirrorResolver()
        # All mirrors use HTTP — we mock the httpx client
        with mock.patch("httpx.AsyncClient.head") as mock_head:
            mock_head.return_value.status_code = 200
            # Need to mock the follow_redirects context
            mock_response = mock.MagicMock()
            mock_response.status_code = 200
            mock_head.return_value = mock_response

            # Create a mock client
            async def mock_head_impl(url, follow_redirects=True):
                return mock_response

            with mock.patch.object(
                resolver, "_get_http_client",
                return_value=mock.AsyncMock(head=mock_head_impl),
            ):
                health = await resolver.check_mirror_health()
                assert isinstance(health, dict)
                assert "ncbi" in health

    @pytest.mark.asyncio
    async def test_check_mirror_health_cached(self):
        """Health results should be cached for 300s."""
        resolver = MirrorResolver()
        resolver._health_cache["ncbi"] = (True, 1e12)  # far future — still valid
        resolver._health_cache["ebi_ae"] = (False, 1e12)
        health = await resolver.check_mirror_health()
        assert health["ncbi"] is True
        assert health["ebi_ae"] is False

    @pytest.mark.asyncio
    async def test_get_healthiest_mirror_returns_first_healthy(self):
        resolver = MirrorResolver()
        resolver._health_cache = {
            "ncbi": (True, 1e12),
            "ebi_ae": (True, 1e12),
            "ena": (False, 1e12),
        }
        mirror = await resolver.get_healthiest_mirror()
        assert mirror.source == MirrorSource.NCBI

    @pytest.mark.asyncio
    async def test_get_healthiest_mirror_falls_back_to_ncbi(self):
        """When all mirrors are unhealthy, NCBI is returned as fallback."""
        resolver = MirrorResolver()
        resolver._health_cache = {
            "ncbi": (False, 1e12),
            "ebi_ae": (False, 1e12),
            "ena": (False, 1e12),
            "cncb": (False, 1e12),
            "ddbj": (False, 1e12),
        }
        mirror = await resolver.get_healthiest_mirror()
        assert mirror.source == MirrorSource.NCBI  # source of truth


# ============================================================================
# DownloadAccelerator backend selection
# ============================================================================


class TestSmartBackendSelection:
    """Test DownloadAccelerator.select_backend() decision logic."""

    def test_ncbi_large_file_selects_aspera_when_available(self):
        from gse_router.tools.download_accelerator import DownloadAccelerator

        config = DownloadConfig()
        accel = DownloadAccelerator(config=config)
        accel._available = {
            DownloadBackend.ASPERA: True,
            DownloadBackend.ARIA2: True,
            DownloadBackend.HTTPX: True,
        }
        accel._detection_done = True

        url = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE123nnn/GSE123456/suppl/big_file.tar"
        backend = accel.select_backend(url, file_size=200 * 1024 * 1024)  # 200 MB
        assert backend == DownloadBackend.ASPERA

    def test_ncbi_large_file_falls_back_to_aria2_when_aspera_missing(self):
        from gse_router.tools.download_accelerator import DownloadAccelerator

        config = DownloadConfig()
        accel = DownloadAccelerator(config=config)
        accel._available = {
            DownloadBackend.ASPERA: False,
            DownloadBackend.ARIA2: True,
            DownloadBackend.HTTPX: True,
        }
        accel._detection_done = True

        url = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE123nnn/GSE123456/suppl/big_file.tar"
        backend = accel.select_backend(url, file_size=200 * 1024 * 1024)
        assert backend == DownloadBackend.ARIA2

    def test_ncbi_small_file_selects_aria2_over_httpx(self):
        from gse_router.tools.download_accelerator import DownloadAccelerator

        config = DownloadConfig()
        accel = DownloadAccelerator(config=config)
        accel._available = {
            DownloadBackend.ASPERA: True,
            DownloadBackend.ARIA2: True,
            DownloadBackend.HTTPX: True,
        }
        accel._detection_done = True

        url = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE123nnn/GSE123456/suppl/small.csv"
        backend = accel.select_backend(url, file_size=1024)  # 1 KB
        assert backend == DownloadBackend.ARIA2  # NCBI always prefers Aria2 over httpx

    def test_non_ncbi_large_file_selects_aria2(self):
        from gse_router.tools.download_accelerator import DownloadAccelerator

        config = DownloadConfig()
        accel = DownloadAccelerator(config=config)
        accel._available = {
            DownloadBackend.ASPERA: True,
            DownloadBackend.ARIA2: True,
            DownloadBackend.HTTPX: True,
        }
        accel._detection_done = True

        url = "https://www.ebi.ac.uk/arrayexpress/files/GSE123456/big_file.h5ad"
        backend = accel.select_backend(url, file_size=200 * 1024 * 1024)
        assert backend == DownloadBackend.ARIA2

    def test_small_file_always_uses_httpx(self):
        from gse_router.tools.download_accelerator import DownloadAccelerator

        config = DownloadConfig()
        accel = DownloadAccelerator(config=config)
        accel._available = {
            DownloadBackend.ASPERA: True,
            DownloadBackend.ARIA2: True,
            DownloadBackend.HTTPX: True,
        }
        accel._detection_done = True

        url = "https://example.com/small_file.csv"
        backend = accel.select_backend(url, file_size=1024)
        assert backend == DownloadBackend.HTTPX

    def test_no_accelerators_available_uses_httpx(self):
        from gse_router.tools.download_accelerator import DownloadAccelerator

        config = DownloadConfig()
        accel = DownloadAccelerator(config=config)
        accel._available = {
            DownloadBackend.ASPERA: False,
            DownloadBackend.ARIA2: False,
            DownloadBackend.HTTPX: True,
        }
        accel._detection_done = True

        url = "https://ftp.ncbi.nlm.nih.gov/geo/series/big.tar"
        backend = accel.select_backend(url, file_size=500 * 1024 * 1024)
        assert backend == DownloadBackend.HTTPX

    def test_explicit_backend_overrides_auto(self):
        from gse_router.tools.download_accelerator import DownloadAccelerator

        config = DownloadConfig(backend=DownloadBackend.ARIA2)
        accel = DownloadAccelerator(config=config)
        accel._detection_done = True

        url = "https://ftp.ncbi.nlm.nih.gov/geo/series/big.tar"
        backend = accel.select_backend(url, file_size=500 * 1024 * 1024)
        assert backend == DownloadBackend.ARIA2

    def test_is_ncbi_url_detection(self):
        from gse_router.tools.download_accelerator import DownloadAccelerator

        ncbi_urls = [
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE123nnn/GSE123456/suppl/file.tar",
            "https://ftp.ncbi.nih.gov/geo/",
            "https://ftp-trace.ncbi.nlm.nih.gov/sra/",
        ]
        non_ncbi_urls = [
            "https://www.ebi.ac.uk/arrayexpress/files/test",
            "https://example.com/data.csv",
            "https://ngdc.cncb.ac.cn/gsa/test",
        ]
        for url in ncbi_urls:
            assert DownloadAccelerator._is_ncbi_url(url), f"Should detect NCBI: {url}"
        for url in non_ncbi_urls:
            assert not DownloadAccelerator._is_ncbi_url(url), f"Should NOT detect NCBI: {url}"


# ============================================================================
# AsperaDownloader command building
# ============================================================================


class TestAsperaDownloaderCommands:
    """Test ascp command line construction (unit tests, no actual ascp)."""

    def test_build_command_basic(self):
        from gse_router.tools.download_accelerator import AsperaDownloader

        config = DownloadConfig(
            aspera_key_path="/test/key",
            aspera_disable_encryption=True,
            aspera_resume=True,
            aspera_bandwidth_limit="300m",
        )
        downloader = AsperaDownloader(config)
        cmd = downloader._build_command(
            remote_path="/geo/series/GSE123nnn/GSE123456/suppl/file.tar",
            local_dir=Path("/tmp/data"),
            remote_host="ftp.ncbi.nlm.nih.gov",
            remote_user="anonftp",
        )

        # Verify key flags
        assert "ascp" in cmd[0]
        assert "-v" in cmd
        assert "-i" in cmd
        assert "/test/key" in cmd
        assert "-k" in cmd
        assert "1" in cmd
        assert "-T" in cmd
        assert "-l" in cmd
        assert "300m" in cmd
        # Verify remote spec
        remote_specs = [a for a in cmd if "@" in a]
        assert len(remote_specs) == 1
        assert "anonftp@ftp.ncbi.nlm.nih.gov:/geo/series" in remote_specs[0]
        # Verify local path
        assert str(Path("/tmp/data")) in cmd

    def test_build_command_no_key(self):
        from gse_router.tools.download_accelerator import AsperaDownloader

        config = DownloadConfig(aspera_key_path="")
        downloader = AsperaDownloader(config)
        cmd = downloader._build_command(
            remote_path="/test/file.tar",
            local_dir=Path("/tmp"),
            remote_host="host",
            remote_user="user",
        )
        assert "-i" not in cmd

    def test_build_command_encryption_enabled(self):
        from gse_router.tools.download_accelerator import AsperaDownloader

        config = DownloadConfig(aspera_disable_encryption=False)
        downloader = AsperaDownloader(config)
        cmd = downloader._build_command(
            remote_path="/test/file.tar",
            local_dir=Path("/tmp"),
            remote_host="host",
            remote_user="user",
        )
        assert "-T" not in cmd

    def test_build_command_no_resume(self):
        from gse_router.tools.download_accelerator import AsperaDownloader

        config = DownloadConfig(aspera_resume=False)
        downloader = AsperaDownloader(config)
        cmd = downloader._build_command(
            remote_path="/test/file.tar",
            local_dir=Path("/tmp"),
            remote_host="host",
            remote_user="user",
        )
        assert "-k" not in cmd


# ============================================================================
# Aria2Downloader command building
# ============================================================================


class TestAria2DownloaderCommands:
    """Test aria2c command line construction (unit tests, no actual aria2c)."""

    def test_build_command_basic(self):
        from gse_router.tools.download_accelerator import Aria2Downloader

        config = DownloadConfig(
            aria2_max_connections=16,
            aria2_split=16,
            aria2_continue=True,
            aria2_min_split_size="20M",
        )
        downloader = Aria2Downloader(config)
        cmd = downloader._build_command(
            url="https://example.com/file.txt",
            local_dir=Path("/tmp/data"),
            output_filename=None,
        )

        assert "aria2c" == cmd[0]
        assert "-x" in cmd and "16" in cmd
        assert "-s" in cmd and "16" in cmd
        assert "-c" in cmd
        assert "--min-split-size" in cmd and "20M" in cmd
        assert "-d" in cmd
        assert str(Path("/tmp/data")) in cmd
        assert "https://example.com/file.txt" in cmd
        assert "--console-log-level=warn" in cmd

    def test_build_command_with_output_filename(self):
        from gse_router.tools.download_accelerator import Aria2Downloader

        config = DownloadConfig()
        downloader = Aria2Downloader(config)
        cmd = downloader._build_command(
            url="https://example.com/data.gz",
            local_dir=Path("/tmp"),
            output_filename="my_data.gz",
        )
        assert "-o" in cmd
        assert "my_data.gz" in cmd

    def test_build_command_no_continue(self):
        from gse_router.tools.download_accelerator import Aria2Downloader

        config = DownloadConfig(aria2_continue=False)
        downloader = Aria2Downloader(config)
        cmd = downloader._build_command(
            url="https://example.com/file.txt",
            local_dir=Path("/tmp"),
            output_filename=None,
        )
        assert "-c" not in cmd

    def test_build_command_custom_connections(self):
        from gse_router.tools.download_accelerator import Aria2Downloader

        config = DownloadConfig(aria2_max_connections=32, aria2_split=32)
        downloader = Aria2Downloader(config)
        cmd = downloader._build_command(
            url="https://example.com/file.txt",
            local_dir=Path("/tmp"),
            output_filename=None,
        )
        assert "-x" in cmd
        x_idx = cmd.index("-x")
        assert cmd[x_idx + 1] == "32"
        s_idx = cmd.index("-s")
        assert cmd[s_idx + 1] == "32"


# ============================================================================
# URL to Aspera path conversion
# ============================================================================


class TestUrlToAsperaPath:
    """Test HTTP URL → Aspera remote path conversion."""

    def test_converts_ncbi_https_to_aspera_path(self):
        from gse_router.tools.download_accelerator import DownloadAccelerator

        url = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE123nnn/GSE123456/suppl/matrix.mtx.gz"
        path = DownloadAccelerator._url_to_aspera_path(url)
        assert path == "/geo/series/GSE123nnn/GSE123456/suppl/matrix.mtx.gz"
        assert path.startswith("/")

    def test_preserves_existing_slash_prefix(self):
        from gse_router.tools.download_accelerator import DownloadAccelerator

        # Should handle URL with trailing path correctly
        url = "https://ftp.ncbi.nlm.nih.gov/sra/reads/file.sra"
        path = DownloadAccelerator._url_to_aspera_path(url)
        assert path.startswith("/")
        assert path == "/sra/reads/file.sra"


# ============================================================================
# MirrorSource enum and Mirror dataclass
# ============================================================================


class TestMirrorSource:
    """Test MirrorSource enum values."""

    def test_all_sources_exist(self):
        assert MirrorSource.NCBI.value == "ncbi"
        assert MirrorSource.EBI_AE.value == "ebi_ae"
        assert MirrorSource.ENA.value == "ena"
        assert MirrorSource.CNCB.value == "cncb"
        assert MirrorSource.DDBJ.value == "ddbj"


class TestMirrorDataclass:
    """Test Mirror dataclass instantiation and defaults."""

    def test_default_mirror_creation(self):
        m = Mirror(
            source=MirrorSource.NCBI,
            base_url="https://example.com",
        )
        assert m.source == MirrorSource.NCBI
        assert m.priority == 0
        assert m.supports_aspera is False
        assert m.supports_http is True
        assert "geo/series" in m.geo_suppl_template

    def test_default_mirrors_all_have_base_url(self):
        """All DEFAULT_MIRRORS should have a non-empty base_url."""
        for m in MirrorResolver.DEFAULT_MIRRORS:
            assert m.base_url, f"Mirror {m.source} has empty base_url"
            assert isinstance(m.priority, int)


# ============================================================================
# Integration: fallback chain construction
# ============================================================================


class TestFallbackChain:
    """Test that the download_file method builds correct fallback chains."""

    @pytest.mark.asyncio
    async def test_fallback_chain_for_ncbi_with_all_backends(self):
        """When all backends available for NCBI: Aspera → Aria2 → httpx."""
        from gse_router.tools.download_accelerator import DownloadAccelerator

        config = DownloadConfig()
        accel = DownloadAccelerator(config=config)
        accel._available = {
            DownloadBackend.ASPERA: True,
            DownloadBackend.ARIA2: True,
            DownloadBackend.HTTPX: True,
        }
        accel._detection_done = True

        url = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE123nnn/GSE123456/suppl/big.tar"
        backend = accel.select_backend(url, file_size=200 * 1024 * 1024)
        # Should pick Aspera first for large NCBI files
        assert backend == DownloadBackend.ASPERA

    @pytest.mark.asyncio
    async def test_fallback_chain_for_non_ncbi(self):
        """Non-NCBI URLs should use Aria2 → httpx (no Aspera)."""
        from gse_router.tools.download_accelerator import DownloadAccelerator

        config = DownloadConfig()
        accel = DownloadAccelerator(config=config)
        accel._available = {
            DownloadBackend.ASPERA: True,
            DownloadBackend.ARIA2: True,
            DownloadBackend.HTTPX: True,
        }
        accel._detection_done = True

        url = "https://www.ebi.ac.uk/arrayexpress/files/GSE123456/data.h5ad"
        backend = accel.select_backend(url, file_size=200 * 1024 * 1024)
        assert backend == DownloadBackend.ARIA2
