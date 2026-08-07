"""Mirror URL resolution — GEO/SRA data replicated across global bioinformatics centres.

NCBI data is synchronised to multiple international mirrors. Using a
geographically closer or better-connected mirror can dramatically improve
download speed and reliability compared to hitting ``ftp.ncbi.nlm.nih.gov``
from every location.

This module provides:
- URL construction for 5 major mirrors (NCBI, EBI ArrayExpress, ENA, CNCB, DDBJ)
- Health-checked mirror selection with configurable priority
- Automatic failover when a mirror is unreachable

References:
    EBI ArrayExpress: https://www.ebi.ac.uk/arrayexpress/
    ENA (European Nucleotide Archive): https://www.ebi.ac.uk/ena/
    CNCB / GSA (China National Bioinformation Center): https://ngdc.cncb.ac.cn/gsa/
    DDBJ / DRA (DNA Data Bank of Japan): https://www.ddbj.nig.ac.jp/
"""

from __future__ import annotations

import asyncio
import re
import time
from dataclasses import dataclass, field
from enum import Enum
from typing import ClassVar

import structlog

logger = structlog.get_logger(__name__)


# ---------------------------------------------------------------------------
# Enums
# ---------------------------------------------------------------------------


class MirrorSource(str, Enum):
    """Known database mirror sources for GEO and SRA data."""

    NCBI = "ncbi"        # ftp.ncbi.nlm.nih.gov — primary, supports Aspera
    EBI_AE = "ebi_ae"    # EBI ArrayExpress — mirrors GEO expression data
    ENA = "ena"          # EBI European Nucleotide Archive — mirrors SRA fully
    CNCB = "cncb"        # CNCB / GSA (China) — comprehensive NCBI mirror
    DDBJ = "ddbj"        # DDBJ / DRA (Japan) — mirrors SRA, Aspera-capable


# ---------------------------------------------------------------------------
# Mirror definition
# ---------------------------------------------------------------------------


@dataclass
class Mirror:
    """Definition of a single database mirror.

    Attributes:
        source: Which mirror this is (NCBI, EBI, etc.).
        base_url: Root URL for the mirror's HTTP/FTP interface.
        priority: Lower = preferred (0 = highest priority).
        supports_aspera: Whether this mirror accepts Aspera FASP transfers.
        supports_http: Whether this mirror serves files over HTTP/HTTPS.
        geo_suppl_template: URL template for GEO supplementary directory.
            Variables: ``{base}``, ``{prefix}``, ``{gse_id}``.
        geo_file_template: URL template for a single GEO supplementary file.
            Variables: ``{base}``, ``{prefix}``, ``{gse_id}``, ``{filename}``.
        sra_template: URL template for SRA run data.
    """

    source: MirrorSource
    base_url: str
    priority: int = 0
    supports_aspera: bool = False
    supports_http: bool = True
    geo_suppl_template: str = (
        "{base}/geo/series/{prefix}nnn/{gse_id}/suppl/"
    )
    geo_file_template: str = (
        "{base}/geo/series/{prefix}nnn/{gse_id}/suppl/{filename}"
    )
    sra_template: str = ""
    description: str = ""


# ---------------------------------------------------------------------------
# Mirror resolver
# ---------------------------------------------------------------------------


class MirrorResolver:
    """Resolves GEO and SRA URLs across multiple database mirrors.

    Constructs download URLs for each known mirror and can verify which
    mirrors are currently reachable via lightweight HEAD requests.

    Typical usage::

        resolver = MirrorResolver()
        urls = resolver.resolve_geo_file_url("GSE123456", "GSM3722100_matrix.mtx.gz")
        # urls[0] → NCBI URL, urls[1] → EBI ArrayExpress URL, etc.

        health = await resolver.check_mirror_health(gse_id="GSE123456")
        # health → {"ncbi": True, "ebi_ae": True, "ena": False, ...}

    Health results are cached for 300 seconds to avoid repeated checks
    during a single dataset processing session.
    """

    # Default mirrors in priority order.
    # NCBI is priority 0 because it is the authoritative source and the only
    # one that reliably supports Aspera FASP transfers.
    DEFAULT_MIRRORS: ClassVar[list[Mirror]] = [
        Mirror(
            source=MirrorSource.NCBI,
            base_url="https://ftp.ncbi.nlm.nih.gov",
            priority=0,
            supports_aspera=True,
            supports_http=True,
            geo_suppl_template="{base}/geo/series/{prefix}nnn/{gse_id}/suppl/",
            geo_file_template="{base}/geo/series/{prefix}nnn/{gse_id}/suppl/{filename}",
            sra_template="{base}/sra/sra-instant/reads/ByRun/sra/{srr_prefix}/{srr_prefix2}/{srr}/{srr}.sra",
            description="NCBI GEO FTP (primary, supports Aspera FASP)",
        ),
        Mirror(
            source=MirrorSource.EBI_AE,
            base_url="https://www.ebi.ac.uk/arrayexpress/files",
            priority=1,
            supports_aspera=False,
            supports_http=True,
            geo_suppl_template="{base}/{gse_id}/",
            geo_file_template="{base}/{gse_id}/{filename}",
            description="EBI ArrayExpress (mirrors GEO expression matrices)",
        ),
        Mirror(
            source=MirrorSource.ENA,
            base_url="https://www.ebi.ac.uk/ena",
            priority=2,
            supports_aspera=False,
            supports_http=True,
            geo_suppl_template="",   # ENA does not mirror GEO supplementary directly
            geo_file_template="",    # — use EBI_AE for GEO, ENA for SRA
            sra_template="{base}/browser/api/links/run?accession={srr}&format=fastq&protocol=ftp",
            description="EBI European Nucleotide Archive (full SRA mirror, FTP downloads)",
        ),
        Mirror(
            source=MirrorSource.CNCB,
            base_url="https://ngdc.cncb.ac.cn",
            priority=3,
            supports_aspera=False,
            supports_http=True,
            geo_suppl_template="{base}/gsa/browse/{gse_id}/",
            geo_file_template="{base}/gsa/browse/{gse_id}/suppl/{filename}",
            description="CNCB / GSA — China National Bioinformation Center",
        ),
        Mirror(
            source=MirrorSource.DDBJ,
            base_url="https://ddbj.nig.ac.jp/public",
            priority=4,
            supports_aspera=True,
            supports_http=True,
            geo_suppl_template="{base}/geo/series/{prefix}nnn/{gse_id}/suppl/",
            geo_file_template="{base}/geo/series/{prefix}nnn/{gse_id}/suppl/{filename}",
            sra_template="{base}/dra/{srr_prefix}/{srr}",
            description="DDBJ / DRA (Japan) — supports Aspera + high-speed FTP",
        ),
    ]

    # NCBI GEO organises series into 3-level directory hierarchy.
    # GSE123456 → prefix "GSE123", directory "GSE123nnn".
    _GEO_ID_RE: ClassVar[re.Pattern] = re.compile(r"^(GS[EPM])(\d+)$", re.IGNORECASE)

    # SRR accessions: first 6 chars = prefix (SRR123 → SRR123456)
    _SRR_PREFIX_RE: ClassVar[re.Pattern] = re.compile(r"^(SRR\d{3})")

    # ------------------------------------------------------------------
    # Construction
    # ------------------------------------------------------------------

    def __init__(
        self,
        mirrors: list[Mirror] | None = None,
        mirror_priority: list[str] | None = None,
    ) -> None:
        """Initialise the resolver.

        Args:
            mirrors: Custom mirror list. If None, uses ``DEFAULT_MIRRORS``.
            mirror_priority: Ordered source names (e.g. ``["ncbi", "ebi_ae"]``).
                If None, mirrors are sorted by their ``priority`` field.
        """
        self._mirrors = mirrors or list(self.DEFAULT_MIRRORS)
        if mirror_priority:
            # Reorder to match the requested priority
            priority_map = {name: idx for idx, name in enumerate(mirror_priority)}
            self._mirrors.sort(key=lambda m: priority_map.get(m.source.value, 999))
        else:
            self._mirrors.sort(key=lambda m: m.priority)
        self._health_cache: dict[str, tuple[bool, float]] = {}
        self._health_cache_ttl: float = 300.0  # 5 minutes
        self._logger = logger.bind(component="MirrorResolver")

    # ------------------------------------------------------------------
    # GEO URL resolution
    # ------------------------------------------------------------------

    @staticmethod
    def _parse_geo_id(gse_id: str) -> tuple[str, str]:
        """Extract the GSE prefix and directory prefix from a GEO accession.

        Args:
            gse_id: e.g. ``"GSE123456"``.

        Returns:
            Tuple of ``(geo_prefix, dir_prefix)``, e.g. ``("GSE123", "GSE123nnn")``.

        Raises:
            ValueError: If ``gse_id`` is not a valid GEO accession.
        """
        gse_id = gse_id.strip().upper()
        match = MirrorResolver._GEO_ID_RE.match(gse_id)
        if not match:
            raise ValueError(
                f"Invalid GEO accession '{gse_id}'. "
                f"Expected pattern: GSE/GPL/GSM followed by digits."
            )
        prefix_type = match.group(1)  # GSE, GPL, or GSM
        num_str = match.group(2)
        num = int(num_str)
        geo_prefix = f"{prefix_type}{num_str[:-3]}"  # GSE123
        dir_prefix = f"{geo_prefix}nnn"              # GSE123nnn
        return geo_prefix, dir_prefix

    def resolve_geo_suppl_url(self, gse_id: str) -> list[tuple[Mirror, str]]:
        """Return ordered list of ``(mirror, url)`` for a GEO supplementary directory.

        Only mirrors that have a non-empty ``geo_suppl_template`` are included.

        Args:
            gse_id: GEO Series accession (e.g. ``"GSE123456"``).

        Returns:
            List of ``(Mirror, URL)`` tuples in priority order.
        """
        geo_prefix, _dir_prefix = self._parse_geo_id(gse_id)
        results: list[tuple[Mirror, str]] = []
        for mirror in self._mirrors:
            if not mirror.geo_suppl_template:
                continue
            url = mirror.geo_suppl_template.format(
                base=mirror.base_url.rstrip("/"),
                prefix=geo_prefix,
                gse_id=gse_id.upper(),
            )
            results.append((mirror, url))
        return results

    def resolve_geo_file_url(
        self,
        gse_id: str,
        filename: str,
    ) -> list[tuple[Mirror, str]]:
        """Return ordered list of ``(mirror, url)`` for a specific GEO supplementary file.

        Args:
            gse_id: GEO Series accession.
            filename: File name on the supplementary directory.

        Returns:
            List of ``(Mirror, URL)`` tuples in priority order.
        """
        geo_prefix, _dir_prefix = self._parse_geo_id(gse_id)
        results: list[tuple[Mirror, str]] = []
        for mirror in self._mirrors:
            if not mirror.geo_file_template:
                continue
            url = mirror.geo_file_template.format(
                base=mirror.base_url.rstrip("/"),
                prefix=geo_prefix,
                gse_id=gse_id.upper(),
                filename=filename,
            )
            results.append((mirror, url))
        return results

    # ------------------------------------------------------------------
    # SRA URL resolution
    # ------------------------------------------------------------------

    def resolve_sra_url(self, srr_id: str) -> list[tuple[Mirror, str]]:
        """Return ordered list of ``(mirror, url)`` for an SRA run accession.

        Args:
            srr_id: SRA run accession (e.g. ``"SRR12345678"``).

        Returns:
            List of ``(Mirror, URL)`` tuples in priority order.
        """
        srr_id = srr_id.strip().upper()
        results: list[tuple[Mirror, str]] = []
        prefix_match = self._SRR_PREFIX_RE.match(srr_id)
        srr_prefix = prefix_match.group(1) if prefix_match else srr_id[:6]

        for mirror in self._mirrors:
            if not mirror.sra_template:
                continue
            url = mirror.sra_template.format(
                base=mirror.base_url.rstrip("/"),
                srr=srr_id,
                srr_prefix=srr_prefix,
                srr_prefix2=srr_id[:7] if len(srr_id) > 7 else srr_id,
            )
            results.append((mirror, url))
        return results

    # ------------------------------------------------------------------
    # Health checking
    # ------------------------------------------------------------------

    async def check_mirror_health(
        self,
        gse_id: str | None = None,
        sources: list[MirrorSource] | None = None,
    ) -> dict[str, bool]:
        """Check which mirrors are currently reachable.

        Sends a lightweight HEAD request to each mirror's base URL (or a
        specific dataset directory if ``gse_id`` is provided). Results are
        cached for 300 seconds to avoid repeated checks during a single
        processing session.

        Args:
            gse_id: Optional GEO accession — if provided, checks the specific
                supplementary directory for this dataset on each mirror.
            sources: Specific mirrors to check. If None, checks all.

        Returns:
            Dict mapping source name (e.g. ``"ncbi"``) → reachable (bool).
        """
        now = time.monotonic()

        # Filter mirrors to check
        to_check: list[Mirror] = []
        if sources:
            source_set = {s.value for s in sources}
            to_check = [m for m in self._mirrors if m.source.value in source_set]
        else:
            to_check = list(self._mirrors)

        # Determine check URLs
        check_urls: list[tuple[str, str]] = []  # (source_name, url)
        for mirror in to_check:
            cache_key = mirror.source.value
            cached = self._health_cache.get(cache_key)
            if cached and (now - cached[1]) < self._health_cache_ttl:
                continue  # Use cached result

            if gse_id and mirror.geo_suppl_template:
                try:
                    urls = self.resolve_geo_suppl_url(gse_id)
                    for src, url in urls:
                        if src.source == mirror.source:
                            check_urls.append((mirror.source.value, url))
                            break
                except ValueError:
                    check_urls.append((mirror.source.value, mirror.base_url))
            else:
                check_urls.append((mirror.source.value, mirror.base_url))

        if not check_urls:
            # All results cached — return cached values
            results: dict[str, bool] = {}
            for mirror in to_check:
                cached = self._health_cache.get(mirror.source.value)
                results[mirror.source.value] = cached[0] if cached else False
            return results

        # Run HEAD requests concurrently
        async with self._get_http_client() as client:
            tasks = []
            for source_name, url in check_urls:
                tasks.append(self._check_one(client, source_name, url))
            await asyncio.gather(*tasks, return_exceptions=True)

        # Build full result set
        results = {}
        for mirror in self._mirrors:
            if sources and mirror.source not in (sources if isinstance(sources, set) else set(sources)):
                continue
            cached = self._health_cache.get(mirror.source.value)
            results[mirror.source.value] = cached[0] if cached else False

        return results

    async def _check_one(
        self,
        client: "httpx.AsyncClient",
        source_name: str,
        url: str,
    ) -> None:
        """Check a single mirror URL and cache the result."""
        import httpx
        try:
            response = await client.head(
                url,
                follow_redirects=True,
            )
            reachable = response.status_code < 500
        except (httpx.HTTPError, httpx.TimeoutException, OSError):
            reachable = False
        except Exception:
            reachable = False

        self._health_cache[source_name] = (reachable, time.monotonic())
        self._logger.debug(
            "Mirror health check",
            source=source_name,
            url=url,
            reachable=reachable,
        )

    async def get_healthiest_mirror(
        self,
        gse_id: str | None = None,
    ) -> Mirror:
        """Return the highest-priority mirror that passes a health check.

        Falls back to NCBI if no mirror is reachable (NCBI is the source of
        truth — if it is unreachable too, the caller must handle the error).

        Args:
            gse_id: Optional GEO accession for targeted health check.

        Returns:
            The healthiest ``Mirror`` (never None — always defaults to NCBI).
        """
        health = await self.check_mirror_health(gse_id=gse_id)
        for mirror in self._mirrors:
            if health.get(mirror.source.value, False):
                self._logger.info(
                    "Selected healthiest mirror",
                    source=mirror.source.value,
                    base_url=mirror.base_url,
                )
                return mirror

        # Fallback to NCBI — source of truth
        ncbi = next(
            (m for m in self._mirrors if m.source == MirrorSource.NCBI),
            self._mirrors[0],
        )
        self._logger.warning(
            "No healthy mirror found, falling back to NCBI",
            checked=health,
        )
        return ncbi

    def _get_http_client(self) -> "httpx.AsyncClient":
        """Create a shared httpx client for health checks."""
        import httpx
        return httpx.AsyncClient(timeout=10.0, follow_redirects=True)

    async def aclose(self) -> None:
        """No-op — httpx clients are created per-check and auto-closed."""
        pass


# ---------------------------------------------------------------------------
# Standalone convenience functions
# ---------------------------------------------------------------------------


async def check_mirror_health(
    gse_id: str | None = None,
) -> dict[str, bool]:
    """Convenience: check mirror reachability without managing a resolver.

    Args:
        gse_id: Optional GEO accession for targeted health check.

    Returns:
        Dict mapping mirror source name → reachable.
    """
    resolver = MirrorResolver()
    return await resolver.check_mirror_health(gse_id=gse_id)


def resolve_geo_urls(
    gse_id: str,
    filename: str,
) -> list[dict[str, str]]:
    """Convenience: get all mirror URLs for a GEO file as plain dicts.

    Args:
        gse_id: GEO Series accession.
        filename: Supplementary file name.

    Returns:
        List of ``{"source": "ncbi", "url": "https://..."}`` dicts.
    """
    resolver = MirrorResolver()
    results = resolver.resolve_geo_file_url(gse_id, filename)
    return [
        {"source": mirror.source.value, "url": url}
        for mirror, url in results
    ]
