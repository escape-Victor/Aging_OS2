"""Database metadata probe — queries GEO/GSM/SRA for study context.

Extracts: species, tissue, disease, treatment, platform, protocol,
sample groups, and publication information into structured JSON.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from datetime import datetime, timezone

import structlog

logger = structlog.get_logger(__name__)


@dataclass
class SampleInfo:
    """Structured information about a single GSM sample."""

    gsm: str
    title: str = ""
    organism: str = ""
    source_name: str = ""
    characteristics: dict[str, str] = field(default_factory=dict)
    sra_runs: list[str] = field(default_factory=list)
    platform: str = ""
    molecule: str = ""


@dataclass
class StudyContext:
    """Full study context extracted from GEO/SRA/publication metadata."""

    dataset_id: str
    title: str = ""
    summary: str = ""
    organism: list[str] = field(default_factory=list)
    platform: list[str] = field(default_factory=list)
    technology: str = ""
    tissue_or_organ: list[str] = field(default_factory=list)
    disease: list[str] = field(default_factory=list)
    condition_groups: list[str] = field(default_factory=list)
    samples: list[SampleInfo] = field(default_factory=list)
    publication: dict[str, str] = field(default_factory=dict)
    confidence: dict[str, float] = field(default_factory=dict)
    raw_geo_metadata: dict[str, object] = field(default_factory=dict)
    probed_at: str = field(default_factory=lambda: datetime.now(timezone.utc).isoformat())


class MetadataProbe:
    """Probes GEO, SRA, and publication databases for study context.

    This is the first intelligence layer: it reads web/database metadata
    to understand the biological experiment before touching any data files.
    """

    def __init__(
        self,
        email: str | None = None,
        api_key: str | None = None,
        use_cache: bool = True,
        cache_dir: str | None = None,
        timeout: float = 60.0,
        proxy: str | None = None,
    ) -> None:
        self.email = email
        self.api_key = api_key
        self.use_cache = use_cache
        self.cache_dir = cache_dir
        self.timeout = timeout
        self.proxy = proxy
        self._logger = logger.bind(component="MetadataProbe")

    async def probe_geo(
        self,
        gse_id: str,
        supplementary_filenames: list[str] | None = None,
    ) -> StudyContext:
        """Query GEO for a GSE accession and extract study context.

        Uses GEOparse (or equivalent) to fetch the SOFT record and extract
        structured metadata: species, tissue, disease, platform, samples.

        Args:
            gse_id: GEO Series accession (e.g., 'GSE123456').
            supplementary_filenames: Optional list of supplementary filenames
                from GEO FTP for 3rd-layer technology detection.

        Returns:
            StudyContext with extracted metadata and confidence scores.
        """
        self._logger.info("Probing GEO metadata", gse_id=gse_id)
        try:
            import GEOparse  # type: ignore[import-untyped]

            # Configure network settings before GEOparse downloads
            _old_timeout = None
            try:
                import socket
                _old_timeout = socket.getdefaulttimeout()
                socket.setdefaulttimeout(self.timeout)
            except Exception:
                pass

            if self.proxy:
                import os as _os_proxy
                _os_proxy.environ.setdefault("HTTP_PROXY", self.proxy)
                _os_proxy.environ.setdefault("HTTPS_PROXY", self.proxy)

            try:
                gse = GEOparse.get_GEO(geo=gse_id, silent=True)
            finally:
                if _old_timeout is not None:
                    try:
                        socket.setdefaulttimeout(_old_timeout)
                    except Exception:
                        pass
        except Exception as exc:
            self._logger.warning("GEO query failed", gse_id=gse_id, error=str(exc))
            return StudyContext(
                dataset_id=gse_id,
                title=f"GEO Series {gse_id}",
                summary=f"(GEO query failed: {exc})",
                organism=[],
                confidence={"organism": 0.0, "overall": 0.0},
            )

        # Extract basic metadata
        metadata = gse.metadata
        title = metadata.get("title", [f"GEO Series {gse_id}"])[0]
        summary = metadata.get("summary", [""])[0]
        organism = [o.strip() for o in metadata.get("organism", [])]

        # Fallback 1: use taxid if organism field is empty
        if not organism:
            organism = self._taxids_to_organisms(metadata)

        # Fallback 2: use supplementary_file or summary text hints
        if not organism:
            organism = self._guess_organism_from_text(title, summary)

        # Extract platform info
        platforms = []
        for platform_id in metadata.get("platform_id", []):
            platforms.append(platform_id)

        # Extract samples
        samples: list[SampleInfo] = []
        for gsm_name, gsm in gse.gsms.items():
            gsm_metadata = gsm.metadata
            sample_org = (
                gsm_metadata.get("organism", [""])[0]
                if gsm_metadata.get("organism")
                else ""
            )
            # Fallback: use organism inferred at series level
            if not sample_org and organism:
                sample_org = organism[0]

            sample = SampleInfo(
                gsm=gsm_name,
                title=gsm_metadata.get("title", [""])[0],
                organism=sample_org,
                source_name=(
                    gsm_metadata.get("source_name_ch1", [""])[0]
                    if gsm_metadata.get("source_name_ch1")
                    else ""
                ),
                characteristics=self._parse_characteristics(gsm_metadata),
                platform=(
                    gsm_metadata.get("platform_id", [""])[0]
                    if gsm_metadata.get("platform_id")
                    else ""
                ),
                molecule=(
                    gsm_metadata.get("molecule_ch1", [""])[0]
                    if gsm_metadata.get("molecule_ch1")
                    else ""
                ),
            )
            # Extract SRA runs from relations
            for rel in gsm_metadata.get("relation", []):
                if rel.startswith("SRA:"):
                    srr_id = rel.split("SRA:", 1)[1].strip().split("/")[-1]
                    sample.sra_runs.append(srr_id)
            samples.append(sample)

        # Fallback 3: aggregate organisms from samples
        if not organism:
            sample_orgs: set[str] = set()
            for s in samples:
                if s.organism:
                    sample_orgs.add(s.organism)
            organism = sorted(sample_orgs)

        # Infer technology from platforms + molecule + supplementary filenames
        technology = self._infer_technology(platforms, samples, supplementary_filenames)

        # Build confidence scores
        confidence = {
            "organism": 0.95 if organism else 0.0,
            "platform": 0.90 if platforms else 0.0,
            "samples": 0.98 if samples else 0.0,
            "overall": 0.85 if organism and samples else 0.3,
        }

        return StudyContext(
            dataset_id=gse_id,
            title=title,
            summary=summary,
            organism=organism,
            platform=platforms,
            technology=technology,
            samples=samples,
            confidence=confidence,
            raw_geo_metadata=dict(metadata),
        )

    async def probe_sra(self, srr_ids: list[str]) -> dict[str, object]:
        """Query SRA for run metadata (platform, layout, spots, bases).

        Args:
            srr_ids: List of SRA run accessions.

        Returns:
            Dict mapping SRR -> metadata dict.
        """
        self._logger.info("Probing SRA metadata", count=len(srr_ids))
        if not srr_ids:
            return {}

        try:
            import xml.etree.ElementTree as ET

            import httpx
        except ImportError:
            self._logger.warning("httpx not available — skipping SRA probe")
            return {}

        results: dict[str, object] = {}
        # Query in batches of 20 to avoid overly long URLs
        batch_size = 20
        for i in range(0, len(srr_ids), batch_size):
            batch = srr_ids[i : i + batch_size]
            ids_param = ",".join(batch)
            url = (
                "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
                f"?db=sra&id={ids_param}&rettype=full&retmode=xml"
            )
            try:
                async with httpx.AsyncClient(timeout=30.0) as client:
                    response = await client.get(url)
                    if response.status_code == 200:
                        root = ET.fromstring(response.text)
                        for run in root.iter("RUN"):
                            srr = run.get("accession", "")
                            if not srr:
                                continue
                            run_info: dict[str, object] = {}
                            total_spots = run.get("total_spots", "")
                            total_bases = run.get("total_bases", "")
                            if total_spots:
                                run_info["total_spots"] = int(total_spots)
                            if total_bases:
                                run_info["total_bases"] = int(total_bases)
                            # Extract platform
                            platform_el = run.find(".//PLATFORM")
                            if platform_el is not None and len(platform_el) > 0:
                                run_info["platform"] = platform_el[0].tag
                            # Extract layout
                            layout_el = run.find(".//LIBRARY_LAYOUT")
                            if layout_el is not None and len(layout_el) > 0:
                                run_info["layout"] = layout_el[0].tag  # PAIRED or SINGLE
                            results[srr] = run_info
                    else:
                        self._logger.warning(
                            "SRA query returned non-200", status=response.status_code
                        )
            except Exception as exc:
                self._logger.warning("SRA query failed for batch", error=str(exc))
        return results

    async def infer_condition_groups(self, ctx: StudyContext) -> list[str]:
        """Infer condition/experimental groups from sample characteristics.

        Analyzes sample metadata to detect treatment vs control, disease vs
        normal, time points, doses, etc.

        Args:
            ctx: StudyContext with populated samples.

        Returns:
            List of condition group labels.
        """
        self._logger.info("Inferring condition groups", sample_count=len(ctx.samples))
        groups: list[str] = []
        seen: set[str] = set()

        # Common characteristic keys that encode experimental conditions
        condition_keys = {
            "treatment", "condition", "disease", "genotype", "strain",
            "time point", "time", "dose", "tissue", "cell type",
            "cell line", "agent", "stimulation", "developmental stage",
            "age", "gender", "sex",
        }

        for sample in ctx.samples:
            for key, value in sample.characteristics.items():
                key_lower = key.lower().strip()
                # Check if this key is a condition-relevant one
                is_condition = any(ck in key_lower for ck in condition_keys)
                # Also include keys with known condition patterns
                if is_condition or any(
                    kw in key_lower
                    for kw in ["treatment", "disease", "condition",
                               "genotype", "strain", "time"]
                ):
                    label = f"{key}:{value}"
                    if label not in seen:
                        seen.add(label)
                        groups.append(label)

        # If no explicit groups found, fall back to source_name or title
        if not groups:
            for sample in ctx.samples:
                if sample.source_name and sample.source_name not in seen:
                    seen.add(sample.source_name)
                    groups.append(sample.source_name)

        return groups

    async def resolve_organism(self, ctx: StudyContext) -> dict[str, float]:
        """Resolve organism(s) with confidence scores.

        Cross-references GEO metadata with gene symbols in processed files
        to validate organism assignments.

        Args:
            ctx: StudyContext with raw metadata.

        Returns:
            Dict mapping organism name -> confidence score.
        """
        self._logger.info("Resolving organism")
        if not ctx.organism:
            # Try to infer from sample metadata
            for sample in ctx.samples:
                if sample.organism:
                    ctx.organism = [sample.organism]
                    return {sample.organism: 0.70}
            return {"unknown": 0.0}

        results: dict[str, float] = {}
        for org in ctx.organism:
            org_lower = org.lower().strip()
            if any(t in org_lower for t in ("homo sapiens", "human", "hsapiens")):
                results["Homo sapiens"] = 0.95
            elif any(t in org_lower for t in ("mus musculus", "mouse", "mmusculus")):
                results["Mus musculus"] = 0.95
            elif any(t in org_lower for t in ("danio rerio", "zebrafish")):
                results["Danio rerio"] = 0.90
            elif any(t in org_lower for t in ("rattus", "rat")):
                results[org] = 0.85
            else:
                results[org] = 0.70
        return results

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------

    # NCBI Taxonomy ID → species name mapping (commonly seen in scRNA-seq)
    _TAXID_MAP: dict[str, str] = {
        "9606": "Homo sapiens",
        "10090": "Mus musculus",
        "10116": "Rattus norvegicus",
        "7955": "Danio rerio",
        "6239": "Caenorhabditis elegans",
        "7227": "Drosophila melanogaster",
        "9544": "Macaca mulatta",
        "9913": "Bos taurus",
        "9031": "Gallus gallus",
        "9823": "Sus scrofa",
        "9615": "Canis lupus familiaris",
        "9685": "Felis catus",
        "9598": "Pan troglodytes",
        "8364": "Xenopus tropicalis",
    }

    @staticmethod
    def _taxids_to_organisms(metadata: dict) -> list[str]:
        """Resolve organism from NCBI taxonomy IDs in GEO metadata.

        Checks platform_taxid, sample_taxid fields, and taxid in characteristics.
        """
        organisms: list[str] = []
        for taxid_field in ("platform_taxid", "sample_taxid"):
            for taxid in metadata.get(taxid_field, []):
                org = MetadataProbe._TAXID_MAP.get(str(taxid))
                if org and org not in organisms:
                    organisms.append(org)
        return organisms

    @staticmethod
    def _guess_organism_from_text(title: str, summary: str) -> list[str]:
        """Guess organism from study title/summary text as last-resort fallback."""
        combined = f"{title} {summary}".lower()
        organisms: list[str] = []
        text_hints = [
            ("human", "Homo sapiens"),
            ("homo sapiens", "Homo sapiens"),
            ("patient", "Homo sapiens"),  # "patient" typically implies human
            ("mouse", "Mus musculus"),
            ("mus musculus", "Mus musculus"),
            ("murine", "Mus musculus"),
            ("zebrafish", "Danio rerio"),
            ("danio rerio", "Danio rerio"),
            ("rat", "Rattus norvegicus"),
            ("rattus", "Rattus norvegicus"),
            ("drosophila", "Drosophila melanogaster"),
            ("c. elegans", "Caenorhabditis elegans"),
            ("macaque", "Macaca mulatta"),
            ("xenopus", "Xenopus tropicalis"),
        ]
        for keyword, species in text_hints:
            if keyword in combined and species not in organisms:
                organisms.append(species)
        return organisms

    @staticmethod
    def _parse_characteristics(metadata: dict) -> dict[str, str]:
        """Parse GEO sample characteristics into a key:value dict."""
        characteristics: dict[str, str] = {}
        for key, values in metadata.items():
            if key.startswith("characteristics_ch") and values:
                for entry in values:
                    if ":" in entry:
                        k, v = entry.split(":", 1)
                        characteristics[k.strip()] = v.strip()
                    else:
                        characteristics[key] = entry
        return characteristics

    @staticmethod
    def _infer_technology(
        platforms: list[str],
        samples: list[SampleInfo],
        supplementary_filenames: list[str] | None = None,
    ) -> str:
        """Infer the sequencing technology from platform, sample, and file evidence.

        Uses a three-layer approach:
        1. Platform ID → known technology mapping
        2. Sample molecule field → protocol detection
        3. Supplementary filenames → format evidence (NEW)

        If both layers agree on 10x, returns the specific 10x chemistry.
        If platform says generic Illumina but molecule says 10x, returns 10x
        (many 10x datasets use generic Illumina platform IDs like GPL16791).
        If platform and molecule are both non-specific, checks filenames for
        10x patterns (barcodes.tsv, features.tsv, matrix.mtx, etc.).
        """
        platform_tech_map: dict[str, str] = {
            "GPL24676": "10x Genomics 3' v3",
            "GPL18573": "10x Genomics 3' v2",
            "GPL21103": "10x Genomics 5'",
            "GPL19415": "10x Genomics",
            "GPL20301": "10x Genomics",
            "GPL17021": "Illumina HiSeq",
            "GPL16791": "Illumina HiSeq",
            "GPL11154": "Illumina HiSeq",
            "GPL19057": "Illumina NextSeq",
            "GPL21697": "Illumina NovaSeq",
            "GPL24247": "Illumina NovaSeq",
            "GPL16417": "Illumina HiSeq",
            "GPL18480": "Illumina NextSeq",
            "GPL21273": "Illumina MiSeq",
            "GPL15520": "Illumina MiSeq",
            "GPL21493": "Illumina HiSeq",
            "GPL13112": "Illumina HiSeq",
        }

        # Layer 1: platform-based detection
        platform_tech: str | None = None
        for platform_id in platforms:
            if platform_id in platform_tech_map:
                platform_tech = platform_tech_map[platform_id]
                break

        # Layer 2: molecule-based detection
        molecule_tech: str | None = None
        for sample in samples:
            mol = sample.molecule.lower()
            if "10x" in mol or "chromium" in mol:
                molecule_tech = "10x Genomics"
                break
            if "smart-seq" in mol:
                molecule_tech = "Smart-seq2"
                break
            if "drop-seq" in mol or "dropseq" in mol:
                molecule_tech = "Drop-seq"
                break

        # Combine evidence: molecule evidence overrides generic Illumina labels
        if molecule_tech:
            if platform_tech and "illumina" in platform_tech.lower():
                # Platform says generic Illumina, molecule says 10x/Smart-seq2/Drop-seq
                # → use the molecule-based detection (more specific)
                return molecule_tech
            if platform_tech and "10x" in platform_tech.lower():
                # Both agree on 10x — use the platform-specific chemistry
                return platform_tech
            # No platform info or platform unknown → use molecule
            return molecule_tech

        # Layer 3: Supplementary filename evidence
        # When both platform and molecule are non-specific, check GEO
        # supplementary filenames for 10x format markers.
        if supplementary_filenames:
            tenx_patterns = [
                "_10X_", "_10x_", "barcodes.tsv", "features.tsv",
                "matrix.mtx", "filtered_gene_bc_matrices",
                "filtered_feature_bc_matrix",
            ]
            pattern_hits = sum(
                1 for f in supplementary_filenames
                if any(p in f for p in tenx_patterns)
            )
            if pattern_hits >= 2:
                # Combine with existing evidence
                if platform_tech and "illumina" in platform_tech.lower():
                    return "10x Genomics"
                if not platform_tech:
                    return "10x Genomics"
                # Platform says something else, but filenames strongly indicate 10x
                if "10x" not in platform_tech.lower():
                    return "10x Genomics"

        # No molecule evidence → fall back to platform or generic label
        return platform_tech or "scRNA-seq"

    async def probe_all(
        self,
        gse_id: str,
        supplementary_filenames: list[str] | None = None,
    ) -> StudyContext:
        """Run all metadata probes for a GSE accession.

        This is the main entry point: query GEO for the study record,
        extract sample info, probe SRA for run details, and infer
        experimental groups.

        Args:
            gse_id: GEO Series accession.
            supplementary_filenames: Optional supplementary filenames for
                technology detection.

        Returns:
            Fully populated StudyContext.
        """
        ctx = await self.probe_geo(gse_id, supplementary_filenames)
        if ctx.samples:
            srr_ids = [run for s in ctx.samples for run in s.sra_runs]
            if srr_ids:
                await self.probe_sra(srr_ids)
        if ctx.samples and not ctx.condition_groups:
            ctx.condition_groups = await self.infer_condition_groups(ctx)
        # Optional LLM-powered enhancement: infer experimental design
        if ctx.samples:
            try:
                titles = {s.gsm: s.title for s in ctx.samples if s.title and s.gsm}
                if titles and len(titles) >= 4:  # need min samples for meaningful inference
                    llm_conditions = await self._llm_infer_conditions(titles, gse_id)
                    if llm_conditions:
                        ctx.condition_groups = list(llm_conditions.keys())
                        ctx.raw_geo_metadata["llm_inferred_conditions"] = llm_conditions
            except Exception:
                pass  # LLM is best-effort; never blocks the pipeline
        return ctx

    async def _llm_infer_conditions(
        self,
        sample_titles: dict[str, str],
        dataset_id: str,
    ) -> dict[str, list[str]] | None:
        """Use LLM to infer experimental design from sample titles.

        This is an OPTIONAL enhancement that uses an LLM (e.g., DeepSeek)
        to classify samples into control vs treatment groups based on
        their human-readable titles. Best-effort only — failures are
        silently caught upstream.

        Ported from v1 ``DynamicDesignMatrixGenerator.infer_from_metadata()``.

        Requires: ``DEEPSEEK_API_KEY`` in environment.

        Args:
            sample_titles: Dict mapping sample IDs (GSMxxx) to titles.
            dataset_id: GEO accession for context.

        Returns:
            Dict with ``"control"`` and ``"treatment"`` keys, or None on failure.
        """
        import json
        import os

        api_key = os.getenv("DEEPSEEK_API_KEY")
        if not api_key:
            self._logger.debug("No DEEPSEEK_API_KEY — skipping LLM condition inference")
            return None

        try:
            from openai import OpenAI
        except ImportError:
            self._logger.debug("openai package not installed — skipping LLM inference")
            return None

        system_prompt = (
            "You are a strict bioinformatics experimental design inferencer. "
            "The user will provide a dictionary mapping single-cell sample IDs "
            "to their biological descriptions (titles). "
            "Your task is to identify the baseline \"Control\" (e.g., young, "
            "wild-type, healthy, mock) and the \"Treatment\" (e.g., aged, "
            "diseased, perturbed, knockout) groups.\n\n"
            "RULES:\n"
            "1. Return strictly formatted JSON.\n"
            '2. Schema: {"control": ["GSMxxx", ...], "treatment": ["GSMxxx", ...]}\n'
            "3. Ignore irrelevant samples or technical replicates.\n"
            "4. If a definitive contrast cannot be inferred, return empty arrays."
        )

        user_prompt = (
            f"Dataset: {dataset_id}\n"
            f"Sample Metadata Mapping:\n"
            f"{json.dumps(sample_titles, indent=2)}"
        )

        try:
            client = OpenAI(api_key=api_key, base_url="https://api.deepseek.com")
            response = client.chat.completions.create(
                model="deepseek-chat",
                messages=[
                    {"role": "system", "content": system_prompt},
                    {"role": "user", "content": user_prompt},
                ],
                response_format={"type": "json_object"},
                temperature=0.1,
            )
            design_matrix = json.loads(response.choices[0].message.content)
            control_count = len(design_matrix.get("control", []))
            treat_count = len(design_matrix.get("treatment", []))

            if control_count == 0 or treat_count == 0:
                self._logger.info(
                    "LLM inference found insufficient contrast",
                    dataset_id=dataset_id,
                )
                return None

            self._logger.info(
                "LLM inferred experimental design",
                dataset_id=dataset_id,
                n_control=control_count,
                n_treatment=treat_count,
            )
            return design_matrix

        except Exception as exc:
            self._logger.debug("LLM condition inference failed", error=str(exc))
            return None

    @staticmethod
    def detect_species_from_gene_symbols(
        var_names: list[str],
        sample_size: int = 5000,
    ) -> str:
        """Detect species from gene symbol casing patterns.

        Human genes are typically UPPERCASE (BRCA1, TP53, ACTB).
        Mouse genes are typically Titlecase (Brca1, Tp53, Actb).
        This provides an orthogonal cross-validation to GEO metadata,
        which can occasionally be mislabeled.

        Ported from v1 ``UniversalOmicsIntegrator._detect_species()``.

        Algorithm:
        1. Sample up to ``sample_size`` gene symbols.
        2. Count uppercase-only genes (len > 2, all uppercase letters).
        3. Count titlecase genes (len > 2, first letter uppercase).
        4. Compute ``upper_ratio = uppercase / (uppercase + titlecase)``.
        5. If ``upper_ratio > 0.65`` → ``"Human"``.
        6. If ``upper_ratio < 0.35`` → ``"Mouse"``.
        7. Tiebreakers: MT- prefix counts, RPL prefix counts.
        8. Default → ``"Mixed_or_Unknown"``. Returns ``"Unknown"`` if
           no valid gene symbols are found.

        Args:
            var_names: List of gene names (e.g., from ``adata.var_names``).
            sample_size: Maximum genes to sample for performance.

        Returns:
            One of ``"Human"``, ``"Mouse"``, ``"Mixed_or_Unknown"``,
            or ``"Unknown"``.
        """
        if not var_names:
            return "Unknown"

        # Sample for performance
        sample_genes = [str(g) for g in var_names[:sample_size]]

        # Count casing patterns (only genes longer than 2 chars)
        # NOTE: use isalnum() not isalpha() — many gene symbols contain digits (BRCA1, TP53)
        # NOTE: exclude Ensembl IDs (ENSG..., ENSMUSG...) from casing analysis
        def _is_valid_gene_symbol(g: str) -> bool:
            return (
                not g.upper().startswith("ENS")
                and g.replace("-", "").replace(".", "").isalnum()
            )

        upper_count = sum(
            1 for g in sample_genes
            if len(g) > 2 and g.isupper() and _is_valid_gene_symbol(g)
        )
        title_count = sum(
            1 for g in sample_genes
            if len(g) > 2 and g[0].isupper() and g[1:].islower()
            and _is_valid_gene_symbol(g)
        )

        total_valid = upper_count + title_count
        if total_valid == 0:
            return "Unknown"

        upper_ratio = upper_count / total_valid

        if upper_ratio > 0.65:
            return "Human"
        elif upper_ratio < 0.35:
            return "Mouse"
        else:
            # Tiebreakers
            has_mt_prefix = sum(
                1 for g in sample_genes
                if str(g).upper().startswith("MT-")
            )
            has_rpl_prefix = sum(
                1 for g in sample_genes
                if str(g).upper().startswith("RPL")
            )

            if has_mt_prefix > len(sample_genes) * 0.05:
                return "Human"
            elif has_rpl_prefix > len(sample_genes) * 0.03:
                return "Human"  # RPL genes are human naming convention

            return "Mixed_or_Unknown"
