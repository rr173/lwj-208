from typing import List, Dict, Optional
from sqlalchemy.orm import Session
from app import models, schemas


DOMAIN_RISK_MAP = {
    "catalytic": "高风险",
    "binding": "中风险",
    "structural": "低风险",
    "regulatory": "中风险",
}

RISK_ORDER = {"高风险": 3, "中风险": 2, "低风险": 1, "未知风险": 0}

MISSENSE_CONSEQUENCES = {"missense_variant", "stop_gained", "stop_lost"}


def register_domains(
    db: Session, request: schemas.ProteinDomainBatchRequest
) -> List[schemas.ProteinDomainOut]:
    created = []
    for item in request.domains:
        if item.start_codon > item.end_codon:
            continue
        domain = models.ProteinDomain(
            gene_name=item.gene_name,
            domain_name=item.domain_name,
            start_codon=item.start_codon,
            end_codon=item.end_codon,
            domain_type=item.domain_type,
        )
        db.add(domain)
        created.append(domain)
    db.commit()
    for d in created:
        db.refresh(d)
    return [
        schemas.ProteinDomainOut(
            id=d.id,
            gene_name=d.gene_name,
            domain_name=d.domain_name,
            start_codon=d.start_codon,
            end_codon=d.end_codon,
            domain_type=d.domain_type,
        )
        for d in created
    ]


def query_domains_by_gene(
    db: Session, gene_name: str
) -> List[schemas.ProteinDomainOut]:
    rows = (
        db.query(models.ProteinDomain)
        .filter(models.ProteinDomain.gene_name == gene_name)
        .order_by(models.ProteinDomain.start_codon)
        .all()
    )
    return [
        schemas.ProteinDomainOut(
            id=r.id,
            gene_name=r.gene_name,
            domain_name=r.domain_name,
            start_codon=r.start_codon,
            end_codon=r.end_codon,
            domain_type=r.domain_type,
        )
        for r in rows
    ]


def _compute_codon_position(ref_pos: int, gene_exons: List, strand: str = "+") -> Optional[int]:
    """
    计算变异相对于整个 CDS 的密码子位置（1-based）。
    gene_exons: 该基因所有的 exon GeneAnnotation 记录列表，每条有 start, end。
    算法: 先按坐标排序 (+ strand 升序, - strand 降序), 找到 ref_pos 所在 exon,
    累加之前所有 exon 的密码子数, 再加上当前 exon 内部的 (offset // 3 + 1)。
    """
    if not gene_exons:
        return None

    exon_list = [e for e in gene_exons]
    if strand == "-":
        exon_list.sort(key=lambda e: e.end, reverse=True)
    else:
        exon_list.sort(key=lambda e: e.start)

    target_exon_idx = -1
    for i, exon in enumerate(exon_list):
        if exon.start <= ref_pos <= exon.end:
            target_exon_idx = i
            break
    if target_exon_idx < 0:
        return None

    codon_count_before = 0
    for i in range(target_exon_idx):
        exon = exon_list[i]
        exon_len = exon.end - exon.start + 1
        codon_count_before += exon_len // 3

    target_exon = exon_list[target_exon_idx]
    if strand == "-":
        offset_in_exon = target_exon.end - ref_pos
    else:
        offset_in_exon = ref_pos - target_exon.start
    codon_in_exon = offset_in_exon // 3 + 1

    return codon_count_before + codon_in_exon


def _find_matching_domains(
    codon_pos: int, domains: List[models.ProteinDomain]
) -> List[schemas.HitDomainInfo]:
    hits = []
    for d in domains:
        if d.start_codon <= codon_pos <= d.end_codon:
            hits.append(
                schemas.HitDomainInfo(
                    domain_name=d.domain_name,
                    domain_type=d.domain_type,
                    start_codon=d.start_codon,
                    end_codon=d.end_codon,
                )
            )
    return hits


def _assess_risk(
    consequence: Optional[str],
    hit_domains: List[schemas.HitDomainInfo],
    has_gene_domains: bool,
) -> tuple:
    if not hit_domains:
        if not has_gene_domains:
            return "未知风险", "无结构域信息"
        return "未知风险", "不在任何结构域内"

    is_missense_or_nonsense = (
        consequence is not None and consequence in MISSENSE_CONSEQUENCES
    )

    if not is_missense_or_nonsense:
        return "未知风险", "命中结构域但非错义/无义突变"

    best_risk = "未知风险"
    best_reason_parts = []
    for hd in hit_domains:
        risk = DOMAIN_RISK_MAP.get(hd.domain_type, "未知风险")
        if RISK_ORDER.get(risk, 0) > RISK_ORDER.get(best_risk, 0):
            best_risk = risk
        best_reason_parts.append(f"{hd.domain_name}({hd.domain_type})")

    reason = f"{consequence}, 命中结构域: {', '.join(best_reason_parts)}"
    return best_risk, reason


def map_sample_variants_to_domains(
    db: Session, sample_id: int
) -> schemas.SampleImpactAssessmentOut:
    sample = db.query(models.Sample).filter(models.Sample.id == sample_id).first()
    if not sample:
        return None

    exon_variants = (
        db.query(models.SampleVariantSpectrum)
        .filter(
            models.SampleVariantSpectrum.sample_id == sample_id,
            models.SampleVariantSpectrum.feature_type == "exon",
        )
        .order_by(
            models.SampleVariantSpectrum.ref_pos,
        )
        .all()
    )

    gene_names = list({v.gene_name for v in exon_variants if v.gene_name})

    gene_exons: Dict[str, List[models.GeneAnnotation]] = {}
    gene_strands: Dict[str, str] = {}
    if gene_names:
        ann_rows = (
            db.query(models.GeneAnnotation)
            .filter(models.GeneAnnotation.gene_name.in_(gene_names))
            .all()
        )
        for ann in ann_rows:
            if ann.gene_name not in gene_exons:
                gene_exons[ann.gene_name] = []
                gene_strands[ann.gene_name] = ann.strand or "+"
            if ann.feature_type == "exon":
                gene_exons[ann.gene_name].append(ann)

    gene_domains: Dict[str, List[models.ProteinDomain]] = {}
    if gene_names:
        domain_rows = (
            db.query(models.ProteinDomain)
            .filter(models.ProteinDomain.gene_name.in_(gene_names))
            .all()
        )
        for d in domain_rows:
            if d.gene_name not in gene_domains:
                gene_domains[d.gene_name] = []
            gene_domains[d.gene_name].append(d)

    results: List[schemas.VariantImpactOut] = []

    for v in exon_variants:
        gene_name = v.gene_name
        if not gene_name:
            continue

        exons_for_gene = gene_exons.get(gene_name, [])
        strand = gene_strands.get(gene_name, "+")
        if not exons_for_gene:
            continue

        codon_pos = _compute_codon_position(v.ref_pos, exons_for_gene, strand)
        if codon_pos is None:
            continue

        domains_for_gene = gene_domains.get(gene_name, [])
        has_gene_domains = len(domains_for_gene) > 0
        hit_domains = _find_matching_domains(codon_pos, domains_for_gene)

        aa_ref = v.aa_ref
        aa_alt = v.aa_alt
        amino_acid_change = None
        if aa_ref and aa_alt and (v.consequence in MISSENSE_CONSEQUENCES):
            amino_acid_change = f"{aa_ref}{codon_pos}{aa_alt}"

        risk_level, risk_reason = _assess_risk(
            v.consequence, hit_domains, has_gene_domains
        )

        results.append(
            schemas.VariantImpactOut(
                ref_pos=v.ref_pos,
                gene_name=gene_name,
                codon_position=codon_pos,
                hit_domains=hit_domains,
                aa_ref=aa_ref,
                aa_alt=aa_alt,
                amino_acid_change=amino_acid_change,
                consequence=v.consequence,
                variant_type=v.variant_type,
                ref_base=v.ref_base,
                alt_base=v.alt_base,
                risk_level=risk_level,
                risk_reason=risk_reason,
            )
        )

    results.sort(key=lambda r: RISK_ORDER.get(r.risk_level, 0), reverse=True)

    risk_summary: Dict[str, int] = {}
    for r in results:
        risk_summary[r.risk_level] = risk_summary.get(r.risk_level, 0) + 1

    variants_with_domain_hit = sum(1 for r in results if len(r.hit_domains) > 0)

    return schemas.SampleImpactAssessmentOut(
        sample_id=sample_id,
        sample_name=sample.name,
        total_exon_variants=len(exon_variants),
        mapped_variants=len(results),
        variants_with_domain_hit=variants_with_domain_hit,
        risk_summary=risk_summary,
        variants=results,
    )


def summarize_domain_variants_across_samples(
    db: Session, gene_name: str, domain_name: str
) -> Optional[schemas.DomainVariantSummaryOut]:
    domain = (
        db.query(models.ProteinDomain)
        .filter(
            models.ProteinDomain.gene_name == gene_name,
            models.ProteinDomain.domain_name == domain_name,
        )
        .first()
    )
    if not domain:
        return None

    gene_exons = (
        db.query(models.GeneAnnotation)
        .filter(
            models.GeneAnnotation.gene_name == gene_name,
            models.GeneAnnotation.feature_type == "exon",
        )
        .all()
    )
    if not gene_exons:
        return None

    strand = gene_exons[0].strand or "+"

    all_exon_variants = (
        db.query(models.SampleVariantSpectrum)
        .filter(
            models.SampleVariantSpectrum.gene_name == gene_name,
            models.SampleVariantSpectrum.feature_type == "exon",
        )
        .all()
    )

    sample_ids = {v.sample_id for v in all_exon_variants}
    samples = db.query(models.Sample).filter(models.Sample.id.in_(sample_ids)).all()
    sample_map = {s.id: s for s in samples}

    gene_domains = (
        db.query(models.ProteinDomain)
        .filter(models.ProteinDomain.gene_name == gene_name)
        .all()
    )
    has_gene_domains = len(gene_domains) > 0

    variant_groups: Dict[tuple, schemas.DomainVariantSummaryEntry] = {}

    for v in all_exon_variants:
        codon_pos = _compute_codon_position(v.ref_pos, gene_exons, strand)
        if codon_pos is None:
            continue
        if not (domain.start_codon <= codon_pos <= domain.end_codon):
            continue

        hit_domains = _find_matching_domains(codon_pos, gene_domains)

        aa_ref = v.aa_ref
        aa_alt = v.aa_alt
        amino_acid_change = None
        if aa_ref and aa_alt and (v.consequence in MISSENSE_CONSEQUENCES):
            amino_acid_change = f"{aa_ref}{codon_pos}{aa_alt}"

        risk_level, risk_reason = _assess_risk(
            v.consequence, hit_domains, has_gene_domains
        )

        variant_key = (
            v.ref_pos,
            v.variant_type,
            v.ref_base,
            v.alt_base,
            v.consequence,
        )

        if variant_key not in variant_groups:
            variant_groups[variant_key] = schemas.DomainVariantSummaryEntry(
                codon_position=codon_pos,
                ref_pos=v.ref_pos,
                ref_base=v.ref_base,
                alt_base=v.alt_base,
                variant_type=v.variant_type,
                consequence=v.consequence,
                aa_ref=aa_ref,
                aa_alt=aa_alt,
                amino_acid_change=amino_acid_change,
                risk_level=risk_level,
                risk_reason=risk_reason,
                carrier_samples=[],
            )

        sample = sample_map.get(v.sample_id)
        if sample:
            already_present = any(
                cs.sample_id == sample.id
                for cs in variant_groups[variant_key].carrier_samples
            )
            if not already_present:
                variant_groups[variant_key].carrier_samples.append(
                    schemas.DomainVariantCarrier(
                        sample_id=sample.id,
                        sample_name=sample.name,
                    )
                )

    variants_list = list(variant_groups.values())
    variants_list.sort(key=lambda r: RISK_ORDER.get(r.risk_level, 0), reverse=True)

    all_sample_ids: set = set()
    high = medium = low = unknown = 0
    for v in variants_list:
        for cs in v.carrier_samples:
            all_sample_ids.add(cs.sample_id)
        if v.risk_level == "高风险":
            high += 1
        elif v.risk_level == "中风险":
            medium += 1
        elif v.risk_level == "低风险":
            low += 1
        else:
            unknown += 1

    summary_stats = schemas.DomainSummaryStats(
        total_variants=len(variants_list),
        unique_samples=len(all_sample_ids),
        high_risk_count=high,
        medium_risk_count=medium,
        low_risk_count=low,
        unknown_risk_count=unknown,
    )

    return schemas.DomainVariantSummaryOut(
        gene_name=gene_name,
        domain_name=domain_name,
        domain_type=domain.domain_type,
        start_codon=domain.start_codon,
        end_codon=domain.end_codon,
        summary_stats=summary_stats,
        variants=variants_list,
    )
