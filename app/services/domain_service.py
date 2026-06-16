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


def _compute_codon_position(ref_pos: int, gene_start: int) -> int:
    return (ref_pos - gene_start) // 3 + 1


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

    gene_annotations: Dict[str, models.GeneAnnotation] = {}
    if gene_names:
        ann_rows = (
            db.query(models.GeneAnnotation)
            .filter(models.GeneAnnotation.gene_name.in_(gene_names))
            .all()
        )
        for ann in ann_rows:
            if ann.gene_name not in gene_annotations:
                gene_annotations[ann.gene_name] = ann

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

        ann = gene_annotations.get(gene_name)
        if not ann:
            continue

        codon_pos = _compute_codon_position(v.ref_pos, ann.start)

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
