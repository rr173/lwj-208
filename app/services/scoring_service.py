from typing import List, Optional, Dict, Any
from sqlalchemy.orm import Session
from sqlalchemy import func
from app import models, schemas


DEFAULT_RULES = [
    {
        "name": "nonsense_mutation",
        "description": "stop_gained (nonsense) variant",
        "condition_field": "consequence",
        "condition_operator": "eq",
        "condition_value": "stop_gained",
        "weight": 8.0,
    },
    {
        "name": "missense_mutation",
        "description": "Missense variant",
        "condition_field": "consequence",
        "condition_operator": "eq",
        "condition_value": "missense_variant",
        "weight": 4.0,
    },
    {
        "name": "synonymous_mutation",
        "description": "Synonymous variant",
        "condition_field": "consequence",
        "condition_operator": "eq",
        "condition_value": "synonymous_variant",
        "weight": -3.0,
    },
    {
        "name": "intron_variant",
        "description": "Intronic variant",
        "condition_field": "feature_type",
        "condition_operator": "eq",
        "condition_value": "intron",
        "weight": -2.0,
    },
    {
        "name": "intergenic_variant",
        "description": "Intergenic variant",
        "condition_field": "feature_type",
        "condition_operator": "eq",
        "condition_value": "intergenic",
        "weight": -4.0,
    },
    {
        "name": "high_population_frequency",
        "description": "Population frequency > 0.05",
        "condition_field": "population_frequency",
        "condition_operator": "gt",
        "condition_value": "0.05",
        "weight": -6.0,
    },
]


def _score_to_classification(score: float) -> str:
    if score >= 6:
        return "pathogenic"
    if score >= 3:
        return "likely_pathogenic"
    if score >= -2:
        return "uncertain_significance"
    if score > -6:
        return "likely_benign"
    return "benign"


def _get_or_create_rule_version(db: Session) -> models.RuleSetVersion:
    rsv = db.query(models.RuleSetVersion).first()
    if not rsv:
        rsv = models.RuleSetVersion(version=1)
        db.add(rsv)
        db.commit()
        db.refresh(rsv)
    return rsv


def get_current_rule_version(db: Session) -> int:
    rsv = _get_or_create_rule_version(db)
    return rsv.version


def _bump_rule_version(db: Session) -> int:
    rsv = _get_or_create_rule_version(db)
    rsv.version += 1
    db.commit()
    db.refresh(rsv)
    return rsv.version


def seed_default_rules(db: Session) -> None:
    existing = db.query(models.ScoringRule).first()
    if existing:
        return
    for rule_data in DEFAULT_RULES:
        rule = models.ScoringRule(**rule_data, enabled=True)
        db.add(rule)
    _get_or_create_rule_version(db)
    db.commit()


def create_rule(db: Session, data: schemas.ScoringRuleCreate) -> models.ScoringRule:
    rule = models.ScoringRule(
        name=data.name,
        description=data.description or "",
        condition_field=data.condition_field,
        condition_operator=data.condition_operator,
        condition_value=data.condition_value,
        weight=data.weight,
        enabled=data.enabled,
    )
    db.add(rule)
    _bump_rule_version(db)
    db.commit()
    db.refresh(rule)
    return rule


def get_rule(db: Session, rule_id: int) -> Optional[models.ScoringRule]:
    return db.query(models.ScoringRule).filter(models.ScoringRule.id == rule_id).first()


def list_rules(db: Session, enabled_only: bool = False) -> List[models.ScoringRule]:
    q = db.query(models.ScoringRule).order_by(models.ScoringRule.id)
    if enabled_only:
        q = q.filter(models.ScoringRule.enabled == True)
    return q.all()


def update_rule(
    db: Session, rule: models.ScoringRule, data: schemas.ScoringRuleUpdate
) -> models.ScoringRule:
    update_data = data.dict(exclude_unset=True)
    changed = False
    for field, value in update_data.items():
        if getattr(rule, field) != value:
            setattr(rule, field, value)
            changed = True
    if changed:
        _bump_rule_version(db)
    db.commit()
    db.refresh(rule)
    return rule


def delete_rule(db: Session, rule: models.ScoringRule) -> None:
    db.delete(rule)
    _bump_rule_version(db)
    db.commit()


def _evaluate_condition(
    operator: str, actual_value: Any, condition_value: str
) -> bool:
    try:
        if actual_value is None:
            return False
        cv_typed: Any
        if isinstance(actual_value, (int, float)):
            cv_typed = float(condition_value)
        else:
            cv_typed = condition_value
            actual_value = str(actual_value)
        if operator == "eq":
            return actual_value == cv_typed
        if operator == "ne":
            return actual_value != cv_typed
        if operator == "gt":
            return float(actual_value) > float(cv_typed)
        if operator == "gte":
            return float(actual_value) >= float(cv_typed)
        if operator == "lt":
            return float(actual_value) < float(cv_typed)
        if operator == "lte":
            return float(actual_value) <= float(cv_typed)
    except (ValueError, TypeError):
        return False
    return False


def _get_variant_population_frequency(
    db: Session, reference_id: int, ref_pos: int,
    variant_type: str, ref_base: str, alt_base: str,
) -> Optional[float]:
    db.commit()
    total_samples = db.query(func.count(models.Sample.id)).scalar() or 0
    if total_samples == 0:
        return None
    sample_count = db.query(
        func.count(func.distinct(models.SampleVariantSpectrum.sample_id))
    ).filter(
        models.SampleVariantSpectrum.reference_id == reference_id,
        models.SampleVariantSpectrum.ref_pos == ref_pos,
        models.SampleVariantSpectrum.variant_type == variant_type,
        models.SampleVariantSpectrum.ref_base == ref_base,
        models.SampleVariantSpectrum.alt_base == alt_base,
    ).scalar() or 0
    return sample_count / total_samples


def _score_variant(
    db: Session,
    variant: models.SampleVariantSpectrum,
    rules: List[models.ScoringRule],
) -> Dict[str, Any]:
    variant_attrs: Dict[str, Any] = {
        "impact": variant.impact,
        "feature_type": variant.feature_type,
        "variant_type": variant.variant_type,
        "consequence": variant.consequence,
    }

    pop_freq = _get_variant_population_frequency(
        db,
        variant.reference_id,
        variant.ref_pos,
        variant.variant_type,
        variant.ref_base,
        variant.alt_base,
    )
    variant_attrs["population_frequency"] = pop_freq

    total = 0.0
    matched = []

    for rule in rules:
        field = rule.condition_field
        actual_value = variant_attrs.get(field)
        if field == "population_frequency" and actual_value is None:
            continue
        if _evaluate_condition(rule.condition_operator, actual_value, rule.condition_value):
            total += rule.weight
            matched.append({
                "rule_id": rule.id,
                "rule_name": rule.name,
                "weight": rule.weight,
            })

    return {
        "reference_id": variant.reference_id,
        "ref_pos": variant.ref_pos,
        "variant_type": variant.variant_type,
        "ref_base": variant.ref_base,
        "alt_base": variant.alt_base,
        "feature_type": variant.feature_type,
        "gene_name": variant.gene_name,
        "impact": variant.impact,
        "consequence": variant.consequence,
        "population_frequency": pop_freq,
        "matched_rules": matched,
        "variant_score": total,
        "variant_classification": _score_to_classification(total),
    }


def score_sample(db: Session, sample_id: int) -> schemas.SampleScoreOut:
    sample = db.query(models.Sample).filter(models.Sample.id == sample_id).first()
    if not sample:
        raise ValueError(f"Sample {sample_id} not found")

    spectrum = db.query(models.SampleVariantSpectrum).filter(
        models.SampleVariantSpectrum.sample_id == sample_id
    ).order_by(
        models.SampleVariantSpectrum.reference_id,
        models.SampleVariantSpectrum.ref_pos,
    ).all()

    rules = list_rules(db, enabled_only=True)
    current_version = get_current_rule_version(db)

    variant_details = []
    total_score = 0.0

    for v in spectrum:
        vd = _score_variant(db, v, rules)
        variant_details.append(vd)
        total_score += vd["variant_score"]

    classification = _score_to_classification(total_score)

    existing = db.query(models.SamplePathogenicityScore).filter(
        models.SamplePathogenicityScore.sample_id == sample_id
    ).first()

    if existing:
        existing.total_score = total_score
        existing.classification = classification
        existing.rule_version = current_version
        existing.variant_scores = variant_details
    else:
        existing = models.SamplePathogenicityScore(
            sample_id=sample_id,
            total_score=total_score,
            classification=classification,
            rule_version=current_version,
            variant_scores=variant_details,
        )
        db.add(existing)

    db.commit()
    db.refresh(existing)

    return schemas.SampleScoreOut(
        sample_id=sample_id,
        sample_name=sample.name,
        total_score=total_score,
        classification=classification,
        rule_version=existing.rule_version,
        current_rule_version=current_version,
        score_may_be_outdated=existing.rule_version < current_version,
        total_variants=len(variant_details),
        variant_details=[schemas.VariantScoreDetail(**vd) for vd in variant_details],
        scored_at=existing.scored_at,
    )


def get_sample_score(db: Session, sample_id: int) -> Optional[schemas.SampleScoreSummaryOut]:
    existing = db.query(models.SamplePathogenicityScore).filter(
        models.SamplePathogenicityScore.sample_id == sample_id
    ).first()
    if not existing:
        return None

    sample = db.query(models.Sample).filter(models.Sample.id == sample_id).first()
    current_version = get_current_rule_version(db)

    return schemas.SampleScoreSummaryOut(
        sample_id=sample_id,
        sample_name=sample.name if sample else f"sample_{sample_id}",
        total_score=existing.total_score,
        classification=existing.classification,
        rule_version=existing.rule_version,
        current_rule_version=current_version,
        score_may_be_outdated=existing.rule_version < current_version,
        total_variants=len(existing.variant_scores) if existing.variant_scores else 0,
        scored_at=existing.scored_at,
    )


def get_sample_score_detail(db: Session, sample_id: int) -> Optional[schemas.SampleScoreOut]:
    existing = db.query(models.SamplePathogenicityScore).filter(
        models.SamplePathogenicityScore.sample_id == sample_id
    ).first()
    if not existing:
        return None

    sample = db.query(models.Sample).filter(models.Sample.id == sample_id).first()
    current_version = get_current_rule_version(db)

    variant_details = []
    if existing.variant_scores:
        variant_details = [schemas.VariantScoreDetail(**vd) for vd in existing.variant_scores]

    return schemas.SampleScoreOut(
        sample_id=sample_id,
        sample_name=sample.name if sample else f"sample_{sample_id}",
        total_score=existing.total_score,
        classification=existing.classification,
        rule_version=existing.rule_version,
        current_rule_version=current_version,
        score_may_be_outdated=existing.rule_version < current_version,
        total_variants=len(variant_details),
        variant_details=variant_details,
        scored_at=existing.scored_at,
    )


def filter_variants_by_classification(
    db: Session,
    sample_id: int,
    classification: Optional[str] = None,
    gene_name: Optional[str] = None,
) -> List[schemas.FilteredVariantOut]:
    existing = db.query(models.SamplePathogenicityScore).filter(
        models.SamplePathogenicityScore.sample_id == sample_id
    ).first()
    if not existing or not existing.variant_scores:
        return []

    results = []
    for vd in existing.variant_scores:
        if classification and vd.get("variant_classification") != classification:
            continue
        if gene_name and vd.get("gene_name") != gene_name:
            continue
        results.append(schemas.FilteredVariantOut(**vd))

    return results


def get_score_stats(db: Session, sample_id: int) -> Optional[schemas.ScoreStatsOut]:
    existing = db.query(models.SamplePathogenicityScore).filter(
        models.SamplePathogenicityScore.sample_id == sample_id
    ).first()
    if not existing:
        return None

    sample = db.query(models.Sample).filter(models.Sample.id == sample_id).first()

    distribution: Dict[str, int] = {c: 0 for c in schemas.VALID_CLASSIFICATIONS}
    if existing.variant_scores:
        for vd in existing.variant_scores:
            cls = vd.get("variant_classification", "uncertain_significance")
            if cls in distribution:
                distribution[cls] += 1

    return schemas.ScoreStatsOut(
        sample_id=sample_id,
        sample_name=sample.name if sample else f"sample_{sample_id}",
        total_variants=len(existing.variant_scores) if existing.variant_scores else 0,
        classification_distribution=[
            schemas.ClassificationCount(classification=cls, count=count)
            for cls, count in distribution.items()
        ],
    )
