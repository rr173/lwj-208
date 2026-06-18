from typing import List, Optional, Dict, Any, Tuple, Set, NamedTuple
from sqlalchemy.orm import Session
from sqlalchemy import func
from datetime import datetime
from app import models, schemas


class _QCInternalResult(NamedTuple):
    sample: models.Sample
    active_rule_set: models.QCRuleSet
    rules: List[models.QCRule]
    spectrum_rows: List[models.SampleVariantSpectrum]
    pass_count: int
    fail_count: int
    qc_results: List[schemas.VariantQCOut]
    failed_rule_types: Dict[str, int]


def create_rule_set(db: Session, data: schemas.QCRuleSetCreate) -> models.QCRuleSet:
    if data.activate:
        _deactivate_all_rule_sets(db)

    rule_set = models.QCRuleSet(
        name=data.name,
        description=data.description or "",
        is_active=data.activate,
    )
    db.add(rule_set)
    db.flush()

    for idx, rule_data in enumerate(data.rules):
        rule = models.QCRule(
            rule_set_id=rule_set.id,
            rule_index=idx,
            rule_type=rule_data.rule_type,
            min_alignment_score=rule_data.min_alignment_score,
            min_coverage_depth=rule_data.min_coverage_depth,
            edge_margin=rule_data.edge_margin,
        )
        db.add(rule)

    if data.activate:
        _mark_existing_qc_results_stale(db)

    db.commit()
    db.refresh(rule_set)
    return rule_set


def get_rule_set(db: Session, rule_set_id: int) -> Optional[models.QCRuleSet]:
    return db.query(models.QCRuleSet).filter(models.QCRuleSet.id == rule_set_id).first()


def list_rule_sets(db: Session) -> List[models.QCRuleSet]:
    return db.query(models.QCRuleSet).order_by(models.QCRuleSet.created_at.desc()).all()


def update_rule_set(
    db: Session, rule_set: models.QCRuleSet, data: schemas.QCRuleSetUpdate
) -> models.QCRuleSet:
    update_data = data.dict(exclude_unset=True)
    rules_data = update_data.pop("rules", None)

    for field, value in update_data.items():
        setattr(rule_set, field, value)

    if rules_data is not None:
        db.query(models.QCRule).filter(
            models.QCRule.rule_set_id == rule_set.id
        ).delete(synchronize_session=False)

        for idx, rule_data in enumerate(rules_data):
            rule = models.QCRule(
                rule_set_id=rule_set.id,
                rule_index=idx,
                rule_type=rule_data.rule_type,
                min_alignment_score=rule_data.min_alignment_score,
                min_coverage_depth=rule_data.min_coverage_depth,
                edge_margin=rule_data.edge_margin,
            )
            db.add(rule)

        if rule_set.is_active:
            _mark_existing_qc_results_stale(db)

    rule_set.updated_at = func.now()
    db.commit()
    db.refresh(rule_set)
    return rule_set


def delete_rule_set(db: Session, rule_set: models.QCRuleSet) -> None:
    was_active = rule_set.is_active
    db.delete(rule_set)
    if was_active:
        _mark_existing_qc_results_stale(db)
    db.commit()


def activate_rule_set(db: Session, rule_set_id: int) -> models.QCRuleSet:
    rule_set = db.query(models.QCRuleSet).filter(models.QCRuleSet.id == rule_set_id).first()
    if not rule_set:
        raise ValueError(f"Rule set {rule_set_id} not found")

    _deactivate_all_rule_sets(db)
    rule_set.is_active = True
    _mark_existing_qc_results_stale(db)
    db.commit()
    db.refresh(rule_set)
    return rule_set


def deactivate_rule_set(db: Session, rule_set_id: int) -> models.QCRuleSet:
    rule_set = db.query(models.QCRuleSet).filter(models.QCRuleSet.id == rule_set_id).first()
    if not rule_set:
        raise ValueError(f"Rule set {rule_set_id} not found")

    rule_set.is_active = False
    _mark_existing_qc_results_stale(db)
    db.commit()
    db.refresh(rule_set)
    return rule_set


def get_active_rule_set(db: Session) -> Optional[models.QCRuleSet]:
    return db.query(models.QCRuleSet).filter(models.QCRuleSet.is_active == True).first()


def _deactivate_all_rule_sets(db: Session) -> None:
    db.query(models.QCRuleSet).filter(
        models.QCRuleSet.is_active == True
    ).update({"is_active": False})


def _mark_existing_qc_results_stale(db: Session) -> None:
    db.query(models.VariantQCResult).filter(
        models.VariantQCResult.is_stale == False
    ).update({"is_stale": True})


def _get_coverage_depth_for_variant(
    db: Session,
    sample_id: int,
    reference_id: int,
    ref_pos: int,
    variant_type: str,
    ref_base: str,
    alt_base: str,
) -> int:
    count = db.query(func.count(func.distinct(models.Variant.alignment_id))).join(
        models.SampleAlignmentLink,
        models.SampleAlignmentLink.alignment_id == models.Variant.alignment_id
    ).filter(
        models.SampleAlignmentLink.sample_id == sample_id,
        models.Variant.reference_id == reference_id,
        models.Variant.ref_pos == ref_pos,
        models.Variant.variant_type == variant_type,
        models.Variant.ref_base == ref_base,
        models.Variant.alt_base == alt_base,
    ).scalar() or 0
    return count


def _get_best_alignment_score_for_variant(
    db: Session,
    sample_id: int,
    reference_id: int,
    ref_pos: int,
    variant_type: str,
    ref_base: str,
    alt_base: str,
) -> Optional[int]:
    max_score = db.query(func.max(models.AlignmentResult.score)).join(
        models.Variant,
        models.Variant.alignment_id == models.AlignmentResult.id
    ).join(
        models.SampleAlignmentLink,
        models.SampleAlignmentLink.alignment_id == models.AlignmentResult.id
    ).filter(
        models.SampleAlignmentLink.sample_id == sample_id,
        models.Variant.reference_id == reference_id,
        models.Variant.ref_pos == ref_pos,
        models.Variant.variant_type == variant_type,
        models.Variant.ref_base == ref_base,
        models.Variant.alt_base == alt_base,
    ).scalar()

    return max_score


def _get_variant_position_relative_to_ends(
    db: Session, reference_id: int, ref_pos: int
) -> Tuple[int, int]:
    ref = db.query(models.ReferenceSequence).filter(
        models.ReferenceSequence.id == reference_id
    ).first()
    if not ref:
        return ref_pos, ref_pos

    dist_from_start = ref_pos
    dist_from_end = ref.length - ref_pos + 1
    return dist_from_start, dist_from_end


def _evaluate_variant_against_rules(
    db: Session,
    sample_id: int,
    reference_id: int,
    ref_pos: int,
    variant_type: str,
    ref_base: str,
    alt_base: str,
    rules: List[models.QCRule],
) -> Tuple[bool, List[Dict[str, Any]]]:
    failed_rules = []

    for rule in rules:
        if rule.rule_type == "alignment_score":
            best_score = _get_best_alignment_score_for_variant(
                db, sample_id, reference_id, ref_pos, variant_type, ref_base, alt_base
            )
            if best_score is None or best_score < rule.min_alignment_score:
                failed_rules.append({
                    "rule_id": rule.id,
                    "rule_type": rule.rule_type,
                    "rule_description": f"Alignment score {best_score} below threshold {rule.min_alignment_score}",
                    "actual_value": best_score,
                    "threshold": rule.min_alignment_score,
                })

        elif rule.rule_type == "coverage_depth":
            depth = _get_coverage_depth_for_variant(
                db, sample_id, reference_id, ref_pos, variant_type, ref_base, alt_base
            )
            if depth < rule.min_coverage_depth:
                failed_rules.append({
                    "rule_id": rule.id,
                    "rule_type": rule.rule_type,
                    "rule_description": f"Coverage depth {depth} below threshold {rule.min_coverage_depth}",
                    "actual_value": depth,
                    "threshold": rule.min_coverage_depth,
                })

        elif rule.rule_type == "edge_position":
            dist_start, dist_end = _get_variant_position_relative_to_ends(
                db, reference_id, ref_pos
            )
            min_dist = min(dist_start, dist_end)
            if min_dist <= rule.edge_margin:
                failed_rules.append({
                    "rule_id": rule.id,
                    "rule_type": rule.rule_type,
                    "rule_description": f"Position {ref_pos} is within {rule.edge_margin} bases of sequence end (min distance: {min_dist})",
                    "actual_value": min_dist,
                    "threshold": rule.edge_margin,
                })

    passed = len(failed_rules) == 0
    return passed, failed_rules


def _run_qc_internal(
    db: Session, sample_id: int
) -> _QCInternalResult:
    sample = db.query(models.Sample).filter(models.Sample.id == sample_id).first()
    if not sample:
        raise ValueError(f"Sample {sample_id} not found")

    active_rule_set = get_active_rule_set(db)
    if not active_rule_set:
        raise ValueError("No active QC rule set. Activate a rule set before running QC.")

    rules = db.query(models.QCRule).filter(
        models.QCRule.rule_set_id == active_rule_set.id
    ).order_by(models.QCRule.rule_index).all()

    if not rules:
        raise ValueError(f"Active rule set '{active_rule_set.name}' has no rules.")

    db.query(models.VariantQCResult).filter(
        models.VariantQCResult.sample_id == sample_id
    ).delete(synchronize_session=False)
    db.flush()

    spectrum_rows = db.query(models.SampleVariantSpectrum).filter(
        models.SampleVariantSpectrum.sample_id == sample_id
    ).order_by(
        models.SampleVariantSpectrum.reference_id,
        models.SampleVariantSpectrum.ref_pos,
    ).all()

    pass_count = 0
    fail_count = 0
    qc_results = []
    failed_rule_types: Dict[str, int] = {}

    for row in spectrum_rows:
        passed, failed_rules = _evaluate_variant_against_rules(
            db,
            sample_id,
            row.reference_id,
            row.ref_pos,
            row.variant_type,
            row.ref_base,
            row.alt_base,
            rules,
        )

        status = "PASS" if passed else "FAIL"
        if passed:
            pass_count += 1
        else:
            fail_count += 1
            seen_types = set()
            for fr in failed_rules:
                rt = fr.get("rule_type", "unknown")
                if rt not in seen_types:
                    seen_types.add(rt)
                    failed_rule_types[rt] = failed_rule_types.get(rt, 0) + 1

        qc_result = models.VariantQCResult(
            sample_id=sample_id,
            reference_id=row.reference_id,
            ref_pos=row.ref_pos,
            variant_type=row.variant_type,
            ref_base=row.ref_base,
            alt_base=row.alt_base,
            rule_set_id=active_rule_set.id,
            qc_status=status,
            failed_rules=failed_rules,
            is_stale=False,
        )
        db.add(qc_result)

        qc_results.append(schemas.VariantQCOut(
            reference_id=row.reference_id,
            ref_pos=row.ref_pos,
            variant_type=row.variant_type,
            ref_base=row.ref_base,
            alt_base=row.alt_base,
            qc_status=status,
            failed_rules=[schemas.FailedRuleDetail(**fr) for fr in failed_rules],
            is_stale=False,
        ))

    return _QCInternalResult(
        sample=sample,
        active_rule_set=active_rule_set,
        rules=rules,
        spectrum_rows=spectrum_rows,
        pass_count=pass_count,
        fail_count=fail_count,
        qc_results=qc_results,
        failed_rule_types=failed_rule_types,
    )


def execute_qc_for_sample(db: Session, sample_id: int) -> schemas.SampleQCResultOut:
    result = _run_qc_internal(db, sample_id)
    db.commit()
    return schemas.SampleQCResultOut(
        sample_id=sample_id,
        sample_name=result.sample.name,
        rule_set_id=result.active_rule_set.id,
        rule_set_name=result.active_rule_set.name,
        total_variants=len(result.spectrum_rows),
        pass_count=result.pass_count,
        fail_count=result.fail_count,
        stale_count=0,
        qc_results=result.qc_results,
        executed_at=datetime.now(),
    )


def get_sample_qc_summary(db: Session, sample_id: int) -> schemas.SampleQCSummaryOut:
    sample = db.query(models.Sample).filter(models.Sample.id == sample_id).first()
    if not sample:
        raise ValueError(f"Sample {sample_id} not found")

    total_variants = db.query(func.count(models.SampleVariantSpectrum.id)).filter(
        models.SampleVariantSpectrum.sample_id == sample_id
    ).scalar() or 0

    qc_results = db.query(models.VariantQCResult).filter(
        models.VariantQCResult.sample_id == sample_id
    ).all()

    has_qc = len(qc_results) > 0

    if not has_qc:
        return schemas.SampleQCSummaryOut(
            sample_id=sample_id,
            sample_name=sample.name,
            has_qc_results=False,
            is_stale=False,
            total_variants=total_variants,
            pass_count=0,
            fail_count=0,
            stale_count=0,
        )

    pass_count = sum(1 for r in qc_results if r.qc_status == "PASS" and not r.is_stale)
    fail_count = sum(1 for r in qc_results if r.qc_status == "FAIL" and not r.is_stale)
    stale_count = sum(1 for r in qc_results if r.is_stale)

    rule_set_id = qc_results[0].rule_set_id if qc_results else None
    rule_set_name = None
    is_stale = any(r.is_stale for r in qc_results)

    if rule_set_id:
        rs = db.query(models.QCRuleSet).filter(models.QCRuleSet.id == rule_set_id).first()
        if rs:
            rule_set_name = rs.name

    return schemas.SampleQCSummaryOut(
        sample_id=sample_id,
        sample_name=sample.name,
        rule_set_id=rule_set_id,
        rule_set_name=rule_set_name,
        has_qc_results=True,
        is_stale=is_stale,
        total_variants=total_variants,
        pass_count=pass_count,
        fail_count=fail_count,
        stale_count=stale_count,
    )


def get_sample_qc_results(
    db: Session, sample_id: int, qc_status: Optional[str] = None
) -> List[schemas.VariantQCOut]:
    sample = db.query(models.Sample).filter(models.Sample.id == sample_id).first()
    if not sample:
        raise ValueError(f"Sample {sample_id} not found")

    query = db.query(models.VariantQCResult).filter(
        models.VariantQCResult.sample_id == sample_id
    )

    if qc_status:
        query = query.filter(models.VariantQCResult.qc_status == qc_status)

    results = query.order_by(
        models.VariantQCResult.reference_id,
        models.VariantQCResult.ref_pos,
    ).all()

    return [
        schemas.VariantQCOut(
            reference_id=r.reference_id,
            ref_pos=r.ref_pos,
            variant_type=r.variant_type,
            ref_base=r.ref_base,
            alt_base=r.alt_base,
            qc_status=r.qc_status,
            failed_rules=[schemas.FailedRuleDetail(**fr) for fr in (r.failed_rules or [])],
            is_stale=r.is_stale,
        )
        for r in results
    ]


def get_pass_variant_keys_for_sample(
    db: Session, sample_id: int, include_failed_qc: bool = False
) -> Optional[Set[Tuple[int, int, str, str, str]]]:
    if include_failed_qc:
        return None

    active_rule_set = get_active_rule_set(db)
    if not active_rule_set:
        return None

    qc_results = db.query(models.VariantQCResult).filter(
        models.VariantQCResult.sample_id == sample_id,
        models.VariantQCResult.is_stale == False,
    ).all()

    if not qc_results:
        return None

    pass_keys = set()
    for r in qc_results:
        if r.qc_status == "PASS":
            pass_keys.add((r.reference_id, r.ref_pos, r.variant_type, r.ref_base, r.alt_base))
    return pass_keys


def filter_spectrum_by_qc(
    db: Session,
    spectrum_rows: List[models.SampleVariantSpectrum],
    sample_id: int,
    include_failed_qc: bool = False,
) -> List[models.SampleVariantSpectrum]:
    if include_failed_qc:
        return spectrum_rows

    pass_keys = get_pass_variant_keys_for_sample(db, sample_id, include_failed_qc)
    if pass_keys is None:
        return spectrum_rows

    filtered = []
    for row in spectrum_rows:
        key = (row.reference_id, row.ref_pos, row.variant_type, row.ref_base, row.alt_base)
        if key in pass_keys:
            filtered.append(row)
    return filtered


def execute_batch_qc(db: Session, sample_ids: List[int]) -> schemas.BatchQCResultOut:
    active_rule_set = get_active_rule_set(db)
    if not active_rule_set:
        raise ValueError("No active QC rule set. Activate a rule set before running QC.")

    rules = db.query(models.QCRule).filter(
        models.QCRule.rule_set_id == active_rule_set.id
    ).order_by(models.QCRule.rule_index).all()

    if not rules:
        raise ValueError(f"Active rule set '{active_rule_set.name}' has no rules.")

    sample_results: List[schemas.BatchQCSampleResult] = []
    rule_type_fail_counts: Dict[str, int] = {}
    processed_samples = 0
    failed_samples = 0
    total_pass_rate_sum = 0.0

    for sample_id in sample_ids:
        sample = db.query(models.Sample).filter(models.Sample.id == sample_id).first()
        if not sample:
            sample_results.append(schemas.BatchQCSampleResult(
                sample_id=sample_id,
                sample_name="NOT_FOUND",
                total_variants=0,
                pass_count=0,
                fail_count=0,
                pass_rate=0.0,
                error_message=f"Sample {sample_id} not found",
            ))
            failed_samples += 1
            continue

        try:
            savepoint = db.begin_nested()
            try:
                result = _run_qc_internal(db, sample_id)
                savepoint.commit()
            except Exception as inner_exc:
                savepoint.rollback()
                raise inner_exc

            for rt, count in result.failed_rule_types.items():
                rule_type_fail_counts[rt] = rule_type_fail_counts.get(rt, 0) + count

            total_variants = len(result.spectrum_rows)
            pass_rate = result.pass_count / total_variants if total_variants > 0 else 0.0

            sample_results.append(schemas.BatchQCSampleResult(
                sample_id=sample_id,
                sample_name=sample.name,
                total_variants=total_variants,
                pass_count=result.pass_count,
                fail_count=result.fail_count,
                pass_rate=round(pass_rate, 4),
                error_message=None,
            ))

            processed_samples += 1
            total_pass_rate_sum += pass_rate

        except Exception as e:
            failed_samples += 1
            sample_results.append(schemas.BatchQCSampleResult(
                sample_id=sample_id,
                sample_name=sample.name,
                total_variants=0,
                pass_count=0,
                fail_count=0,
                pass_rate=0.0,
                error_message=str(e),
            ))

    db.commit()

    average_pass_rate = total_pass_rate_sum / processed_samples if processed_samples > 0 else 0.0

    rule_fail_counts = [
        schemas.RuleFailCount(rule_type=rt, fail_count=fc)
        for rt, fc in sorted(rule_type_fail_counts.items(), key=lambda x: -x[1])
    ]

    return schemas.BatchQCResultOut(
        rule_set_id=active_rule_set.id,
        rule_set_name=active_rule_set.name,
        total_samples=len(sample_ids),
        processed_samples=processed_samples,
        failed_samples=failed_samples,
        average_pass_rate=round(average_pass_rate, 4),
        sample_results=sample_results,
        rule_fail_counts=rule_fail_counts,
        executed_at=datetime.now(),
    )


def generate_qc_report_by_reference(db: Session, reference_name: str) -> schemas.QCReportOut:
    ref = db.query(models.ReferenceSequence).filter(
        models.ReferenceSequence.name == reference_name
    ).first()
    if not ref:
        raise ValueError(f"Reference sequence '{reference_name}' not found")

    reference_length = ref.length

    all_results = db.query(models.VariantQCResult).filter(
        models.VariantQCResult.reference_id == ref.id,
        models.VariantQCResult.is_stale == False,
    ).all()

    if not all_results:
        return schemas.QCReportOut(
            reference_name=reference_name,
            reference_length=reference_length,
            total_samples=0,
            total_variants=0,
            total_pass=0,
            total_fail=0,
            overall_pass_rate=0.0,
            bucket_distribution=[],
            hotspot_regions=[],
            rule_type_stats=[],
            generated_at=datetime.now(),
        )

    total_variants = len(all_results)
    total_pass = sum(1 for r in all_results if r.qc_status == "PASS")
    total_fail = sum(1 for r in all_results if r.qc_status == "FAIL")
    overall_pass_rate = total_pass / total_variants if total_variants > 0 else 0.0

    sample_ids = set(r.sample_id for r in all_results)
    total_samples = len(sample_ids)

    BUCKET_SIZE = 500
    num_buckets = (reference_length + BUCKET_SIZE - 1) // BUCKET_SIZE
    bucket_data: List[Dict[str, int]] = [{"pass": 0, "fail": 0} for _ in range(num_buckets)]

    rule_type_counts: Dict[str, int] = {}
    total_rule_fail_events = 0

    for r in all_results:
        bucket_idx = min((r.ref_pos - 1) // BUCKET_SIZE, num_buckets - 1)
        if r.qc_status == "PASS":
            bucket_data[bucket_idx]["pass"] += 1
        else:
            bucket_data[bucket_idx]["fail"] += 1

        if r.qc_status == "FAIL" and r.failed_rules:
            rule_types_in_variant = set()
            for fr in r.failed_rules:
                rt = fr.get("rule_type")
                if rt and rt not in rule_types_in_variant:
                    rule_types_in_variant.add(rt)
            for rt in rule_types_in_variant:
                rule_type_counts[rt] = rule_type_counts.get(rt, 0) + 1
                total_rule_fail_events += 1

    bucket_distribution: List[schemas.QCBucketDistribution] = []
    for idx, bd in enumerate(bucket_data):
        start = idx * BUCKET_SIZE + 1
        end = min((idx + 1) * BUCKET_SIZE, reference_length)
        total = bd["pass"] + bd["fail"]
        fail_density = bd["fail"] / BUCKET_SIZE if BUCKET_SIZE > 0 else 0.0
        bucket_distribution.append(schemas.QCBucketDistribution(
            bucket_start=start,
            bucket_end=end,
            pass_count=bd["pass"],
            fail_count=bd["fail"],
            fail_density=round(fail_density, 6),
        ))

    sorted_buckets = sorted(
        bucket_distribution,
        key=lambda b: (b.fail_density, b.fail_count),
        reverse=True,
    )
    hotspot_regions: List[schemas.HotspotRegion] = []
    for rank, b in enumerate(sorted_buckets[:10], 1):
        if b.fail_count > 0:
            hotspot_regions.append(schemas.HotspotRegion(
                bucket_start=b.bucket_start,
                bucket_end=b.bucket_end,
                fail_count=b.fail_count,
                fail_density=b.fail_density,
                rank=rank,
            ))

    rule_type_stats: List[schemas.RuleTypeStats] = []
    for rt, count in sorted(rule_type_counts.items(), key=lambda x: -x[1]):
        pct = (count / total_rule_fail_events * 100) if total_rule_fail_events > 0 else 0.0
        rule_type_stats.append(schemas.RuleTypeStats(
            rule_type=rt,
            fail_count=count,
            fail_percentage=round(pct, 2),
        ))

    return schemas.QCReportOut(
        reference_name=reference_name,
        reference_length=reference_length,
        total_samples=total_samples,
        total_variants=total_variants,
        total_pass=total_pass,
        total_fail=total_fail,
        overall_pass_rate=round(overall_pass_rate, 4),
        bucket_distribution=bucket_distribution,
        hotspot_regions=hotspot_regions,
        rule_type_stats=rule_type_stats,
        generated_at=datetime.now(),
    )
