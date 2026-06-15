from typing import List, Optional, Dict, Tuple, Set
from sqlalchemy.orm import Session
from sqlalchemy import func, and_, or_
from datetime import datetime
from app import models, schemas
from app.services.phylogeny_service import invalidate_tree_results_for_samples


def _variant_five_tuple_key(v) -> Tuple:
    return (
        v.reference_id,
        v.ref_pos,
        v.variant_type,
        v.ref_base,
        v.alt_base,
    )


def create_sample(db: Session, data: schemas.SampleCreate) -> models.Sample:
    sample = models.Sample(
        name=data.name,
        species=data.species,
        collection_date=data.collection_date,
        notes=data.notes or "",
    )
    db.add(sample)
    db.commit()
    db.refresh(sample)
    return sample


def get_sample_by_id(db: Session, sample_id: int) -> Optional[models.Sample]:
    return db.query(models.Sample).filter(models.Sample.id == sample_id).first()


def get_sample_by_name(db: Session, name: str) -> Optional[models.Sample]:
    return db.query(models.Sample).filter(models.Sample.name == name).first()


def list_samples(db: Session, skip: int = 0, limit: int = 100) -> List[models.Sample]:
    return db.query(models.Sample).order_by(models.Sample.created_at.desc()).offset(skip).limit(limit).all()


def update_sample(db: Session, sample: models.Sample, data: schemas.SampleUpdate) -> models.Sample:
    update_data = data.dict(exclude_unset=True)
    for field, value in update_data.items():
        setattr(sample, field, value)
    db.commit()
    db.refresh(sample)
    return sample


def delete_sample(db: Session, sample: models.Sample) -> None:
    sample_id = sample.id
    db.delete(sample)
    db.commit()
    invalidate_tree_results_for_samples(db, [sample_id])


def sample_to_out(db: Session, sample: models.Sample, with_detail: bool = False) -> schemas.SampleOut:
    align_count = len(sample.alignment_links)
    variant_count = db.query(func.count(models.SampleVariantSpectrum.id)).filter(
        models.SampleVariantSpectrum.sample_id == sample.id
    ).scalar() or 0

    if with_detail:
        linked_ids = [link.alignment_id for link in sample.alignment_links]
        return schemas.SampleDetail(
            id=sample.id,
            name=sample.name,
            species=sample.species,
            collection_date=sample.collection_date,
            notes=sample.notes,
            created_at=sample.created_at,
            updated_at=sample.updated_at,
            alignment_count=align_count,
            variant_count=variant_count,
            linked_alignment_ids=linked_ids,
        )
    return schemas.SampleOut(
        id=sample.id,
        name=sample.name,
        species=sample.species,
        collection_date=sample.collection_date,
        notes=sample.notes,
        created_at=sample.created_at,
        updated_at=sample.updated_at,
        alignment_count=align_count,
        variant_count=variant_count,
    )


def _collect_variants_from_alignments(db: Session, alignment_ids: List[int]) -> Dict[Tuple, dict]:
    """Collect unique variants from given alignment IDs, keyed by five-tuple."""
    variants = db.query(models.Variant).filter(
        models.Variant.alignment_id.in_(alignment_ids)
    ).all()

    result: Dict[Tuple, dict] = {}
    for v in variants:
        key = _variant_five_tuple_key(v)
        if key not in result:
            result[key] = {
                "reference_id": v.reference_id,
                "ref_pos": v.ref_pos,
                "variant_type": v.variant_type,
                "ref_base": v.ref_base,
                "alt_base": v.alt_base,
                "feature_type": v.feature_type,
                "gene_name": v.gene_name,
                "impact": v.impact,
                "consequence": v.consequence,
                "aa_ref": v.aa_ref,
                "aa_alt": v.aa_alt,
                "source_alignment_ids": set(),
            }
        result[key]["source_alignment_ids"].add(v.alignment_id)
    return result


def rebuild_sample_spectrum(db: Session, sample_id: int) -> None:
    """Rebuild the variant spectrum cache for a sample from scratch."""
    db.query(models.SampleVariantSpectrum).filter(
        models.SampleVariantSpectrum.sample_id == sample_id
    ).delete(synchronize_session=False)

    alignment_ids = [
        link.alignment_id
        for link in db.query(models.SampleAlignmentLink).filter(
            models.SampleAlignmentLink.sample_id == sample_id
        ).all()
    ]

    if not alignment_ids:
        db.commit()
        _invalidate_frequency_cache_for_sample(db, sample_id)
        return

    variants_by_key = _collect_variants_from_alignments(db, alignment_ids)

    for key, vdata in variants_by_key.items():
        spec = models.SampleVariantSpectrum(
            sample_id=sample_id,
            reference_id=vdata["reference_id"],
            ref_pos=vdata["ref_pos"],
            variant_type=vdata["variant_type"],
            ref_base=vdata["ref_base"],
            alt_base=vdata["alt_base"],
            feature_type=vdata["feature_type"],
            gene_name=vdata["gene_name"],
            impact=vdata["impact"],
            consequence=vdata["consequence"],
            aa_ref=vdata["aa_ref"],
            aa_alt=vdata["aa_alt"],
            source_alignment_ids=sorted(list(vdata["source_alignment_ids"])),
        )
        db.add(spec)

    db.commit()
    _invalidate_frequency_cache_for_sample(db, sample_id)
    invalidate_tree_results_for_samples(db, [sample_id])


def link_alignments_to_sample(
    db: Session, sample: models.Sample, alignment_ids: List[int]
) -> Tuple[int, List[int]]:
    """
    Link alignment results to a sample. Returns (added_count, invalid_alignment_ids).
    Invalid IDs are those that don't exist.
    """
    sample_id = sample.id
    sample_db = db.query(models.Sample).filter(models.Sample.id == sample_id).first()
    if not sample_db:
        return 0, list(alignment_ids)

    existing_ids = {link.alignment_id for link in sample_db.alignment_links}
    invalid_ids = []
    to_add = []

    for aid in alignment_ids:
        alignment = db.query(models.AlignmentResult).filter(
            models.AlignmentResult.id == aid
        ).first()
        if not alignment:
            invalid_ids.append(aid)
            continue
        if aid in existing_ids:
            continue
        to_add.append(aid)

    for aid in to_add:
        link = models.SampleAlignmentLink(
            sample_id=sample_id,
            alignment_id=aid,
        )
        db.add(link)

    if to_add or invalid_ids:
        db.commit()

    if to_add:
        rebuild_sample_spectrum(db, sample_id)

    return len(to_add), invalid_ids


def unlink_alignment_from_sample(
    db: Session, sample: models.Sample, alignment_id: int
) -> bool:
    sample_id = sample.id
    link = db.query(models.SampleAlignmentLink).filter(
        models.SampleAlignmentLink.sample_id == sample_id,
        models.SampleAlignmentLink.alignment_id == alignment_id,
    ).first()
    if not link:
        return False
    db.delete(link)
    db.commit()
    rebuild_sample_spectrum(db, sample_id)
    return True


def get_sample_spectrum(
    db: Session,
    sample: models.Sample,
    reference_name: Optional[str] = None,
    impact_filter: Optional[str] = None,
    gene_filter: Optional[str] = None,
    variant_type_filter: Optional[str] = None,
) -> schemas.SampleSpectrumOut:
    sample_id = sample.id
    query = db.query(models.SampleVariantSpectrum).filter(
        models.SampleVariantSpectrum.sample_id == sample_id
    )

    if reference_name:
        ref = db.query(models.ReferenceSequence).filter(
            models.ReferenceSequence.name == reference_name
        ).first()
        if ref:
            query = query.filter(models.SampleVariantSpectrum.reference_id == ref.id)

    if impact_filter:
        query = query.filter(models.SampleVariantSpectrum.impact == impact_filter)
    if gene_filter:
        query = query.filter(models.SampleVariantSpectrum.gene_name == gene_filter)
    if variant_type_filter:
        query = query.filter(models.SampleVariantSpectrum.variant_type == variant_type_filter)

    spectrum_rows = query.order_by(
        models.SampleVariantSpectrum.reference_id,
        models.SampleVariantSpectrum.ref_pos,
    ).all()

    ref_cache: Dict[int, str] = {}

    def get_ref_name(ref_id: int) -> str:
        if ref_id not in ref_cache:
            ref = db.query(models.ReferenceSequence).filter(
                models.ReferenceSequence.id == ref_id
            ).first()
            ref_cache[ref_id] = ref.name if ref else f"ref_{ref_id}"
        return ref_cache[ref_id]

    variants_out = [
        schemas.SpectrumVariantOut(
            reference_name=get_ref_name(r.reference_id),
            ref_pos=r.ref_pos,
            variant_type=r.variant_type,
            ref_base=r.ref_base,
            alt_base=r.alt_base,
            feature_type=r.feature_type,
            gene_name=r.gene_name,
            impact=r.impact,
            consequence=r.consequence,
            aa_ref=r.aa_ref,
            aa_alt=r.aa_alt,
            source_alignment_ids=r.source_alignment_ids or [],
        )
        for r in spectrum_rows
    ]

    sample_db = db.query(models.Sample).filter(models.Sample.id == sample_id).first()
    return schemas.SampleSpectrumOut(
        sample_id=sample_id,
        sample_name=sample_db.name if sample_db else f"sample_{sample_id}",
        total_variants=len(variants_out),
        variants=variants_out,
    )


def _get_affected_reference_ids_for_sample(db: Session, sample_id: int) -> Set[int]:
    rows = db.query(models.SampleVariantSpectrum.reference_id).filter(
        models.SampleVariantSpectrum.sample_id == sample_id
    ).distinct().all()
    return {r[0] for r in rows}


def _invalidate_frequency_cache_for_sample(db: Session, sample_id: int) -> None:
    ref_ids = _get_affected_reference_ids_for_sample(db, sample_id)
    for ref_id in ref_ids:
        meta = db.query(models.PopulationFrequencyMeta).filter(
            models.PopulationFrequencyMeta.reference_id == ref_id
        ).first()
        if meta:
            db.delete(meta)
    db.query(models.PopulationFrequencyCache).filter(
        models.PopulationFrequencyCache.reference_id.in_(list(ref_ids))
    ).delete(synchronize_session=False)
    db.commit()


def _is_frequency_cache_valid(db: Session, reference_id: int, total_samples: int) -> bool:
    meta = db.query(models.PopulationFrequencyMeta).filter(
        models.PopulationFrequencyMeta.reference_id == reference_id
    ).first()
    if not meta:
        return False
    return meta.total_samples_analyzed == total_samples


def _build_frequency_cache(db: Session, reference_id: int) -> None:
    """Build frequency cache for a given reference across all samples."""
    db.query(models.PopulationFrequencyCache).filter(
        models.PopulationFrequencyCache.reference_id == reference_id
    ).delete(synchronize_session=False)

    sample_ids = [
        s.id for s in db.query(models.Sample).all()
    ]
    total_samples = len(sample_ids)

    variant_sample_count: Dict[Tuple, Dict] = {}

    spectrum_rows = db.query(models.SampleVariantSpectrum).filter(
        models.SampleVariantSpectrum.reference_id == reference_id
    ).all()

    for row in spectrum_rows:
        key = (
            row.reference_id,
            row.ref_pos,
            row.variant_type,
            row.ref_base,
            row.alt_base,
        )
        if key not in variant_sample_count:
            variant_sample_count[key] = {
                "sample_ids": set(),
                "impact": row.impact,
                "gene_name": row.gene_name,
            }
        variant_sample_count[key]["sample_ids"].add(row.sample_id)

    for key, data in variant_sample_count.items():
        sample_count = len(data["sample_ids"])
        freq = sample_count / total_samples if total_samples > 0 else 0.0
        cache_entry = models.PopulationFrequencyCache(
            reference_id=key[0],
            ref_pos=key[1],
            variant_type=key[2],
            ref_base=key[3],
            alt_base=key[4],
            sample_count=sample_count,
            total_samples=total_samples,
            frequency=freq,
            impact=data["impact"],
            gene_name=data["gene_name"],
        )
        db.add(cache_entry)

    meta = db.query(models.PopulationFrequencyMeta).filter(
        models.PopulationFrequencyMeta.reference_id == reference_id
    ).first()
    if meta:
        meta.total_samples_analyzed = total_samples
        meta.last_updated = func.now()
    else:
        meta = models.PopulationFrequencyMeta(
            reference_id=reference_id,
            total_samples_analyzed=total_samples,
        )
        db.add(meta)

    db.commit()


def get_population_frequency(
    db: Session,
    reference_name: str,
    min_frequency: float = 0.0,
    max_frequency: float = 1.0,
) -> schemas.PopulationFrequencyOut:
    ref = db.query(models.ReferenceSequence).filter(
        models.ReferenceSequence.name == reference_name
    ).first()
    if not ref:
        return schemas.PopulationFrequencyOut(
            reference_name=reference_name,
            total_samples_analyzed=0,
            total_variant_sites=0,
            entries=[],
        )

    total_samples = db.query(func.count(models.Sample.id)).scalar() or 0

    if not _is_frequency_cache_valid(db, ref.id, total_samples):
        _build_frequency_cache(db, ref.id)

    query = db.query(models.PopulationFrequencyCache).filter(
        models.PopulationFrequencyCache.reference_id == ref.id,
        models.PopulationFrequencyCache.frequency >= min_frequency,
        models.PopulationFrequencyCache.frequency <= max_frequency,
    ).order_by(models.PopulationFrequencyCache.ref_pos)

    entries = [
        schemas.PopulationFrequencyEntry(
            ref_pos=r.ref_pos,
            variant_type=r.variant_type,
            ref_base=r.ref_base,
            alt_base=r.alt_base,
            frequency=r.frequency,
            sample_count=r.sample_count,
            total_samples=r.total_samples,
            impact=r.impact,
            gene_name=r.gene_name,
        )
        for r in query.all()
    ]

    return schemas.PopulationFrequencyOut(
        reference_name=reference_name,
        total_samples_analyzed=total_samples,
        total_variant_sites=len(entries),
        entries=entries,
    )


def _spectrum_to_dict(
    db: Session, sample: models.Sample
) -> Dict[Tuple, schemas.SpectrumVariantOut]:
    spectrum = get_sample_spectrum(db, sample)
    result = {}
    ref_cache: Dict[int, str] = {}
    for v in spectrum.variants:
        key = (v.reference_name, v.ref_pos, v.variant_type, v.ref_base, v.alt_base)
        result[key] = v
    return result


def _spectrum_by_position(
    spectrum: Dict[Tuple, schemas.SpectrumVariantOut]
) -> Dict[Tuple[str, int], List[schemas.SpectrumVariantOut]]:
    result = {}
    for key, variant in spectrum.items():
        ref_name, pos, *_ = key
        pos_key = (ref_name, pos)
        if pos_key not in result:
            result[pos_key] = []
        result[pos_key].append(variant)
    return result


def compare_samples(
    db: Session, sample_a: models.Sample, sample_b: models.Sample
) -> schemas.CompareResultOut:
    sample_a_db = db.query(models.Sample).filter(models.Sample.id == sample_a.id).first()
    sample_b_db = db.query(models.Sample).filter(models.Sample.id == sample_b.id).first()
    if not sample_a_db or not sample_b_db:
        return schemas.CompareResultOut(
            sample_a_id=sample_a.id,
            sample_a_name="unknown",
            sample_b_id=sample_b.id,
            sample_b_name="unknown",
            summary=schemas.CompareSummary(
                only_a_count=0, only_b_count=0, shared_count=0,
                same_site_diff_alt_count=0, jaccard_similarity=0.0,
            ),
        )

    spec_a = _spectrum_to_dict(db, sample_a_db)
    spec_b = _spectrum_to_dict(db, sample_b_db)

    keys_a = set(spec_a.keys())
    keys_b = set(spec_b.keys())

    only_a_keys = keys_a - keys_b
    only_b_keys = keys_b - keys_a
    shared_keys = keys_a & keys_b

    positions_a = _spectrum_by_position(spec_a)
    positions_b = _spectrum_by_position(spec_b)

    same_site_diff_alt = []
    all_positions = set(positions_a.keys()) | set(positions_b.keys())
    for pos_key in all_positions:
        vars_a = positions_a.get(pos_key, [])
        vars_b = positions_b.get(pos_key, [])
        if len(vars_a) > 0 and len(vars_b) > 0:
            alts_a = {(v.variant_type, v.ref_base, v.alt_base) for v in vars_a}
            alts_b = {(v.variant_type, v.ref_base, v.alt_base) for v in vars_b}
            if alts_a != alts_b and not (alts_a & alts_b):
                same_site_diff_alt.append({
                    "reference_name": pos_key[0],
                    "ref_pos": pos_key[1],
                    "sample_a_variants": [v.dict() for v in vars_a],
                    "sample_b_variants": [v.dict() for v in vars_b],
                })

    shared_exact = sorted([spec_a[k] for k in shared_keys], key=lambda v: (v.reference_name, v.ref_pos))
    only_a = sorted([spec_a[k] for k in only_a_keys], key=lambda v: (v.reference_name, v.ref_pos))
    only_b = sorted([spec_b[k] for k in only_b_keys], key=lambda v: (v.reference_name, v.ref_pos))

    union_count = len(keys_a | keys_b)
    jaccard = len(shared_keys) / union_count if union_count > 0 else 0.0

    summary = schemas.CompareSummary(
        only_a_count=len(only_a),
        only_b_count=len(only_b),
        shared_count=len(shared_exact),
        same_site_diff_alt_count=len(same_site_diff_alt),
        jaccard_similarity=round(jaccard, 6),
    )

    return schemas.CompareResultOut(
        sample_a_id=sample_a_db.id,
        sample_a_name=sample_a_db.name,
        sample_b_id=sample_b_db.id,
        sample_b_name=sample_b_db.name,
        summary=summary,
        only_in_a=only_a,
        only_in_b=only_b,
        shared=shared_exact,
        same_site_different_alt=same_site_diff_alt,
    )


def find_hotspots(
    db: Session,
    reference_name: str,
    threshold_per_100bp: float = 3.0,
) -> schemas.HotspotAnalysisOut:
    ref = db.query(models.ReferenceSequence).filter(
        models.ReferenceSequence.name == reference_name
    ).first()
    if not ref:
        return schemas.HotspotAnalysisOut(
            reference_name=reference_name,
            total_samples_analyzed=0,
            threshold_per_100bp=threshold_per_100bp,
            hotspot_count=0,
            hotspots=[],
        )

    total_samples = db.query(func.count(models.Sample.id)).scalar() or 0

    if not _is_frequency_cache_valid(db, ref.id, total_samples):
        _build_frequency_cache(db, ref.id)

    freq_entries = db.query(models.PopulationFrequencyCache).filter(
        models.PopulationFrequencyCache.reference_id == ref.id
    ).order_by(models.PopulationFrequencyCache.ref_pos).all()

    annotations = db.query(models.GeneAnnotation).filter(
        models.GeneAnnotation.reference_id == ref.id
    ).order_by(models.GeneAnnotation.start).all()

    freq_by_pos = {}
    for fe in freq_entries:
        pos_key = fe.ref_pos
        if pos_key not in freq_by_pos:
            freq_by_pos[pos_key] = []
        freq_by_pos[pos_key].append(fe)

    hotspots: List[schemas.HotspotEntry] = []

    for ann in annotations:
        region_start = ann.start
        region_end = ann.end
        region_length = region_end - region_start + 1
        if region_length <= 0:
            continue

        variant_positions_in_region = [
            pos for pos in freq_by_pos.keys()
            if region_start <= pos <= region_end
        ]
        variant_count = len(variant_positions_in_region)
        if variant_count == 0:
            continue

        density = (variant_count / region_length) * 100.0
        if density < threshold_per_100bp:
            continue

        all_variants_in_region = []
        sample_set = set()
        top_freq_entry = None
        top_freq = -1.0
        for pos in variant_positions_in_region:
            for fe in freq_by_pos[pos]:
                all_variants_in_region.append(fe)
                if fe.frequency > top_freq:
                    top_freq = fe.frequency
                    top_freq_entry = fe
                variant_samples = db.query(models.SampleVariantSpectrum.sample_id).filter(
                    models.SampleVariantSpectrum.reference_id == ref.id,
                    models.SampleVariantSpectrum.ref_pos == fe.ref_pos,
                    models.SampleVariantSpectrum.variant_type == fe.variant_type,
                    models.SampleVariantSpectrum.ref_base == fe.ref_base,
                    models.SampleVariantSpectrum.alt_base == fe.alt_base,
                ).all()
                for (sid,) in variant_samples:
                    sample_set.add(sid)

        top_variant = None
        if top_freq_entry:
            top_variant = schemas.PopulationFrequencyEntry(
                ref_pos=top_freq_entry.ref_pos,
                variant_type=top_freq_entry.variant_type,
                ref_base=top_freq_entry.ref_base,
                alt_base=top_freq_entry.alt_base,
                frequency=top_freq_entry.frequency,
                sample_count=top_freq_entry.sample_count,
                total_samples=top_freq_entry.total_samples,
                impact=top_freq_entry.impact,
                gene_name=top_freq_entry.gene_name,
            )

        hotspot = schemas.HotspotEntry(
            region_start=region_start,
            region_end=region_end,
            region_length=region_length,
            gene_name=ann.gene_name,
            feature_type=ann.feature_type,
            variant_count=variant_count,
            sample_count=len(sample_set),
            density_per_100bp=round(density, 4),
            top_variant=top_variant,
        )
        hotspots.append(hotspot)

    hotspots.sort(key=lambda h: h.density_per_100bp, reverse=True)

    return schemas.HotspotAnalysisOut(
        reference_name=reference_name,
        total_samples_analyzed=total_samples,
        threshold_per_100bp=threshold_per_100bp,
        hotspot_count=len(hotspots),
        hotspots=hotspots,
    )
