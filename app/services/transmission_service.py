from typing import List, Dict, Set, Tuple, Optional
from sqlalchemy.orm import Session
from sqlalchemy import func
from datetime import datetime

from app import models, schemas
from app.transmission import (
    jaccard_similarity,
    compute_jaccard_matrix,
    compute_distance_matrix_from_jaccard,
    upgma_cluster,
    transitive_reduction,
)


def _variant_five_tuple(v) -> Tuple:
    return (v.reference_id, v.ref_pos, v.variant_type, v.ref_base, v.alt_base)


def _ref_name_cache(db: Session) -> Dict[int, str]:
    refs = db.query(models.ReferenceSequence).all()
    return {r.id: r.name for r in refs}


def variant_tracing(
    db: Session,
    request: schemas.VariantTracingRequest,
) -> schemas.VariantTracingOut:
    ref = db.query(models.ReferenceSequence).filter(
        models.ReferenceSequence.name == request.reference_name
    ).first()
    if not ref:
        raise ValueError(f"Reference sequence '{request.reference_name}' not found")

    spectrum_rows = db.query(models.SampleVariantSpectrum).filter(
        models.SampleVariantSpectrum.reference_id == ref.id,
        models.SampleVariantSpectrum.ref_pos == request.ref_pos,
        models.SampleVariantSpectrum.variant_type == request.variant_type,
        models.SampleVariantSpectrum.ref_base == request.ref_base,
        models.SampleVariantSpectrum.alt_base == request.alt_base,
    ).all()

    if not spectrum_rows:
        return schemas.VariantTracingOut(
            reference_name=request.reference_name,
            ref_pos=request.ref_pos,
            variant_type=request.variant_type,
            ref_base=request.ref_base,
            alt_base=request.alt_base,
            total_carrier_samples=0,
            timeline=[],
        )

    sample_ids = list({r.sample_id for r in spectrum_rows})
    samples = db.query(models.Sample).filter(
        models.Sample.id.in_(sample_ids)
    ).all()
    sample_map: Dict[int, models.Sample] = {s.id: s for s in samples}

    entries: List[schemas.VariantTracingSampleEntry] = []
    with_date: List[schemas.VariantTracingSampleEntry] = []
    without_date: List[schemas.VariantTracingSampleEntry] = []

    for sid in sample_ids:
        sample = sample_map.get(sid)
        if not sample:
            continue
        entry = schemas.VariantTracingSampleEntry(
            sample_id=sample.id,
            sample_name=sample.name,
            collection_date=sample.collection_date,
            is_earliest=False,
            is_candidate_source=False,
            date_unknown=sample.collection_date is None,
        )
        if sample.collection_date is not None:
            with_date.append(entry)
        else:
            without_date.append(entry)

    with_date.sort(key=lambda e: e.collection_date)

    if with_date:
        earliest_date = with_date[0].collection_date
        for entry in with_date:
            if entry.collection_date == earliest_date:
                entry.is_earliest = True
                entry.is_candidate_source = True

    timeline = with_date + without_date

    return schemas.VariantTracingOut(
        reference_name=request.reference_name,
        ref_pos=request.ref_pos,
        variant_type=request.variant_type,
        ref_base=request.ref_base,
        alt_base=request.alt_base,
        total_carrier_samples=len(timeline),
        timeline=timeline,
    )


def _get_variant_sets_for_samples(
    db: Session,
    sample_ids: List[int],
    reference_id: Optional[int] = None,
) -> Tuple[Dict[int, Set[Tuple]], Dict[int, Dict[Tuple, models.SampleVariantSpectrum]]]:
    query = db.query(models.SampleVariantSpectrum).filter(
        models.SampleVariantSpectrum.sample_id.in_(sample_ids)
    )
    if reference_id is not None:
        query = query.filter(models.SampleVariantSpectrum.reference_id == reference_id)

    rows = query.all()

    variant_sets: Dict[int, Set[Tuple]] = {sid: set() for sid in sample_ids}
    variant_details: Dict[int, Dict[Tuple, models.SampleVariantSpectrum]] = {sid: {} for sid in sample_ids}

    for r in rows:
        key = (r.reference_id, r.ref_pos, r.variant_type, r.ref_base, r.alt_base)
        variant_sets[r.sample_id].add(key)
        variant_details[r.sample_id][key] = r

    return variant_sets, variant_details


def transmission_chain(
    db: Session,
    request: schemas.TransmissionChainRequest,
) -> schemas.TransmissionChainOut:
    sample_ids = list(set(request.sample_ids))

    for sid in sample_ids:
        sample = db.query(models.Sample).filter(models.Sample.id == sid).first()
        if not sample:
            raise ValueError(f"Sample {sid} not found")

    samples = db.query(models.Sample).filter(
        models.Sample.id.in_(sample_ids)
    ).all()
    sample_map: Dict[int, models.Sample] = {s.id: s for s in samples}

    ref_ids = db.query(models.SampleVariantSpectrum.reference_id).filter(
        models.SampleVariantSpectrum.sample_id.in_(sample_ids)
    ).distinct().all()
    ref_id_list = [r[0] for r in ref_ids]

    variant_sets: Dict[int, Set[Tuple]] = {sid: set() for sid in sample_ids}
    variant_details: Dict[int, Dict[Tuple, models.SampleVariantSpectrum]] = {sid: {} for sid in sample_ids}

    for rid in ref_id_list:
        vs, vd = _get_variant_sets_for_samples(db, sample_ids, rid)
        for sid in sample_ids:
            variant_sets[sid].update(vs[sid])
            variant_details[sid].update(vd[sid])

    dated_ids = []
    undated_ids = []
    for sid in sample_ids:
        s = sample_map[sid]
        if s.collection_date is not None:
            dated_ids.append(sid)
        else:
            undated_ids.append(sid)

    edges: List[Tuple[int, int]] = []
    for i, sid_a in enumerate(dated_ids):
        for j, sid_b in enumerate(dated_ids):
            if i == j:
                continue
            set_a = variant_sets[sid_a]
            set_b = variant_sets[sid_b]
            if len(set_a) == 0:
                continue
            if set_a < set_b:
                date_a = sample_map[sid_a].collection_date
                date_b = sample_map[sid_b].collection_date
                if date_a < date_b:
                    edges.append((sid_a, sid_b))

    reduced_edges = transitive_reduction(edges, set(dated_ids))

    ref_name_map = _ref_name_cache(db)

    nodes_out: List[schemas.TransmissionNodeOut] = []
    for sid in dated_ids:
        s = sample_map[sid]
        nodes_out.append(schemas.TransmissionNodeOut(
            sample_id=sid,
            sample_name=s.name,
            collection_date=s.collection_date,
            variant_count=len(variant_sets[sid]),
            date_unknown=False,
        ))
    for sid in undated_ids:
        s = sample_map[sid]
        nodes_out.append(schemas.TransmissionNodeOut(
            sample_id=sid,
            sample_name=s.name,
            collection_date=None,
            variant_count=len(variant_sets[sid]),
            date_unknown=True,
        ))

    edges_out: List[schemas.TransmissionEdgeOut] = []
    for src_id, tgt_id in reduced_edges:
        new_variant_keys = variant_sets[tgt_id] - variant_sets[src_id]
        new_variants_out = []
        for key in sorted(new_variant_keys, key=lambda k: (k[1], k[2])):
            detail = variant_details[tgt_id].get(key)
            if detail:
                ref_name = ref_name_map.get(detail.reference_id, f"ref_{detail.reference_id}")
                new_variants_out.append(schemas.SpectrumVariantOut(
                    reference_name=ref_name,
                    ref_pos=detail.ref_pos,
                    variant_type=detail.variant_type,
                    ref_base=detail.ref_base,
                    alt_base=detail.alt_base,
                    feature_type=detail.feature_type,
                    gene_name=detail.gene_name,
                    impact=detail.impact,
                    consequence=detail.consequence,
                    aa_ref=detail.aa_ref,
                    aa_alt=detail.aa_alt,
                    source_alignment_ids=detail.source_alignment_ids or [],
                ))
        edges_out.append(schemas.TransmissionEdgeOut(
            source_id=src_id,
            source_name=sample_map[src_id].name,
            target_id=tgt_id,
            target_name=sample_map[tgt_id].name,
            new_variants=new_variants_out,
        ))

    return schemas.TransmissionChainOut(
        sample_ids=sample_ids,
        total_nodes=len(nodes_out),
        total_edges=len(edges_out),
        nodes=nodes_out,
        edges=edges_out,
        excluded_no_date=undated_ids,
    )


def shared_variant_clustering(
    db: Session,
    request: schemas.ClusteringRequest,
) -> schemas.ClusteringOut:
    ref = db.query(models.ReferenceSequence).filter(
        models.ReferenceSequence.name == request.reference_name
    ).first()
    if not ref:
        raise ValueError(f"Reference sequence '{request.reference_name}' not found")

    sample_ids_rows = db.query(models.SampleVariantSpectrum.sample_id).filter(
        models.SampleVariantSpectrum.reference_id == ref.id
    ).distinct().all()
    sample_ids = [r[0] for r in sample_ids_rows]

    if len(sample_ids) < 2:
        return schemas.ClusteringOut(
            reference_name=request.reference_name,
            total_samples=len(sample_ids),
            similarity_threshold=request.similarity_threshold,
            cluster_count=len(sample_ids),
            clusters=[
                schemas.ClusterOut(
                    cluster_index=0,
                    member_count=1,
                    members=[
                        schemas.ClusterMemberOut(
                            sample_id=sid,
                            sample_name=_get_sample_name(db, sid),
                            collection_date=_get_sample_date(db, sid),
                            variant_count=_get_variant_count(db, sid, ref.id),
                        )
                    ],
                    shared_variants=[],
                    avg_jaccard=1.0,
                )
                for sid in sample_ids
            ],
        )

    variant_sets, variant_details = _get_variant_sets_for_samples(
        db, sample_ids, ref.id
    )

    jaccard_matrix = compute_jaccard_matrix(variant_sets, sample_ids)
    dist_matrix = compute_distance_matrix_from_jaccard(jaccard_matrix)

    labels = [str(sid) for sid in sample_ids]

    clusters_indices = upgma_cluster(dist_matrix, labels, request.similarity_threshold)

    samples = db.query(models.Sample).filter(
        models.Sample.id.in_(sample_ids)
    ).all()
    sample_map: Dict[int, models.Sample] = {s.id: s for s in samples}
    ref_name_map = _ref_name_cache(db)

    clusters_out: List[schemas.ClusterOut] = []
    for idx, cluster_idx_list in enumerate(clusters_indices):
        member_ids = [sample_ids[i] for i in cluster_idx_list]

        members = []
        for sid in member_ids:
            s = sample_map.get(sid)
            members.append(schemas.ClusterMemberOut(
                sample_id=sid,
                sample_name=s.name if s else f"sample_{sid}",
                collection_date=s.collection_date if s else None,
                variant_count=len(variant_sets.get(sid, set())),
            ))

        if len(member_ids) > 0:
            shared_keys = variant_sets.get(member_ids[0], set()).copy()
            for sid in member_ids[1:]:
                shared_keys &= variant_sets.get(sid, set())
        else:
            shared_keys = set()

        shared_variants_out = []
        for key in sorted(shared_keys, key=lambda k: (k[1], k[2])):
            any_detail = None
            for sid in member_ids:
                detail = variant_details.get(sid, {}).get(key)
                if detail:
                    any_detail = detail
                    break
            if any_detail:
                ref_name = ref_name_map.get(any_detail.reference_id, f"ref_{any_detail.reference_id}")
                shared_variants_out.append(schemas.SpectrumVariantOut(
                    reference_name=ref_name,
                    ref_pos=any_detail.ref_pos,
                    variant_type=any_detail.variant_type,
                    ref_base=any_detail.ref_base,
                    alt_base=any_detail.alt_base,
                    feature_type=any_detail.feature_type,
                    gene_name=any_detail.gene_name,
                    impact=any_detail.impact,
                    consequence=any_detail.consequence,
                    aa_ref=any_detail.aa_ref,
                    aa_alt=any_detail.aa_alt,
                    source_alignment_ids=any_detail.source_alignment_ids or [],
                ))

        pair_sims = []
        for i in range(len(member_ids)):
            for j in range(i + 1, len(member_ids)):
                mi = member_ids[i]
                mj = member_ids[j]
                ii = sample_ids.index(mi)
                ij = sample_ids.index(mj)
                pair_sims.append(jaccard_matrix[ii][ij])

        avg_j = sum(pair_sims) / len(pair_sims) if pair_sims else 1.0

        clusters_out.append(schemas.ClusterOut(
            cluster_index=idx,
            member_count=len(members),
            members=members,
            shared_variants=shared_variants_out,
            avg_jaccard=round(avg_j, 6),
        ))

    return schemas.ClusteringOut(
        reference_name=request.reference_name,
        total_samples=len(sample_ids),
        similarity_threshold=request.similarity_threshold,
        cluster_count=len(clusters_out),
        clusters=clusters_out,
    )


def _get_sample_name(db: Session, sample_id: int) -> str:
    s = db.query(models.Sample).filter(models.Sample.id == sample_id).first()
    return s.name if s else f"sample_{sample_id}"


def _get_sample_date(db: Session, sample_id: int) -> Optional[datetime]:
    s = db.query(models.Sample).filter(models.Sample.id == sample_id).first()
    return s.collection_date if s else None


def _get_variant_count(db: Session, sample_id: int, reference_id: int) -> int:
    return db.query(func.count(models.SampleVariantSpectrum.id)).filter(
        models.SampleVariantSpectrum.sample_id == sample_id,
        models.SampleVariantSpectrum.reference_id == reference_id,
    ).scalar() or 0
