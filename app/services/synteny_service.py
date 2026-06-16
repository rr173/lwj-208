import hashlib
from typing import Optional, List
from sqlalchemy.orm import Session
from app import models, schemas
from app.synteny.synteny_analysis import (
    analyze_synteny,
    sequence_pair_hash,
    SyntenyBlock,
    SyntenyAnchor,
    RearrangementEvent,
    PairwiseSyntenyResult,
    infer_conserved_core_regions,
)
from app.alignment.smith_waterman import AlignmentTimeoutError


def _compute_data_hash(ref_a: models.ReferenceSequence, ref_b: models.ReferenceSequence) -> str:
    combined = f"{ref_a.name}:{ref_a.sequence}:{ref_b.name}:{ref_b.sequence}"
    return hashlib.md5(combined.encode()).hexdigest()


def _build_cache_key(
    ref_a_name: str,
    ref_b_name: str,
    anchor_length: int,
    score_threshold_ratio: float,
    max_gap_ratio: float,
) -> str:
    sorted_names = sorted([ref_a_name, ref_b_name])
    key_str = f"{sorted_names[0]}__{sorted_names[1]}__{anchor_length}__{score_threshold_ratio}__{max_gap_ratio}"
    return hashlib.md5(key_str.encode()).hexdigest()


def get_cached_synteny(
    db: Session,
    cache_key: str,
) -> Optional[models.SyntenyResult]:
    return db.query(models.SyntenyResult).filter(
        models.SyntenyResult.cache_key == cache_key
    ).first()


def _block_to_dict(block: SyntenyBlock) -> dict:
    anchors = []
    for anchor in block.anchors:
        anchors.append({
            "a_start": anchor.a_start,
            "a_end": anchor.a_end,
            "b_start": anchor.b_start,
            "b_end": anchor.b_end,
            "score": anchor.score,
            "direction": anchor.direction,
        })
    return {
        "a_start": block.a_start,
        "a_end": block.a_end,
        "b_start": block.b_start,
        "b_end": block.b_end,
        "anchor_count": block.anchor_count,
        "direction": block.direction,
        "anchors": anchors,
    }


def _rearrangement_to_dict(event: RearrangementEvent) -> dict:
    return {
        "event_type": event.event_type,
        "a_start": event.a_start,
        "a_end": event.a_end,
        "b_start": event.b_start,
        "b_end": event.b_end,
        "anchor_count": event.anchor_count,
    }


def perform_synteny_analysis(
    db: Session,
    ref_a_name: str,
    ref_b_name: str,
    anchor_length: int = 50,
    score_threshold_ratio: float = 1.5,
    max_gap_ratio: float = 3.0,
) -> schemas.SyntenyAnalysisOut:
    all_refs = db.query(models.ReferenceSequence).all()
    if len(all_refs) < 2:
        raise ValueError(
            "At least two reference sequences are required for synteny analysis. "
            f"Currently only {len(all_refs)} reference(s) exist."
        )

    ref_a = db.query(models.ReferenceSequence).filter(
        models.ReferenceSequence.name == ref_a_name
    ).first()
    if not ref_a:
        raise ValueError(f"Reference sequence '{ref_a_name}' not found")

    ref_b = db.query(models.ReferenceSequence).filter(
        models.ReferenceSequence.name == ref_b_name
    ).first()
    if not ref_b:
        raise ValueError(f"Reference sequence '{ref_b_name}' not found")

    if ref_a.id == ref_b.id:
        raise ValueError("ref_a_name and ref_b_name must be different reference sequences")

    cache_key = _build_cache_key(ref_a_name, ref_b_name, anchor_length, score_threshold_ratio, max_gap_ratio)
    cached = get_cached_synteny(db, cache_key)

    data_hash = _compute_data_hash(ref_a, ref_b)

    if cached and not cached.is_stale and cached.data_hash == data_hash:
        return schemas.SyntenyAnalysisOut(
            ref_a_name=ref_a.name,
            ref_b_name=ref_b.name,
            seq_a_length=cached.seq_a_length,
            seq_b_length=cached.seq_b_length,
            anchor_length=cached.anchor_length,
            score_threshold_ratio=cached.score_threshold_ratio,
            max_gap_ratio=cached.max_gap_ratio,
            total_anchors_aligned=cached.total_anchors_aligned,
            synteny_blocks=[schemas.SyntenyBlockOut(**b) for b in cached.synteny_blocks],
            rearrangements=[schemas.RearrangementEventOut(**r) for r in cached.rearrangements],
            b_non_monotonic_transitions=cached.b_non_monotonic_transitions,
            is_stale=cached.is_stale,
            cached=True,
            created_at=cached.created_at,
        )

    result = analyze_synteny(
        ref_a.sequence,
        ref_b.sequence,
        anchor_length=anchor_length,
        score_threshold_ratio=score_threshold_ratio,
        max_gap_ratio=max_gap_ratio,
    )

    blocks_dict = [_block_to_dict(b) for b in result["synteny_blocks"]]
    rearrangements_dict = [_rearrangement_to_dict(r) for r in result["rearrangements"]]

    if cached:
        cached.seq_a_length = result["seq_a_length"]
        cached.seq_b_length = result["seq_b_length"]
        cached.anchor_length = anchor_length
        cached.score_threshold_ratio = score_threshold_ratio
        cached.max_gap_ratio = max_gap_ratio
        cached.total_anchors_aligned = result["total_anchors_aligned"]
        cached.synteny_blocks = blocks_dict
        cached.rearrangements = rearrangements_dict
        cached.b_non_monotonic_transitions = result["b_non_monotonic_transitions"]
        cached.data_hash = data_hash
        cached.is_stale = False
        db.commit()
        db.refresh(cached)
    else:
        synteny_result = models.SyntenyResult(
            cache_key=cache_key,
            ref_a_name=ref_a.name,
            ref_b_name=ref_b.name,
            ref_a_id=ref_a.id,
            ref_b_id=ref_b.id,
            seq_a_length=result["seq_a_length"],
            seq_b_length=result["seq_b_length"],
            anchor_length=anchor_length,
            score_threshold_ratio=score_threshold_ratio,
            max_gap_ratio=max_gap_ratio,
            total_anchors_aligned=result["total_anchors_aligned"],
            synteny_blocks=blocks_dict,
            rearrangements=rearrangements_dict,
            b_non_monotonic_transitions=result["b_non_monotonic_transitions"],
            data_hash=data_hash,
            is_stale=False,
        )
        db.add(synteny_result)
        db.commit()
        db.refresh(synteny_result)
        cached = synteny_result

    return schemas.SyntenyAnalysisOut(
        ref_a_name=ref_a.name,
        ref_b_name=ref_b.name,
        seq_a_length=cached.seq_a_length,
        seq_b_length=cached.seq_b_length,
        anchor_length=cached.anchor_length,
        score_threshold_ratio=cached.score_threshold_ratio,
        max_gap_ratio=cached.max_gap_ratio,
        total_anchors_aligned=cached.total_anchors_aligned,
        synteny_blocks=[schemas.SyntenyBlockOut(**b) for b in cached.synteny_blocks],
        rearrangements=[schemas.RearrangementEventOut(**r) for r in cached.rearrangements],
        b_non_monotonic_transitions=cached.b_non_monotonic_transitions,
        is_stale=cached.is_stale,
        cached=False,
        created_at=cached.created_at,
    )


def list_synteny_results(
    db: Session,
    skip: int = 0,
    limit: int = 100,
) -> List[models.SyntenyResult]:
    return db.query(models.SyntenyResult).order_by(
        models.SyntenyResult.created_at.desc()
    ).offset(skip).limit(limit).all()


def get_synteny_by_id(
    db: Session,
    result_id: int,
) -> Optional[models.SyntenyResult]:
    return db.query(models.SyntenyResult).filter(
        models.SyntenyResult.id == result_id
    ).first()


def mark_synteny_stale_for_reference(
    db: Session,
    reference_id: int,
):
    results = db.query(models.SyntenyResult).filter(
        (models.SyntenyResult.ref_a_id == reference_id)
        | (models.SyntenyResult.ref_b_id == reference_id)
    ).all()
    for r in results:
        r.is_stale = True
    db.commit()


def _dict_to_block(block_dict: dict) -> Optional[SyntenyBlock]:
    if not block_dict:
        return None

    anchors = []
    for anchor_dict in block_dict.get("anchors", []):
        anchors.append(SyntenyAnchor(
            a_start=anchor_dict["a_start"],
            a_end=anchor_dict["a_end"],
            b_start=anchor_dict["b_start"],
            b_end=anchor_dict["b_end"],
            score=anchor_dict["score"],
            direction=anchor_dict.get("direction", "+"),
        ))

    return SyntenyBlock(
        a_start=block_dict["a_start"],
        a_end=block_dict["a_end"],
        b_start=block_dict["b_start"],
        b_end=block_dict["b_end"],
        anchor_count=block_dict["anchor_count"],
        direction=block_dict["direction"],
        anchors=anchors,
    )


def _get_pairwise_with_anchors(
    db: Session,
    ref_a_name: str,
    ref_b_name: str,
    anchor_length: int,
    score_threshold_ratio: float,
    max_gap_ratio: float,
) -> PairwiseSyntenyResult:
    ref_a = db.query(models.ReferenceSequence).filter(
        models.ReferenceSequence.name == ref_a_name
    ).first()
    if not ref_a:
        raise ValueError(f"Reference sequence '{ref_a_name}' not found")

    ref_b = db.query(models.ReferenceSequence).filter(
        models.ReferenceSequence.name == ref_b_name
    ).first()
    if not ref_b:
        raise ValueError(f"Reference sequence '{ref_b_name}' not found")

    if ref_a.id == ref_b.id:
        raise ValueError("ref_a_name and ref_b_name must be different reference sequences")

    cache_key = _build_cache_key(ref_a_name, ref_b_name, anchor_length, score_threshold_ratio, max_gap_ratio)
    cached = get_cached_synteny(db, cache_key)
    data_hash = _compute_data_hash(ref_a, ref_b)

    if cached and not cached.is_stale and cached.data_hash == data_hash:
        has_anchors = all(
            "anchors" in b for b in cached.synteny_blocks
        ) if cached.synteny_blocks else False

        if has_anchors:
            blocks = []
            for block_dict in cached.synteny_blocks:
                block = _dict_to_block(block_dict)
                if block:
                    blocks.append(block)
            return PairwiseSyntenyResult(
                ref_a_name=ref_a.name,
                ref_b_name=ref_b.name,
                blocks=blocks,
                total_anchors=cached.total_anchors_aligned,
                cached=True,
            )

    result = analyze_synteny(
        ref_a.sequence,
        ref_b.sequence,
        anchor_length=anchor_length,
        score_threshold_ratio=score_threshold_ratio,
        max_gap_ratio=max_gap_ratio,
    )

    blocks_dict = [_block_to_dict(b) for b in result["synteny_blocks"]]
    rearrangements_dict = [_rearrangement_to_dict(r) for r in result["rearrangements"]]

    if cached:
        cached.seq_a_length = result["seq_a_length"]
        cached.seq_b_length = result["seq_b_length"]
        cached.anchor_length = anchor_length
        cached.score_threshold_ratio = score_threshold_ratio
        cached.max_gap_ratio = max_gap_ratio
        cached.total_anchors_aligned = result["total_anchors_aligned"]
        cached.synteny_blocks = blocks_dict
        cached.rearrangements = rearrangements_dict
        cached.b_non_monotonic_transitions = result["b_non_monotonic_transitions"]
        cached.data_hash = data_hash
        cached.is_stale = False
        db.commit()
    else:
        synteny_result = models.SyntenyResult(
            cache_key=cache_key,
            ref_a_name=ref_a.name,
            ref_b_name=ref_b.name,
            ref_a_id=ref_a.id,
            ref_b_id=ref_b.id,
            seq_a_length=result["seq_a_length"],
            seq_b_length=result["seq_b_length"],
            anchor_length=anchor_length,
            score_threshold_ratio=score_threshold_ratio,
            max_gap_ratio=max_gap_ratio,
            total_anchors_aligned=result["total_anchors_aligned"],
            synteny_blocks=blocks_dict,
            rearrangements=rearrangements_dict,
            b_non_monotonic_transitions=result["b_non_monotonic_transitions"],
            data_hash=data_hash,
            is_stale=False,
        )
        db.add(synteny_result)
        db.commit()

    return PairwiseSyntenyResult(
        ref_a_name=ref_a.name,
        ref_b_name=ref_b.name,
        blocks=result["synteny_blocks"],
        total_anchors=result["total_anchors_aligned"],
        cached=False,
    )


def perform_multi_synteny_analysis(
    db: Session,
    reference_names: List[str],
    anchor_length: int = 50,
    score_threshold_ratio: float = 1.5,
    max_gap_ratio: float = 3.0,
) -> schemas.MultiSyntenyOut:
    if len(reference_names) < 3:
        raise ValueError("At least 3 reference sequences are required for multi-sequence synteny analysis")
    if len(reference_names) > 10:
        raise ValueError("Maximum 10 reference sequences are allowed for multi-sequence synteny analysis")

    unique_names = list(dict.fromkeys(reference_names))
    if len(unique_names) != len(reference_names):
        raise ValueError("Duplicate reference names are not allowed")

    ref_map = {}
    ref_lengths = {}
    for name in reference_names:
        ref = db.query(models.ReferenceSequence).filter(
            models.ReferenceSequence.name == name
        ).first()
        if not ref:
            raise ValueError(f"Reference sequence '{name}' not found")
        ref_map[name] = ref
        ref_lengths[name] = ref.length

    primary_name = reference_names[0]

    pairwise_results = []
    pairwise_summaries = []
    cached_count = 0

    from itertools import combinations
    for name_a, name_b in combinations(reference_names, 2):
        if name_a == primary_name or name_b == primary_name:
            result = _get_pairwise_with_anchors(
                db,
                name_a,
                name_b,
                anchor_length,
                score_threshold_ratio,
                max_gap_ratio,
            )

            if name_a != primary_name:
                swapped_blocks = []
                for block in result.blocks:
                    swapped = SyntenyBlock(
                        a_start=block.b_start,
                        a_end=block.b_end,
                        b_start=block.a_start,
                        b_end=block.a_end,
                        anchor_count=block.anchor_count,
                        direction=block.direction,
                        anchors=[
                            SyntenyAnchor(
                                a_start=a.b_start,
                                a_end=a.b_end,
                                b_start=a.a_start,
                                b_end=a.a_end,
                                score=a.score,
                                direction=a.direction,
                            )
                            for a in block.anchors
                        ],
                    )
                    swapped_blocks.append(swapped)
                result = PairwiseSyntenyResult(
                    ref_a_name=primary_name,
                    ref_b_name=name_a,
                    blocks=swapped_blocks,
                    total_anchors=result.total_anchors,
                    cached=result.cached,
                )

            pairwise_results.append(result)
            pairwise_summaries.append(schemas.PairwiseSyntenySummary(
                ref_a_name=result.ref_a_name,
                ref_b_name=result.ref_b_name,
                block_count=len(result.blocks),
                total_anchors=result.total_anchors,
                cached=result.cached,
            ))
            if result.cached:
                cached_count += 1

    conserved_result = infer_conserved_core_regions(
        pairwise_results,
        reference_names,
        ref_lengths,
        anchor_length,
    )

    core_regions = []
    for region_dict in conserved_result["conserved_core_regions"]:
        intervals = []
        for interval_dict in region_dict["intervals"]:
            intervals.append(schemas.SequenceInterval(**interval_dict))
        core_regions.append(schemas.ConservedCoreRegion(
            region_index=region_dict["region_index"],
            intervals=intervals,
            avg_anchor_score=region_dict["avg_anchor_score"],
            sequence_count=region_dict["sequence_count"],
            total_length=region_dict["total_length"],
        ))

    return schemas.MultiSyntenyOut(
        reference_names=reference_names,
        first_reference_name=primary_name,
        first_reference_length=ref_lengths[primary_name],
        conserved_core_regions=core_regions,
        total_conserved_length=conserved_result["total_conserved_length"],
        conservation_ratio=conserved_result["conservation_ratio"],
        pairwise_summaries=pairwise_summaries,
        anchor_length=anchor_length,
        score_threshold_ratio=score_threshold_ratio,
        max_gap_ratio=max_gap_ratio,
        total_pairwise_comparisons=len(pairwise_summaries),
        cached_pairwise_count=cached_count,
    )
