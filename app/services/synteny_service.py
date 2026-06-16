import hashlib
from typing import Optional, List
from sqlalchemy.orm import Session
from app import models, schemas
from app.synteny.synteny_analysis import (
    analyze_synteny,
    sequence_pair_hash,
    SyntenyBlock,
    RearrangementEvent,
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
) -> str:
    sorted_names = sorted([ref_a_name, ref_b_name])
    key_str = f"{sorted_names[0]}__{sorted_names[1]}__{anchor_length}__{score_threshold_ratio}"
    return hashlib.md5(key_str.encode()).hexdigest()


def get_cached_synteny(
    db: Session,
    cache_key: str,
) -> Optional[models.SyntenyResult]:
    return db.query(models.SyntenyResult).filter(
        models.SyntenyResult.cache_key == cache_key
    ).first()


def _block_to_dict(block: SyntenyBlock) -> dict:
    return {
        "a_start": block.a_start,
        "a_end": block.a_end,
        "b_start": block.b_start,
        "b_end": block.b_end,
        "anchor_count": block.anchor_count,
        "direction": block.direction,
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

    cache_key = _build_cache_key(ref_a_name, ref_b_name, anchor_length, score_threshold_ratio)
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
            total_anchors_aligned=cached.total_anchors_aligned,
            synteny_blocks=[schemas.SyntenyBlockOut(**b) for b in cached.synteny_blocks],
            rearrangements=[schemas.RearrangementEventOut(**r) for r in cached.rearrangements],
            is_stale=cached.is_stale,
            cached=True,
            created_at=cached.created_at,
        )

    result = analyze_synteny(
        ref_a.sequence,
        ref_b.sequence,
        anchor_length=anchor_length,
        score_threshold_ratio=score_threshold_ratio,
    )

    blocks_dict = [_block_to_dict(b) for b in result["synteny_blocks"]]
    rearrangements_dict = [_rearrangement_to_dict(r) for r in result["rearrangements"]]

    if cached:
        cached.seq_a_length = result["seq_a_length"]
        cached.seq_b_length = result["seq_b_length"]
        cached.anchor_length = anchor_length
        cached.score_threshold_ratio = score_threshold_ratio
        cached.total_anchors_aligned = result["total_anchors_aligned"]
        cached.synteny_blocks = blocks_dict
        cached.rearrangements = rearrangements_dict
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
            total_anchors_aligned=result["total_anchors_aligned"],
            synteny_blocks=blocks_dict,
            rearrangements=rearrangements_dict,
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
        total_anchors_aligned=cached.total_anchors_aligned,
        synteny_blocks=[schemas.SyntenyBlockOut(**b) for b in cached.synteny_blocks],
        rearrangements=[schemas.RearrangementEventOut(**r) for r in cached.rearrangements],
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
