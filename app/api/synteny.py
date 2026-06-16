from fastapi import APIRouter, Depends, HTTPException
from sqlalchemy.orm import Session
from typing import List
from app import schemas
from app.database import get_db
from app.services import synteny_service
from app.alignment.smith_waterman import AlignmentTimeoutError

router = APIRouter(prefix="/api/synteny", tags=["synteny"])


def _result_to_out(result) -> schemas.SyntenyAnalysisOut:
    return schemas.SyntenyAnalysisOut(
        ref_a_name=result.ref_a_name,
        ref_b_name=result.ref_b_name,
        seq_a_length=result.seq_a_length,
        seq_b_length=result.seq_b_length,
        anchor_length=result.anchor_length,
        score_threshold_ratio=result.score_threshold_ratio,
        max_gap_ratio=result.max_gap_ratio,
        total_anchors_aligned=result.total_anchors_aligned,
        synteny_blocks=[schemas.SyntenyBlockOut(**b) for b in result.synteny_blocks],
        rearrangements=[schemas.RearrangementEventOut(**r) for r in result.rearrangements],
        b_non_monotonic_transitions=result.b_non_monotonic_transitions,
        is_stale=result.is_stale,
        cached=True,
        created_at=result.created_at,
    )


@router.post("/analyze", response_model=schemas.SyntenyAnalysisOut)
def analyze_synteny(request: schemas.SyntenyAnalysisRequest, db: Session = Depends(get_db)):
    try:
        result = synteny_service.perform_synteny_analysis(
            db,
            ref_a_name=request.ref_a_name,
            ref_b_name=request.ref_b_name,
            anchor_length=request.anchor_length,
            score_threshold_ratio=request.score_threshold_ratio,
            max_gap_ratio=request.max_gap_ratio,
        )
        return result
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))
    except AlignmentTimeoutError as e:
        raise HTTPException(status_code=408, detail=str(e))


@router.get("/results", response_model=List[schemas.SyntenyAnalysisOut])
def list_synteny_results(skip: int = 0, limit: int = 100, db: Session = Depends(get_db)):
    results = synteny_service.list_synteny_results(db, skip=skip, limit=limit)
    return [_result_to_out(r) for r in results]


@router.get("/results/{result_id}", response_model=schemas.SyntenyAnalysisOut)
def get_synteny_result(result_id: int, db: Session = Depends(get_db)):
    result = synteny_service.get_synteny_by_id(db, result_id)
    if not result:
        raise HTTPException(status_code=404, detail="Synteny result not found")
    return _result_to_out(result)
