from fastapi import APIRouter, Depends, HTTPException
from sqlalchemy.orm import Session
from typing import List
from app import schemas
from app.database import get_db
from app.services import genome_service

router = APIRouter(prefix="/api/alignment", tags=["alignment"])


def alignment_to_out(alignment) -> schemas.AlignmentResultOut:
    return schemas.AlignmentResultOut(
        id=alignment.id,
        reference_name=alignment.reference.name if alignment.reference else "",
        score=alignment.score,
        ref_start=alignment.ref_start,
        ref_end=alignment.ref_end,
        query_start=alignment.query_start,
        query_end=alignment.query_end,
        cigar=alignment.cigar,
        alignment=schemas.AlignmentDetail(
            query=alignment.alignment_query,
            match=alignment.alignment_match,
            reference=alignment.alignment_ref,
        ),
        variants=alignment.variants,
    )


@router.post("/align", response_model=List[schemas.AlignmentResultOut])
def align_sequence(request: schemas.AlignRequest, db: Session = Depends(get_db)):
    if len(request.query_sequence) > 1000:
        raise HTTPException(status_code=400, detail="Query sequence exceeds 1000 base limit")

    results = genome_service.align_query_all_references(db, request.query_sequence)
    return [alignment_to_out(r) for r in results]


@router.get("/{alignment_id}", response_model=schemas.AlignmentResultOut)
def get_alignment(alignment_id: int, db: Session = Depends(get_db)):
    alignment = genome_service.get_alignment_by_id(db, alignment_id)
    if not alignment:
        raise HTTPException(status_code=404, detail="Alignment result not found")
    return alignment_to_out(alignment)
