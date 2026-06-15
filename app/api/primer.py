from fastapi import APIRouter, Depends, HTTPException
from sqlalchemy.orm import Session
from app import schemas
from app.database import get_db
from app.services.genome_service import get_reference_by_name
from app.primer.primer_design import design_primers, get_cached_primer_result, save_primer_result

router = APIRouter(prefix="/api/primer", tags=["primer"])


@router.post("/design", response_model=schemas.PrimerDesignResult)
def design_primer_pair(request: schemas.PrimerDesignRequest, db: Session = Depends(get_db)):
    if request.target_start >= request.target_end:
        raise HTTPException(
            status_code=400,
            detail="target_start must be less than target_end",
        )

    reference = get_reference_by_name(db, request.reference_name)
    if not reference:
        raise HTTPException(
            status_code=404,
            detail=f"Reference sequence '{request.reference_name}' not found",
        )

    if request.target_end >= reference.length:
        raise HTTPException(
            status_code=400,
            detail=f"target_end ({request.target_end}) exceeds reference length ({reference.length})",
        )

    cached = get_cached_primer_result(db, request.reference_name, request.target_start, request.target_end)
    if cached is not None:
        cached["cached"] = True
        return cached

    result = design_primers(
        reference_sequence=reference.sequence,
        ref_name=request.reference_name,
        target_start=request.target_start,
        target_end=request.target_end,
    )

    if "error" in result:
        raise HTTPException(status_code=422, detail=result["error"])

    result["cached"] = False

    save_primer_result(db, request.reference_name, request.target_start, request.target_end, result)

    return result
