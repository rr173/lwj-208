from typing import Optional
from fastapi import APIRouter, Depends, HTTPException, Query
from sqlalchemy.orm import Session
from app import schemas
from app.database import get_db
from app.services import ld_service

router = APIRouter(prefix="/api/ld", tags=["linkage_disequilibrium"])


@router.get("/{reference_name}/matrix", response_model=schemas.LDMatrixOut)
def get_ld_matrix(
    reference_name: str,
    start: Optional[int] = None,
    end: Optional[int] = None,
    min_maf: float = 0.0,
    include_failed_qc: bool = Query(False, description="Include variants that failed QC"),
    db: Session = Depends(get_db)
):
    try:
        result = ld_service.compute_ld_matrix_endpoint(
            db, reference_name, start=start, end=end, min_maf=min_maf,
            include_failed_qc=include_failed_qc
        )
        return result
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))


@router.get("/{reference_name}/haplotype-blocks", response_model=schemas.HaplotypeBlockOut)
def get_haplotype_blocks(
    reference_name: str,
    r2_threshold: float = 0.8,
    start: Optional[int] = None,
    end: Optional[int] = None,
    min_maf: float = 0.0,
    include_failed_qc: bool = Query(False, description="Include variants that failed QC"),
    db: Session = Depends(get_db)
):
    try:
        result = ld_service.compute_haplotype_blocks_endpoint(
            db, reference_name, r2_threshold, start=start, end=end, min_maf=min_maf,
            include_failed_qc=include_failed_qc
        )
        return result
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))
