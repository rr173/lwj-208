from fastapi import APIRouter, Depends, HTTPException
from sqlalchemy.orm import Session
from app import schemas
from app.database import get_db
from app.services import transmission_service

router = APIRouter(prefix="/api/transmission", tags=["transmission"])


@router.post("/trace", response_model=schemas.VariantTracingOut)
def trace_variant(
    request: schemas.VariantTracingRequest,
    db: Session = Depends(get_db),
):
    try:
        return transmission_service.variant_tracing(db, request)
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))


@router.post("/chain", response_model=schemas.TransmissionChainOut)
def infer_transmission_chain(
    request: schemas.TransmissionChainRequest,
    db: Session = Depends(get_db),
):
    if len(set(request.sample_ids)) < 3:
        raise HTTPException(
            status_code=400,
            detail="At least 3 distinct sample IDs are required for transmission chain inference",
        )
    try:
        return transmission_service.transmission_chain(db, request)
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))


@router.post("/cluster", response_model=schemas.ClusteringOut)
def cluster_shared_variants(
    request: schemas.ClusteringRequest,
    db: Session = Depends(get_db),
):
    try:
        return transmission_service.shared_variant_clustering(db, request)
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))
