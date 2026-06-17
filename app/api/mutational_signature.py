from fastapi import APIRouter, Depends, HTTPException, Query
from sqlalchemy.orm import Session
from typing import List

from app import schemas
from app.database import get_db
from app.services import mutational_signature_service

router = APIRouter(prefix="/api/mutational-signatures", tags=["mutational-signatures"])


@router.post("/spectrum", response_model=schemas.MutationalSpectrumOut)
def build_spectrum(
    request: schemas.MutationalSpectrumRequest,
    include_failed_qc: bool = Query(False, description="Include variants that failed QC"),
    db: Session = Depends(get_db),
):
    try:
        return mutational_signature_service.build_mutational_spectrum(
            db,
            sample_ids=request.sample_ids,
            reference_name=request.reference_name,
            include_failed_qc=include_failed_qc,
        )
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))


@router.post("/extract", response_model=schemas.NMFSignatureOut)
def extract_signatures(
    request: schemas.NMFSignatureRequest,
    include_failed_qc: bool = Query(False, description="Include variants that failed QC"),
    db: Session = Depends(get_db),
):
    try:
        return mutational_signature_service.extract_signatures(
            db,
            sample_ids=request.sample_ids,
            reference_name=request.reference_name,
            k_value=request.k_value,
            include_failed_qc=include_failed_qc,
        )
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))


@router.post("/optimal-k", response_model=schemas.OptimalKOut)
def find_optimal_k(
    request: schemas.OptimalKRequest,
    include_failed_qc: bool = Query(False, description="Include variants that failed QC"),
    db: Session = Depends(get_db),
):
    try:
        return mutational_signature_service.find_optimal_k_value(
            db,
            sample_ids=request.sample_ids,
            reference_name=request.reference_name,
            include_failed_qc=include_failed_qc,
        )
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))


@router.post("/match", response_model=schemas.SignatureComparisonOut)
def match_signature(
    request: schemas.SignatureMatchRequest,
    db: Session = Depends(get_db),
):
    try:
        return mutational_signature_service.match_signature_to_reference(
            db,
            query_signature=request.signature,
            top_n=3,
            similarity_threshold=0.8,
        )
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))


@router.get("/reference-signatures", response_model=schemas.ReferenceSignatureListOut)
def list_reference_signatures(
    db: Session = Depends(get_db),
):
    return mutational_signature_service.list_reference_signatures(db)


@router.post("/temporal-analysis", response_model=schemas.SignatureTemporalAnalysisOut)
def temporal_analysis(
    request: schemas.TemporalAnalysisRequest,
    include_failed_qc: bool = Query(False, description="Include variants that failed QC"),
    db: Session = Depends(get_db),
):
    try:
        return mutational_signature_service.analyze_signature_temporal_trends(
            db,
            sample_ids=request.sample_ids,
            reference_name=request.reference_name,
            k_value=request.k_value,
            p_value_threshold=request.p_value_threshold,
            include_failed_qc=include_failed_qc,
        )
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))


@router.post("/exposure-change", response_model=schemas.ExposureChangeOut)
def exposure_change(
    request: schemas.ExposureChangeRequest,
    include_failed_qc: bool = Query(False, description="Include variants that failed QC"),
    db: Session = Depends(get_db),
):
    try:
        return mutational_signature_service.calculate_exposure_change_between_dates(
            db,
            sample_ids=request.sample_ids,
            reference_name=request.reference_name,
            start_date=request.start_date,
            end_date=request.end_date,
            k_value=request.k_value,
            aggregation=request.aggregation,
            include_failed_qc=include_failed_qc,
        )
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))
