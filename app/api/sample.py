from fastapi import APIRouter, Depends, HTTPException, Query
from sqlalchemy.orm import Session
from typing import List, Optional
from app import schemas
from app.database import get_db
from app.services import sample_service

router = APIRouter(prefix="/api/samples", tags=["samples"])


@router.post("", response_model=schemas.SampleOut, status_code=201)
def create_sample(request: schemas.SampleCreate, db: Session = Depends(get_db)):
    existing = sample_service.get_sample_by_name(db, request.name)
    if existing:
        raise HTTPException(status_code=400, detail=f"Sample with name '{request.name}' already exists")
    sample = sample_service.create_sample(db, request)
    return sample_service.sample_to_out(db, sample)


@router.get("", response_model=List[schemas.SampleOut])
def list_samples(
    skip: int = Query(0, ge=0),
    limit: int = Query(100, ge=1, le=500),
    db: Session = Depends(get_db),
):
    samples = sample_service.list_samples(db, skip=skip, limit=limit)
    return [sample_service.sample_to_out(db, s) for s in samples]


@router.get("/{sample_id}", response_model=schemas.SampleDetail)
def get_sample(sample_id: int, db: Session = Depends(get_db)):
    sample = sample_service.get_sample_by_id(db, sample_id)
    if not sample:
        raise HTTPException(status_code=404, detail="Sample not found")
    return sample_service.sample_to_out(db, sample, with_detail=True)


@router.put("/{sample_id}", response_model=schemas.SampleOut)
def update_sample(sample_id: int, request: schemas.SampleUpdate, db: Session = Depends(get_db)):
    sample = sample_service.get_sample_by_id(db, sample_id)
    if not sample:
        raise HTTPException(status_code=404, detail="Sample not found")
    if request.name and request.name != sample.name:
        existing = sample_service.get_sample_by_name(db, request.name)
        if existing:
            raise HTTPException(status_code=400, detail=f"Sample with name '{request.name}' already exists")
    updated = sample_service.update_sample(db, sample, request)
    return sample_service.sample_to_out(db, updated)


@router.delete("/{sample_id}", status_code=204)
def delete_sample(sample_id: int, db: Session = Depends(get_db)):
    sample = sample_service.get_sample_by_id(db, sample_id)
    if not sample:
        raise HTTPException(status_code=404, detail="Sample not found")
    sample_service.delete_sample(db, sample)
    return None


@router.post("/{sample_id}/alignments", status_code=200)
def link_alignments(
    sample_id: int,
    request: schemas.LinkAlignmentsRequest,
    db: Session = Depends(get_db),
):
    sample = sample_service.get_sample_by_id(db, sample_id)
    if not sample:
        raise HTTPException(status_code=404, detail="Sample not found")
    added, invalid = sample_service.link_alignments_to_sample(db, sample, request.alignment_ids)
    return {
        "sample_id": sample_id,
        "added_count": added,
        "invalid_alignment_ids": invalid,
        "message": f"Linked {added} alignment(s) to sample"
    }


@router.delete("/{sample_id}/alignments/{alignment_id}", status_code=204)
def unlink_alignment(
    sample_id: int,
    alignment_id: int,
    db: Session = Depends(get_db),
):
    sample = sample_service.get_sample_by_id(db, sample_id)
    if not sample:
        raise HTTPException(status_code=404, detail="Sample not found")
    removed = sample_service.unlink_alignment_from_sample(db, sample, alignment_id)
    if not removed:
        raise HTTPException(status_code=404, detail="Alignment link not found")
    return None


@router.get("/{sample_id}/spectrum", response_model=schemas.SampleSpectrumOut)
def get_sample_spectrum(
    sample_id: int,
    reference_name: Optional[str] = Query(None, description="Filter by reference sequence name"),
    impact: Optional[str] = Query(None, description="Filter by impact level (HIGH, MODERATE, LOW, MODIFIER)"),
    gene_name: Optional[str] = Query(None, description="Filter by gene name"),
    variant_type: Optional[str] = Query(None, description="Filter by variant type (SNP, INS, DEL)"),
    db: Session = Depends(get_db),
):
    sample = sample_service.get_sample_by_id(db, sample_id)
    if not sample:
        raise HTTPException(status_code=404, detail="Sample not found")
    return sample_service.get_sample_spectrum(
        db,
        sample,
        reference_name=reference_name,
        impact_filter=impact,
        gene_filter=gene_name,
        variant_type_filter=variant_type,
    )


@router.get("/frequency/{reference_name}", response_model=schemas.PopulationFrequencyOut)
def get_population_frequency(
    reference_name: str,
    min_frequency: float = Query(0.0, ge=0.0, le=1.0, description="Minimum frequency threshold"),
    max_frequency: float = Query(1.0, ge=0.0, le=1.0, description="Maximum frequency threshold"),
    db: Session = Depends(get_db),
):
    if min_frequency > max_frequency:
        raise HTTPException(status_code=400, detail="min_frequency cannot be greater than max_frequency")
    return sample_service.get_population_frequency(
        db,
        reference_name,
        min_frequency=min_frequency,
        max_frequency=max_frequency,
    )


@router.get("/compare/{sample_a_id}/{sample_b_id}", response_model=schemas.CompareResultOut)
def compare_samples(
    sample_a_id: int,
    sample_b_id: int,
    db: Session = Depends(get_db),
):
    if sample_a_id == sample_b_id:
        raise HTTPException(status_code=400, detail="Cannot compare a sample with itself")
    sample_a = sample_service.get_sample_by_id(db, sample_a_id)
    if not sample_a:
        raise HTTPException(status_code=404, detail=f"Sample {sample_a_id} not found")
    sample_b = sample_service.get_sample_by_id(db, sample_b_id)
    if not sample_b:
        raise HTTPException(status_code=404, detail=f"Sample {sample_b_id} not found")
    return sample_service.compare_samples(db, sample_a, sample_b)


@router.get("/hotspots/{reference_name}", response_model=schemas.HotspotAnalysisOut)
def get_hotspots(
    reference_name: str,
    threshold_per_100bp: float = Query(3.0, ge=0.1, description="Variants per 100bp threshold for hotspot detection"),
    db: Session = Depends(get_db),
):
    return sample_service.find_hotspots(
        db,
        reference_name,
        threshold_per_100bp=threshold_per_100bp,
    )
