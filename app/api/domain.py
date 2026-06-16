from fastapi import APIRouter, Depends, HTTPException
from sqlalchemy.orm import Session
from typing import List
from app import schemas
from app.database import get_db
from app.services import domain_service

router = APIRouter(prefix="/api/domains", tags=["domains"])


@router.post("", response_model=List[schemas.ProteinDomainOut], status_code=201)
def register_domains(
    request: schemas.ProteinDomainBatchRequest,
    db: Session = Depends(get_db),
):
    for item in request.domains:
        if item.start_codon > item.end_codon:
            raise HTTPException(
                status_code=400,
                detail=f"Domain '{item.domain_name}': start_codon ({item.start_codon}) must not exceed end_codon ({item.end_codon})",
            )
    return domain_service.register_domains(db, request)


@router.get("/{gene_name}", response_model=List[schemas.ProteinDomainOut])
def query_domains_by_gene(gene_name: str, db: Session = Depends(get_db)):
    return domain_service.query_domains_by_gene(db, gene_name)


@router.post("/impact/{sample_id}", response_model=schemas.SampleImpactAssessmentOut)
def assess_variant_impact(sample_id: int, db: Session = Depends(get_db)):
    result = domain_service.map_sample_variants_to_domains(db, sample_id)
    if result is None:
        raise HTTPException(status_code=404, detail="Sample not found")
    return result
