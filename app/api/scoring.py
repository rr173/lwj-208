from fastapi import APIRouter, Depends, HTTPException, Query
from sqlalchemy.orm import Session
from typing import List, Optional
from app import schemas
from app.database import get_db
from app.services import scoring_service

router = APIRouter(prefix="/api/scoring", tags=["pathogenicity_scoring"])


@router.get("/rules", response_model=List[schemas.ScoringRuleOut])
def list_rules(
    enabled_only: bool = Query(False, description="Only return enabled rules"),
    db: Session = Depends(get_db),
):
    return scoring_service.list_rules(db, enabled_only=enabled_only)


@router.post("/rules", response_model=schemas.ScoringRuleOut, status_code=201)
def create_rule(request: schemas.ScoringRuleCreate, db: Session = Depends(get_db)):
    return scoring_service.create_rule(db, request)


@router.get("/rules/version")
def get_rule_version(db: Session = Depends(get_db)):
    return {"version": scoring_service.get_current_rule_version(db)}


@router.get("/rules/{rule_id}", response_model=schemas.ScoringRuleOut)
def get_rule(rule_id: int, db: Session = Depends(get_db)):
    rule = scoring_service.get_rule(db, rule_id)
    if not rule:
        raise HTTPException(status_code=404, detail="Rule not found")
    return rule


@router.put("/rules/{rule_id}", response_model=schemas.ScoringRuleOut)
def update_rule(
    rule_id: int,
    request: schemas.ScoringRuleUpdate,
    db: Session = Depends(get_db),
):
    rule = scoring_service.get_rule(db, rule_id)
    if not rule:
        raise HTTPException(status_code=404, detail="Rule not found")
    return scoring_service.update_rule(db, rule, request)


@router.delete("/rules/{rule_id}", status_code=204)
def delete_rule(rule_id: int, db: Session = Depends(get_db)):
    rule = scoring_service.get_rule(db, rule_id)
    if not rule:
        raise HTTPException(status_code=404, detail="Rule not found")
    scoring_service.delete_rule(db, rule)
    return None


@router.post("/samples/{sample_id}/score", response_model=schemas.SampleScoreOut)
def score_sample(sample_id: int, db: Session = Depends(get_db)):
    try:
        return scoring_service.score_sample(db, sample_id)
    except ValueError as e:
        raise HTTPException(status_code=404, detail=str(e))


@router.get("/samples/{sample_id}/score", response_model=schemas.SampleScoreSummaryOut)
def get_sample_score(sample_id: int, db: Session = Depends(get_db)):
    result = scoring_service.get_sample_score(db, sample_id)
    if not result:
        raise HTTPException(status_code=404, detail="Score not found for this sample. Run scoring first.")
    return result


@router.get("/samples/{sample_id}/score/detail", response_model=schemas.SampleScoreOut)
def get_sample_score_detail(sample_id: int, db: Session = Depends(get_db)):
    result = scoring_service.get_sample_score_detail(db, sample_id)
    if not result:
        raise HTTPException(status_code=404, detail="Score not found for this sample. Run scoring first.")
    return result


@router.get("/samples/{sample_id}/variants", response_model=List[schemas.FilteredVariantOut])
def filter_variants(
    sample_id: int,
    classification: Optional[str] = Query(None, description="Filter by classification"),
    gene_name: Optional[str] = Query(None, description="Filter by gene name"),
    db: Session = Depends(get_db),
):
    if classification and classification not in schemas.VALID_CLASSIFICATIONS:
        raise HTTPException(
            status_code=400,
            detail=f"Invalid classification. Valid values: {schemas.VALID_CLASSIFICATIONS}",
        )
    return scoring_service.filter_variants_by_classification(
        db, sample_id, classification=classification, gene_name=gene_name
    )


@router.get("/samples/{sample_id}/stats", response_model=schemas.ScoreStatsOut)
def get_score_stats(sample_id: int, db: Session = Depends(get_db)):
    result = scoring_service.get_score_stats(db, sample_id)
    if not result:
        raise HTTPException(status_code=404, detail="Score not found for this sample. Run scoring first.")
    return result
