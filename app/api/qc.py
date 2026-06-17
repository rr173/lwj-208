from fastapi import APIRouter, Depends, HTTPException, Query
from sqlalchemy.orm import Session
from typing import List, Optional
from app import schemas
from app.database import get_db
from app.services import qc_service

router = APIRouter(prefix="/api/qc", tags=["quality_control"])


@router.post("/rule-sets", response_model=schemas.QCRuleSetOut, status_code=201)
def create_rule_set(request: schemas.QCRuleSetCreate, db: Session = Depends(get_db)):
    return qc_service.create_rule_set(db, request)


@router.get("/rule-sets", response_model=List[schemas.QCRuleSetListOut])
def list_rule_sets(db: Session = Depends(get_db)):
    rule_sets = qc_service.list_rule_sets(db)
    return [
        schemas.QCRuleSetListOut(
            id=rs.id,
            name=rs.name,
            description=rs.description,
            is_active=rs.is_active,
            rule_count=len(rs.rules),
            created_at=rs.created_at,
            updated_at=rs.updated_at,
        )
        for rs in rule_sets
    ]


@router.get("/rule-sets/active", response_model=schemas.QCRuleSetOut)
def get_active_rule_set(db: Session = Depends(get_db)):
    rule_set = qc_service.get_active_rule_set(db)
    if not rule_set:
        raise HTTPException(status_code=404, detail="No active QC rule set found")
    return rule_set


@router.get("/rule-sets/{rule_set_id}", response_model=schemas.QCRuleSetOut)
def get_rule_set(rule_set_id: int, db: Session = Depends(get_db)):
    rule_set = qc_service.get_rule_set(db, rule_set_id)
    if not rule_set:
        raise HTTPException(status_code=404, detail="Rule set not found")
    return rule_set


@router.put("/rule-sets/{rule_set_id}", response_model=schemas.QCRuleSetOut)
def update_rule_set(
    rule_set_id: int,
    request: schemas.QCRuleSetUpdate,
    db: Session = Depends(get_db),
):
    rule_set = qc_service.get_rule_set(db, rule_set_id)
    if not rule_set:
        raise HTTPException(status_code=404, detail="Rule set not found")
    return qc_service.update_rule_set(db, rule_set, request)


@router.delete("/rule-sets/{rule_set_id}", status_code=204)
def delete_rule_set(rule_set_id: int, db: Session = Depends(get_db)):
    rule_set = qc_service.get_rule_set(db, rule_set_id)
    if not rule_set:
        raise HTTPException(status_code=404, detail="Rule set not found")
    qc_service.delete_rule_set(db, rule_set)
    return None


@router.post("/rule-sets/{rule_set_id}/activate", response_model=schemas.QCRuleSetOut)
def activate_rule_set(rule_set_id: int, db: Session = Depends(get_db)):
    try:
        return qc_service.activate_rule_set(db, rule_set_id)
    except ValueError as e:
        raise HTTPException(status_code=404, detail=str(e))


@router.post("/rule-sets/{rule_set_id}/deactivate", response_model=schemas.QCRuleSetOut)
def deactivate_rule_set(rule_set_id: int, db: Session = Depends(get_db)):
    try:
        return qc_service.deactivate_rule_set(db, rule_set_id)
    except ValueError as e:
        raise HTTPException(status_code=404, detail=str(e))


@router.post("/samples/{sample_id}/execute", response_model=schemas.SampleQCResultOut)
def execute_qc(sample_id: int, db: Session = Depends(get_db)):
    try:
        return qc_service.execute_qc_for_sample(db, sample_id)
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))


@router.get("/samples/{sample_id}/summary", response_model=schemas.SampleQCSummaryOut)
def get_sample_qc_summary(sample_id: int, db: Session = Depends(get_db)):
    try:
        return qc_service.get_sample_qc_summary(db, sample_id)
    except ValueError as e:
        raise HTTPException(status_code=404, detail=str(e))


@router.get("/samples/{sample_id}/results", response_model=List[schemas.VariantQCOut])
def get_sample_qc_results(
    sample_id: int,
    qc_status: Optional[str] = Query(None, description="Filter by QC status (PASS/FAIL)"),
    db: Session = Depends(get_db),
):
    if qc_status and qc_status not in ("PASS", "FAIL"):
        raise HTTPException(status_code=400, detail="qc_status must be 'PASS' or 'FAIL'")
    try:
        return qc_service.get_sample_qc_results(db, sample_id, qc_status=qc_status)
    except ValueError as e:
        raise HTTPException(status_code=404, detail=str(e))


@router.post("/batch/execute", response_model=schemas.BatchQCResultOut)
def execute_batch_qc(
    request: schemas.BatchQCRequest,
    db: Session = Depends(get_db),
):
    try:
        return qc_service.execute_batch_qc(db, request.sample_ids)
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))


@router.get("/report/{reference_name}", response_model=schemas.QCReportOut)
def get_qc_report(
    reference_name: str,
    db: Session = Depends(get_db),
):
    try:
        return qc_service.generate_qc_report_by_reference(db, reference_name)
    except ValueError as e:
        raise HTTPException(status_code=404, detail=str(e))
