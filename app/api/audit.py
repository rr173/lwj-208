from fastapi import APIRouter, Depends, Query
from sqlalchemy.orm import Session
from typing import Optional
from datetime import datetime
from app import schemas
from app.database import get_db
from app.audit import get_audit_logs

router = APIRouter(prefix="/api/audit", tags=["audit_logs"])


@router.get("/logs", response_model=schemas.AuditLogListOut)
def list_audit_logs(
    start_time: Optional[datetime] = Query(None, description="开始时间"),
    end_time: Optional[datetime] = Query(None, description="结束时间"),
    resource_type: Optional[str] = Query(None, description="资源类型"),
    operation_type: Optional[str] = Query(None, description="操作类型 (create/update/delete)"),
    page: int = Query(1, ge=1, description="页码"),
    page_size: int = Query(20, ge=1, le=100, description="每页数量"),
    db: Session = Depends(get_db),
):
    result = get_audit_logs(
        db,
        start_time=start_time,
        end_time=end_time,
        resource_type=resource_type,
        operation_type=operation_type,
        page=page,
        page_size=page_size,
    )
    return result


@router.get("/logs/{log_id}", response_model=schemas.AuditLogOut)
def get_audit_log(
    log_id: int,
    db: Session = Depends(get_db),
):
    from app import models
    log = db.query(models.AuditLog).filter(models.AuditLog.id == log_id).first()
    if not log:
        from fastapi import HTTPException
        raise HTTPException(status_code=404, detail="Audit log not found")
    return log
