from typing import Optional
from fastapi import APIRouter, Depends, HTTPException, Query
from sqlalchemy.orm import Session

from app import schemas
from app.database import get_db
from app.services import pipeline_service

router = APIRouter(prefix="/api/pipelines", tags=["pipelines"])


@router.post("/templates", response_model=schemas.PipelineTemplateOut, status_code=201)
def create_pipeline_template(
    request: schemas.PipelineTemplateCreate,
    db: Session = Depends(get_db),
):
    try:
        existing = pipeline_service.get_pipeline_template_by_name(db, request.name)
        if existing:
            raise HTTPException(
                status_code=400,
                detail=f"Pipeline template name '{request.name}' already exists"
            )
        return pipeline_service.create_pipeline_template(db, request)
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))


@router.get("/templates", response_model=list[schemas.PipelineTemplateOut])
def list_pipeline_templates(
    skip: int = Query(0, ge=0),
    limit: int = Query(100, ge=1, le=1000),
    db: Session = Depends(get_db),
):
    return pipeline_service.list_pipeline_templates(db, skip=skip, limit=limit)


@router.get("/templates/{template_id}", response_model=schemas.PipelineTemplateOut)
def get_pipeline_template(
    template_id: int,
    db: Session = Depends(get_db),
):
    template = pipeline_service.get_pipeline_template(db, template_id)
    if not template:
        raise HTTPException(status_code=404, detail="Pipeline template not found")
    return template


@router.put("/templates/{template_id}", response_model=schemas.PipelineTemplateOut)
def update_pipeline_template(
    template_id: int,
    request: schemas.PipelineTemplateUpdate,
    db: Session = Depends(get_db),
):
    template = pipeline_service.get_pipeline_template(db, template_id)
    if not template:
        raise HTTPException(status_code=404, detail="Pipeline template not found")
    try:
        return pipeline_service.update_pipeline_template(db, template, request)
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))


@router.delete("/templates/{template_id}", status_code=204)
def delete_pipeline_template(
    template_id: int,
    db: Session = Depends(get_db),
):
    template = pipeline_service.get_pipeline_template(db, template_id)
    if not template:
        raise HTTPException(status_code=404, detail="Pipeline template not found")
    pipeline_service.delete_pipeline_template(db, template)
    return None


@router.post("/templates/{template_id}/execute", response_model=schemas.PipelineExecutionOut)
def execute_pipeline(
    template_id: int,
    request: schemas.PipelineExecutionTrigger,
    db: Session = Depends(get_db),
):
    template = pipeline_service.get_pipeline_template(db, template_id)
    if not template:
        raise HTTPException(status_code=404, detail="Pipeline template not found")
    try:
        execution = pipeline_service.execute_pipeline(
            db, template, request.initial_params
        )
        result = schemas.PipelineExecutionOut(
            id=execution.id,
            template_id=execution.template_id,
            template_name=execution.template.name,
            status=execution.status,
            initial_params=execution.initial_params,
            total_duration_seconds=execution.total_duration_seconds,
            error_message=execution.error_message,
            triggered_at=execution.triggered_at,
            completed_at=execution.completed_at,
            template_snapshot=execution.template_snapshot,
            resume_count=execution.resume_count or 0,
            step_executions=[
                schemas.PipelineStepExecutionOut(
                    step_index=se.step_index,
                    step_name=se.step_name,
                    step_type=se.step_type,
                    status=se.status,
                    input_params=se.input_params,
                    output_summary=se.output_summary,
                    error_message=se.error_message,
                    started_at=se.started_at,
                    completed_at=se.completed_at,
                    duration_seconds=se.duration_seconds,
                    retry_count=se.retry_count or 0,
                    last_error=se.last_error,
                )
                for se in execution.step_executions
            ],
        )
        return result
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))


@router.post("/executions/{execution_id}/resume", response_model=schemas.PipelineResumeResponse)
def resume_pipeline_execution(
    execution_id: int,
    db: Session = Depends(get_db),
):
    execution = pipeline_service.get_pipeline_execution(db, execution_id)
    if not execution:
        raise HTTPException(status_code=404, detail="Pipeline execution not found")
    try:
        resumed_execution, skipped_count, resumed_from = pipeline_service.resume_pipeline_execution(
            db, execution
        )
        result = schemas.PipelineResumeResponse(
            id=resumed_execution.id,
            template_id=resumed_execution.template_id,
            template_name=resumed_execution.template.name,
            status=resumed_execution.status,
            initial_params=resumed_execution.initial_params,
            total_duration_seconds=resumed_execution.total_duration_seconds,
            error_message=resumed_execution.error_message,
            triggered_at=resumed_execution.triggered_at,
            completed_at=resumed_execution.completed_at,
            resume_count=resumed_execution.resume_count or 0,
            resumed_from_step=resumed_from,
            skipped_completed_steps=skipped_count,
            step_executions=[
                schemas.PipelineStepExecutionOut(
                    step_index=se.step_index,
                    step_name=se.step_name,
                    step_type=se.step_type,
                    status=se.status,
                    input_params=se.input_params,
                    output_summary=se.output_summary,
                    error_message=se.error_message,
                    started_at=se.started_at,
                    completed_at=se.completed_at,
                    duration_seconds=se.duration_seconds,
                    retry_count=se.retry_count or 0,
                    last_error=se.last_error,
                )
                for se in resumed_execution.step_executions
            ],
        )
        return result
    except pipeline_service.PipelineResumeError as e:
        raise HTTPException(status_code=409, detail=str(e))
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))


@router.get("/templates/{template_id}/executions", response_model=dict)
def get_pipeline_executions(
    template_id: int,
    status: Optional[str] = Query(None, description="Filter by execution status"),
    skip: int = Query(0, ge=0),
    limit: int = Query(100, ge=1, le=1000),
    db: Session = Depends(get_db),
):
    template = pipeline_service.get_pipeline_template(db, template_id)
    if not template:
        raise HTTPException(status_code=404, detail="Pipeline template not found")
    try:
        executions, total = pipeline_service.get_pipeline_executions(
            db, template_id, status=status, skip=skip, limit=limit
        )
        return {
            "total": total,
            "skip": skip,
            "limit": limit,
            "items": [
                schemas.PipelineExecutionSummaryOut(
                    id=e.id,
                    template_id=e.template_id,
                    template_name=template.name,
                    status=e.status,
                    initial_params=e.initial_params,
                    total_duration_seconds=e.total_duration_seconds,
                    error_message=e.error_message,
                    triggered_at=e.triggered_at,
                    completed_at=e.completed_at,
                    resume_count=e.resume_count or 0,
                )
                for e in executions
            ],
        }
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))


@router.get("/executions/{execution_id}", response_model=schemas.PipelineExecutionOut)
def get_pipeline_execution(
    execution_id: int,
    db: Session = Depends(get_db),
):
    execution = pipeline_service.get_pipeline_execution(db, execution_id)
    if not execution:
        raise HTTPException(status_code=404, detail="Pipeline execution not found")
    return schemas.PipelineExecutionOut(
        id=execution.id,
        template_id=execution.template_id,
        template_name=execution.template.name,
        status=execution.status,
        initial_params=execution.initial_params,
        total_duration_seconds=execution.total_duration_seconds,
        error_message=execution.error_message,
        triggered_at=execution.triggered_at,
        completed_at=execution.completed_at,
        template_snapshot=execution.template_snapshot,
        resume_count=execution.resume_count or 0,
        step_executions=[
            schemas.PipelineStepExecutionOut(
                step_index=se.step_index,
                step_name=se.step_name,
                step_type=se.step_type,
                status=se.status,
                input_params=se.input_params,
                output_summary=se.output_summary,
                error_message=se.error_message,
                started_at=se.started_at,
                completed_at=se.completed_at,
                duration_seconds=se.duration_seconds,
                retry_count=se.retry_count or 0,
                last_error=se.last_error,
            )
            for se in execution.step_executions
        ],
    )


@router.get("/step-types", response_model=list[dict])
def list_available_step_types():
    from app.services.pipeline_service import STEP_SERVICE_MAP
    return [
        {
            "step_type": step_type,
            "required_params": info["required_params"],
            "optional_params": info["optional_params"],
        }
        for step_type, info in STEP_SERVICE_MAP.items()
    ]
