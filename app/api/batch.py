from fastapi import APIRouter, Depends, HTTPException
from sqlalchemy.orm import Session
from typing import List
from app import schemas
from app.database import get_db
from app.services import genome_service, batch_service

router = APIRouter(prefix="/api/batch", tags=["batch"])


def alignment_to_out(alignment) -> schemas.AlignmentResultOut:
    return schemas.AlignmentResultOut(
        id=alignment.id,
        reference_name=alignment.reference.name if alignment.reference else "",
        score=alignment.score,
        ref_start=alignment.ref_start,
        ref_end=alignment.ref_end,
        query_start=alignment.query_start,
        query_end=alignment.query_end,
        cigar=alignment.cigar,
        alignment=schemas.AlignmentDetail(
            query=alignment.alignment_query,
            match=alignment.alignment_match,
            reference=alignment.alignment_ref,
        ),
        variants=alignment.variants,
    )


@router.post("", response_model=schemas.BatchTaskOut)
def create_batch_task(request: schemas.BatchAlignRequest, db: Session = Depends(get_db)):
    if len(request.query_sequences) > 50:
        raise HTTPException(status_code=400, detail="Maximum 50 query sequences per batch")

    for seq in request.query_sequences:
        if len(seq) > 1000:
            raise HTTPException(status_code=400, detail="Query sequence exceeds 1000 base limit")

    task = batch_service.create_batch_task(db, request.query_sequences)
    batch_service.start_batch_task(task.task_id)

    return task


@router.get("/{task_id}", response_model=schemas.BatchTaskResultOut)
def get_batch_task(task_id: str, db: Session = Depends(get_db)):
    task = batch_service.get_batch_task(db, task_id)
    if not task:
        raise HTTPException(status_code=404, detail="Batch task not found")

    results = []
    for item in task.items:
        if item.alignment:
            results.append(alignment_to_out(item.alignment))

    return schemas.BatchTaskResultOut(
        task_id=task.task_id,
        status=task.status,
        total=task.total,
        completed=task.completed,
        created_at=task.created_at,
        completed_at=task.completed_at,
        results=results,
    )
