from fastapi import APIRouter, Depends, HTTPException, Query
from sqlalchemy.orm import Session
from typing import Optional
from app import schemas
from app.database import get_db
from app.services import phylogeny_service

router = APIRouter(prefix="/api/phylogeny", tags=["phylogeny"])


@router.post("/distance-matrix", response_model=schemas.DistanceMatrixOut)
def compute_distance_matrix(
    request: schemas.DistanceMatrixRequest,
    db: Session = Depends(get_db),
):
    try:
        return phylogeny_service.compute_distance_matrix_endpoint(db, request)
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))


@router.post("/tree", response_model=schemas.PhyloTreeOut)
def build_phylogenetic_tree(
    request: schemas.PhyloTreeRequest,
    db: Session = Depends(get_db),
):
    try:
        return phylogeny_service.create_phylogeny_task(db, request)
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))


@router.get("/tree/{task_id}", response_model=schemas.PhyloTreeOut)
def get_tree_result(
    task_id: str,
    db: Session = Depends(get_db),
):
    try:
        result = phylogeny_service.get_task_result(db, task_id)
        if result.is_stale:
            import warnings
            warnings.warn("Tree result may be stale - underlying sample data has changed. Consider rebuilding.")
        return result
    except ValueError as e:
        raise HTTPException(status_code=404, detail=str(e))


@router.get("/tree/{task_id}/status", response_model=schemas.PhyloTaskStatusOut)
def get_tree_task_status(
    task_id: str,
    db: Session = Depends(get_db),
):
    try:
        return phylogeny_service.get_task_status(db, task_id)
    except ValueError as e:
        raise HTTPException(status_code=404, detail=str(e))


@router.post("/tree/compare", response_model=schemas.TreeCompareOut)
def compare_phylogenetic_trees(
    request: schemas.TreeCompareRequest,
    db: Session = Depends(get_db),
):
    try:
        return phylogeny_service.compare_trees_endpoint(db, request)
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))
