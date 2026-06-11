from fastapi import APIRouter, Depends, HTTPException
from sqlalchemy.orm import Session
from app import schemas
from app.database import get_db
from app.services import genome_service

router = APIRouter(prefix="/api/stats", tags=["statistics"])


@router.get("/{reference_name}/variants", response_model=schemas.VariantStatsOut)
def get_variant_stats(reference_name: str, bucket_size: int = 1000, db: Session = Depends(get_db)):
    ref = genome_service.get_reference_by_name(db, reference_name)
    if not ref:
        raise HTTPException(status_code=404, detail="Reference sequence not found")

    stats = genome_service.get_variant_stats(db, ref.id, bucket_size)
    return stats
