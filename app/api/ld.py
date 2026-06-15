from fastapi import APIRouter, Depends, HTTPException
from sqlalchemy.orm import Session
from app import schemas
from app.database import get_db
from app.services import ld_service

router = APIRouter(prefix="/api/ld", tags=["linkage_disequilibrium"])


@router.get("/{reference_name}/matrix", response_model=schemas.LDMatrixOut)
def get_ld_matrix(reference_name: str, db: Session = Depends(get_db)):
    """
    Compute linkage disequilibrium (LD) matrix for all SNP pairs on a reference sequence.
    
    Returns pairwise D, D', and r² values for all SNP loci across all samples.
    Requires at least 5 samples and 2 SNP loci.
    """
    try:
        result = ld_service.compute_ld_matrix_endpoint(db, reference_name)
        return result
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))


@router.get("/{reference_name}/haplotype-blocks", response_model=schemas.HaplotypeBlockOut)
def get_haplotype_blocks(
    reference_name: str,
    r2_threshold: float = 0.8,
    db: Session = Depends(get_db)
):
    """
    Identify haplotype blocks based on linkage disequilibrium.
    
    A haplotype block is a continuous segment where all pairs of SNPs within
    the block have r² >= threshold (not just adjacent pairs).
    
    Args:
        reference_name: Name of the reference sequence
        r2_threshold: Minimum r² value for block membership (0-1, default 0.8)
    """
    try:
        result = ld_service.compute_haplotype_blocks_endpoint(
            db, reference_name, r2_threshold
        )
        return result
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))
