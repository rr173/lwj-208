from typing import Optional
from fastapi import APIRouter, Depends, HTTPException
from sqlalchemy.orm import Session
from app import schemas
from app.database import get_db
from app.services import ld_service

router = APIRouter(prefix="/api/ld", tags=["linkage_disequilibrium"])


@router.get("/{reference_name}/matrix", response_model=schemas.LDMatrixOut)
def get_ld_matrix(
    reference_name: str,
    start: Optional[int] = None,
    end: Optional[int] = None,
    min_maf: float = 0.0,
    db: Session = Depends(get_db)
):
    """
    Compute linkage disequilibrium (LD) matrix for SNP pairs on a reference sequence.
    
    Returns pairwise D, D', and r² values for SNP loci across all samples.
    Requires at least 5 samples and 2 SNP loci.
    
    Args:
        reference_name: Name of the reference sequence
        start: Optional start position (1-based) for region filtering.
               Only SNPs at or after this position are included.
        end: Optional end position (1-based) for region filtering.
             Only SNPs at or before this position are included.
        min_maf: Minimum minor allele frequency (default 0, no filtering).
                 SNPs with MAF < min_maf are excluded.
                 MAF = min(variant_frequency, 1 - variant_frequency).
    """
    try:
        result = ld_service.compute_ld_matrix_endpoint(
            db, reference_name, start=start, end=end, min_maf=min_maf
        )
        return result
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))


@router.get("/{reference_name}/haplotype-blocks", response_model=schemas.HaplotypeBlockOut)
def get_haplotype_blocks(
    reference_name: str,
    r2_threshold: float = 0.8,
    start: Optional[int] = None,
    end: Optional[int] = None,
    min_maf: float = 0.0,
    db: Session = Depends(get_db)
):
    """
    Identify haplotype blocks based on linkage disequilibrium.
    
    A haplotype block is a continuous segment where all pairs of SNPs within
    the block have r² >= threshold (not just adjacent pairs).
    
    Args:
        reference_name: Name of the reference sequence
        r2_threshold: Minimum r² value for block membership (0-1, default 0.8)
        start: Optional start position (1-based) for region filtering.
               Only SNPs at or after this position are included.
        end: Optional end position (1-based) for region filtering.
             Only SNPs at or before this position are included.
        min_maf: Minimum minor allele frequency (default 0, no filtering).
                 SNPs with MAF < min_maf are excluded.
                 MAF = min(variant_frequency, 1 - variant_frequency).
    """
    try:
        result = ld_service.compute_haplotype_blocks_endpoint(
            db, reference_name, r2_threshold, start=start, end=end, min_maf=min_maf
        )
        return result
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))
