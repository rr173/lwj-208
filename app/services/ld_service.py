import hashlib
import json
from typing import List, Dict, Optional, Tuple
from datetime import datetime
from sqlalchemy.orm import Session
from sqlalchemy import func

from app import models, schemas
from app.ld.ld_calculator import compute_ld_matrix
from app.ld.haplotype_blocks import find_haplotype_blocks


def _compute_data_hash(
    db: Session,
    reference_id: int
) -> str:
    """
    Compute a hash of all SNP data across all samples for a given reference.
    Used to detect when underlying data has changed.
    """
    all_snps = db.query(
        models.SampleVariantSpectrum.sample_id,
        models.SampleVariantSpectrum.ref_pos,
        models.SampleVariantSpectrum.ref_base,
        models.SampleVariantSpectrum.alt_base,
    ).filter(
        models.SampleVariantSpectrum.reference_id == reference_id,
        models.SampleVariantSpectrum.variant_type == "SNP",
    ).order_by(
        models.SampleVariantSpectrum.sample_id,
        models.SampleVariantSpectrum.ref_pos,
    ).all()

    hash_input = json.dumps([
        (s.sample_id, s.ref_pos, s.ref_base, s.alt_base)
        for s in all_snps
    ], sort_keys=True)

    return hashlib.md5(hash_input.encode()).hexdigest()


def _get_snp_data_for_reference(
    db: Session,
    reference_id: int
) -> Tuple[List[int], List[int], Dict[int, Dict[int, bool]]]:
    """
    Extract all SNP data for a given reference across all samples.
    
    Returns:
        Tuple of (sample_ids, snp_positions, genotypes_by_sample)
        - sample_ids: sorted list of all sample IDs
        - snp_positions: sorted list of all SNP positions
        - genotypes_by_sample: dict mapping sample_id -> {pos: is_variant}
    """
    samples = db.query(models.Sample).order_by(models.Sample.id).all()
    sample_ids = [s.id for s in samples]

    all_snps = db.query(models.SampleVariantSpectrum).filter(
        models.SampleVariantSpectrum.reference_id == reference_id,
        models.SampleVariantSpectrum.variant_type == "SNP",
    ).all()

    snp_positions_set = set()
    genotypes_by_sample: Dict[int, Dict[int, bool]] = {sid: {} for sid in sample_ids}

    for snp in all_snps:
        snp_positions_set.add(snp.ref_pos)
        if snp.sample_id in genotypes_by_sample:
            genotypes_by_sample[snp.sample_id][snp.ref_pos] = True

    snp_positions = sorted(snp_positions_set)

    return sample_ids, snp_positions, genotypes_by_sample


def _check_ld_matrix_cache(
    db: Session,
    reference_name: str,
    data_hash: str
) -> Optional[models.LDMatrixResult]:
    """
    Check if a valid cached LD matrix result exists.
    Returns the cached result if valid, None otherwise.
    """
    cache_key = hashlib.md5(f"ld_matrix:{reference_name}".encode()).hexdigest()

    result = db.query(models.LDMatrixResult).filter(
        models.LDMatrixResult.cache_key == cache_key
    ).first()

    if not result:
        return None

    result.last_checked_at = func.now()

    if result.is_stale or result.data_hash != data_hash:
        db.delete(result)
        db.commit()
        return None

    db.commit()
    return result


def _save_ld_matrix_result(
    db: Session,
    reference_id: int,
    reference_name: str,
    total_samples: int,
    snp_positions: List[int],
    ld_pairs: List[Dict],
    data_hash: str
) -> models.LDMatrixResult:
    """
    Save LD matrix result to database.
    """
    cache_key = hashlib.md5(f"ld_matrix:{reference_name}".encode()).hexdigest()

    existing = db.query(models.LDMatrixResult).filter(
        models.LDMatrixResult.cache_key == cache_key
    ).first()
    if existing:
        db.delete(existing)
        db.commit()

    result = models.LDMatrixResult(
        cache_key=cache_key,
        reference_name=reference_name,
        reference_id=reference_id,
        total_samples=total_samples,
        snp_count=len(snp_positions),
        snp_positions=snp_positions,
        ld_matrix=ld_pairs,
        data_hash=data_hash,
        is_stale=False,
    )
    db.add(result)
    db.commit()
    db.refresh(result)
    return result


def compute_ld_matrix_endpoint(
    db: Session,
    reference_name: str
) -> schemas.LDMatrixOut:
    """
    Compute or retrieve cached LD matrix for a reference sequence.
    
    Args:
        db: Database session
        reference_name: Name of the reference sequence
    
    Returns:
        LDMatrixOut schema with LD pairs
    
    Raises:
        ValueError: if reference not found, fewer than 5 samples, or no SNPs
    """
    ref = db.query(models.ReferenceSequence).filter(
        models.ReferenceSequence.name == reference_name
    ).first()
    if not ref:
        raise ValueError(f"Reference sequence '{reference_name}' not found")

    total_samples = db.query(func.count(models.Sample.id)).scalar() or 0
    if total_samples < 5:
        raise ValueError(
            f"LD analysis requires at least 5 samples. "
            f"Currently only {total_samples} sample(s) available."
        )

    data_hash = _compute_data_hash(db, ref.id)

    cached = _check_ld_matrix_cache(db, reference_name, data_hash)
    if cached:
        ld_pairs = [
            schemas.LDPair(
                pos_i=p["pos_i"],
                pos_j=p["pos_j"],
                d=p["d"],
                d_prime=p["d_prime"],
                r_squared=p["r_squared"],
            )
            for p in cached.ld_matrix
        ]
        return schemas.LDMatrixOut(
            reference_name=reference_name,
            total_samples=cached.total_samples,
            snp_count=cached.snp_count,
            snp_positions=cached.snp_positions,
            ld_pairs=ld_pairs,
            is_stale=cached.is_stale,
            created_at=cached.created_at,
        )

    sample_ids, snp_positions, genotypes_by_sample = _get_snp_data_for_reference(db, ref.id)

    if len(snp_positions) < 2:
        raise ValueError(
            f"LD analysis requires at least 2 SNP loci. "
            f"Found only {len(snp_positions)} SNP(s) on reference '{reference_name}'."
        )

    ld_pairs_raw = compute_ld_matrix(sample_ids, snp_positions, genotypes_by_sample)

    _save_ld_matrix_result(
        db, ref.id, reference_name, len(sample_ids), snp_positions, ld_pairs_raw, data_hash
    )

    ld_pairs = [
        schemas.LDPair(
            pos_i=p["pos_i"],
            pos_j=p["pos_j"],
            d=p["d"],
            d_prime=p["d_prime"],
            r_squared=p["r_squared"],
        )
        for p in ld_pairs_raw
    ]

    return schemas.LDMatrixOut(
        reference_name=reference_name,
        total_samples=len(sample_ids),
        snp_count=len(snp_positions),
        snp_positions=snp_positions,
        ld_pairs=ld_pairs,
        is_stale=False,
        created_at=datetime.utcnow(),
    )


def _check_haplotype_block_cache(
    db: Session,
    reference_name: str,
    r2_threshold: float,
    data_hash: str
) -> Optional[models.HaplotypeBlockResult]:
    """
    Check if a valid cached haplotype block result exists.
    """
    cache_key = hashlib.md5(
        f"haplotype_blocks:{reference_name}:{r2_threshold}".encode()
    ).hexdigest()

    result = db.query(models.HaplotypeBlockResult).filter(
        models.HaplotypeBlockResult.cache_key == cache_key
    ).first()

    if not result:
        return None

    result.last_checked_at = func.now()

    if result.is_stale or result.data_hash != data_hash:
        db.delete(result)
        db.commit()
        return None

    db.commit()
    return result


def _save_haplotype_block_result(
    db: Session,
    reference_id: int,
    reference_name: str,
    r2_threshold: float,
    total_samples: int,
    snp_count: int,
    blocks: List[Dict],
    data_hash: str
) -> models.HaplotypeBlockResult:
    """
    Save haplotype block result to database.
    """
    cache_key = hashlib.md5(
        f"haplotype_blocks:{reference_name}:{r2_threshold}".encode()
    ).hexdigest()

    existing = db.query(models.HaplotypeBlockResult).filter(
        models.HaplotypeBlockResult.cache_key == cache_key
    ).first()
    if existing:
        db.delete(existing)
        db.commit()

    result = models.HaplotypeBlockResult(
        cache_key=cache_key,
        reference_name=reference_name,
        reference_id=reference_id,
        r2_threshold=r2_threshold,
        total_samples=total_samples,
        snp_count=snp_count,
        block_count=len(blocks),
        blocks=blocks,
        data_hash=data_hash,
        is_stale=False,
    )
    db.add(result)
    db.commit()
    db.refresh(result)
    return result


def _get_ld_matrix_data(
    db: Session,
    reference_id: int,
    reference_name: str
) -> Tuple[List[int], List[Dict]]:
    """
    Get LD matrix data, either from cache or by computing it.
    """
    data_hash = _compute_data_hash(db, reference_id)

    cached = _check_ld_matrix_cache(db, reference_name, data_hash)
    if cached:
        return cached.snp_positions, cached.ld_matrix

    sample_ids, snp_positions, genotypes_by_sample = _get_snp_data_for_reference(db, reference_id)
    ld_pairs_raw = compute_ld_matrix(sample_ids, snp_positions, genotypes_by_sample)

    _save_ld_matrix_result(
        db, reference_id, reference_name, len(sample_ids), snp_positions, ld_pairs_raw, data_hash
    )

    return snp_positions, ld_pairs_raw


def compute_haplotype_blocks_endpoint(
    db: Session,
    reference_name: str,
    r2_threshold: float = 0.8
) -> schemas.HaplotypeBlockOut:
    """
    Compute or retrieve cached haplotype blocks for a reference sequence.
    
    Args:
        db: Database session
        reference_name: Name of the reference sequence
        r2_threshold: Minimum r² threshold for block membership (default 0.8)
    
    Returns:
        HaplotypeBlockOut schema with block information
    
    Raises:
        ValueError: if reference not found, fewer than 5 samples, or no SNPs
    """
    if r2_threshold < 0 or r2_threshold > 1:
        raise ValueError("r² threshold must be between 0 and 1")

    ref = db.query(models.ReferenceSequence).filter(
        models.ReferenceSequence.name == reference_name
    ).first()
    if not ref:
        raise ValueError(f"Reference sequence '{reference_name}' not found")

    total_samples = db.query(func.count(models.Sample.id)).scalar() or 0
    if total_samples < 5:
        raise ValueError(
            f"LD analysis requires at least 5 samples. "
            f"Currently only {total_samples} sample(s) available."
        )

    data_hash = _compute_data_hash(db, ref.id)

    cached = _check_haplotype_block_cache(db, reference_name, r2_threshold, data_hash)
    if cached:
        blocks = [
            schemas.HaplotypeBlock(
                block_index=b["block_index"],
                start_pos=b["start_pos"],
                end_pos=b["end_pos"],
                snp_count=b["snp_count"],
                snp_positions=b["snp_positions"],
                avg_r_squared=b["avg_r_squared"],
            )
            for b in cached.blocks
        ]
        return schemas.HaplotypeBlockOut(
            reference_name=reference_name,
            r2_threshold=r2_threshold,
            total_samples=cached.total_samples,
            snp_count=cached.snp_count,
            block_count=cached.block_count,
            blocks=blocks,
            is_stale=cached.is_stale,
            created_at=cached.created_at,
        )

    snp_positions, ld_pairs_raw = _get_ld_matrix_data(db, ref.id, reference_name)

    if len(snp_positions) < 2:
        raise ValueError(
            f"Haplotype block analysis requires at least 2 SNP loci. "
            f"Found only {len(snp_positions)} SNP(s) on reference '{reference_name}'."
        )

    blocks_raw = find_haplotype_blocks(snp_positions, ld_pairs_raw, r2_threshold)

    total_samples_actual = db.query(func.count(models.Sample.id)).scalar() or 0

    _save_haplotype_block_result(
        db, ref.id, reference_name, r2_threshold,
        total_samples_actual, len(snp_positions), blocks_raw, data_hash
    )

    blocks = [
        schemas.HaplotypeBlock(
            block_index=b["block_index"],
            start_pos=b["start_pos"],
            end_pos=b["end_pos"],
            snp_count=b["snp_count"],
            snp_positions=b["snp_positions"],
            avg_r_squared=b["avg_r_squared"],
        )
        for b in blocks_raw
    ]

    return schemas.HaplotypeBlockOut(
        reference_name=reference_name,
        r2_threshold=r2_threshold,
        total_samples=total_samples_actual,
        snp_count=len(snp_positions),
        block_count=len(blocks),
        blocks=blocks,
        is_stale=False,
        created_at=datetime.utcnow(),
    )


def invalidate_ld_results_for_reference(
    db: Session,
    reference_id: int
):
    """
    Mark all LD and haplotype block results for a given reference as stale.
    """
    ld_results = db.query(models.LDMatrixResult).filter(
        models.LDMatrixResult.reference_id == reference_id,
        models.LDMatrixResult.is_stale == False,
    ).all()
    for r in ld_results:
        r.is_stale = True

    block_results = db.query(models.HaplotypeBlockResult).filter(
        models.HaplotypeBlockResult.reference_id == reference_id,
        models.HaplotypeBlockResult.is_stale == False,
    ).all()
    for r in block_results:
        r.is_stale = True

    db.commit()


def invalidate_ld_results_for_sample(
    db: Session,
    sample_id: int
):
    """
    Mark all LD and haplotype block results as stale if they involve a reference
    that the given sample has variants on.
    """
    ref_ids = db.query(
        models.SampleVariantSpectrum.reference_id
    ).filter(
        models.SampleVariantSpectrum.sample_id == sample_id
    ).distinct().all()

    for (ref_id,) in ref_ids:
        invalidate_ld_results_for_reference(db, ref_id)
