import hashlib
import json
from typing import List, Dict, Optional, Tuple
from datetime import datetime
from sqlalchemy.orm import Session
from sqlalchemy import func

from app import models, schemas
from app.ld.ld_calculator import compute_ld_matrix
from app.ld.haplotype_blocks import find_haplotype_blocks
from app.services import qc_service


def _compute_snp_maf(
    pos: int,
    sample_ids: List[int],
    genotypes_by_sample: Dict[int, Dict[int, bool]]
) -> float:
    """
    Compute Minor Allele Frequency (MAF) for a SNP position.
    
    MAF = min(variant_frequency, 1 - variant_frequency)
    
    Args:
        pos: SNP position
        sample_ids: List of sample IDs
        genotypes_by_sample: Dict mapping sample_id -> {pos: is_variant}
    
    Returns:
        MAF value between 0 and 0.5
    """
    n_samples = len(sample_ids)
    if n_samples == 0:
        return 0.0
    
    variant_count = 0
    for sid in sample_ids:
        sample_geno = genotypes_by_sample.get(sid, {})
        if sample_geno.get(pos, False):
            variant_count += 1
    
    variant_freq = variant_count / n_samples
    return min(variant_freq, 1 - variant_freq)


def _filter_snps(
    sample_ids: List[int],
    snp_positions: List[int],
    genotypes_by_sample: Dict[int, Dict[int, bool]],
    start: Optional[int] = None,
    end: Optional[int] = None,
    min_maf: float = 0.0
) -> Tuple[List[int], int]:
    """
    Filter SNP positions by region and MAF.
    
    Args:
        sample_ids: List of sample IDs
        snp_positions: Original sorted list of SNP positions
        genotypes_by_sample: Dict mapping sample_id -> {pos: is_variant}
        start: Optional start position (1-based, inclusive)
        end: Optional end position (1-based, inclusive)
        min_maf: Minimum minor allele frequency threshold
    
    Returns:
        Tuple of (filtered_snp_positions, filtered_count)
        - filtered_snp_positions: Sorted list of SNP positions that pass all filters
        - filtered_count: Number of SNPs excluded by filtering
    """
    original_count = len(snp_positions)
    filtered_positions = []
    
    for pos in snp_positions:
        if start is not None and pos < start:
            continue
        if end is not None and pos > end:
            continue
        
        if min_maf > 0:
            maf = _compute_snp_maf(pos, sample_ids, genotypes_by_sample)
            if maf < min_maf:
                continue
        
        filtered_positions.append(pos)
    
    filtered_count = original_count - len(filtered_positions)
    return filtered_positions, filtered_count


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
    reference_id: int,
    include_failed_qc: bool = False
) -> Tuple[List[int], List[int], Dict[int, Dict[int, bool]]]:
    samples = db.query(models.Sample).order_by(models.Sample.id).all()
    sample_ids = [s.id for s in samples]

    all_snps_raw = db.query(models.SampleVariantSpectrum).filter(
        models.SampleVariantSpectrum.reference_id == reference_id,
        models.SampleVariantSpectrum.variant_type == "SNP",
    ).all()

    if not include_failed_qc:
        snps_by_sample: Dict[int, List[models.SampleVariantSpectrum]] = {}
        for snp in all_snps_raw:
            snps_by_sample.setdefault(snp.sample_id, []).append(snp)

        all_snps = []
        for sid, snp_list in snps_by_sample.items():
            filtered = qc_service.filter_spectrum_by_qc(
                db, snp_list, sid, include_failed_qc=False
            )
            all_snps.extend(filtered)
    else:
        all_snps = all_snps_raw

    snp_positions_set = set()
    genotypes_by_sample: Dict[int, Dict[int, bool]] = {sid: {} for sid in sample_ids}

    for snp in all_snps:
        snp_positions_set.add(snp.ref_pos)
        if snp.sample_id in genotypes_by_sample:
            genotypes_by_sample[snp.sample_id][snp.ref_pos] = True

    snp_positions = sorted(snp_positions_set)

    return sample_ids, snp_positions, genotypes_by_sample


def _build_ld_cache_key(
    reference_name: str,
    start: Optional[int] = None,
    end: Optional[int] = None,
    min_maf: float = 0.0,
    include_failed_qc: bool = False
) -> str:
    key_parts = [
        "ld_matrix",
        reference_name,
        f"start={start}",
        f"end={end}",
        f"min_maf={min_maf}",
        f"include_failed_qc={include_failed_qc}"
    ]
    return hashlib.md5(":".join(key_parts).encode()).hexdigest()


def _check_ld_matrix_cache(
    db: Session,
    reference_name: str,
    data_hash: str,
    start: Optional[int] = None,
    end: Optional[int] = None,
    min_maf: float = 0.0,
    include_failed_qc: bool = False
) -> Optional[models.LDMatrixResult]:
    cache_key = _build_ld_cache_key(reference_name, start, end, min_maf, include_failed_qc)

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
    data_hash: str,
    start: Optional[int] = None,
    end: Optional[int] = None,
    min_maf: float = 0.0,
    include_failed_qc: bool = False
) -> models.LDMatrixResult:
    cache_key = _build_ld_cache_key(reference_name, start, end, min_maf, include_failed_qc)

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
    reference_name: str,
    start: Optional[int] = None,
    end: Optional[int] = None,
    min_maf: float = 0.0,
    include_failed_qc: bool = False
) -> schemas.LDMatrixOut:
    if min_maf < 0 or min_maf > 0.5:
        raise ValueError("min_maf must be between 0 and 0.5")

    ref = db.query(models.ReferenceSequence).filter(
        models.ReferenceSequence.name == reference_name
    ).first()
    if not ref:
        raise ValueError(f"Reference sequence '{reference_name}' not found")

    if start is not None and start < 1:
        raise ValueError("start must be >= 1")
    if end is not None and end > ref.length:
        raise ValueError(f"end must be <= reference sequence length ({ref.length})")
    if start is not None and end is not None and start > end:
        raise ValueError("start must be <= end")

    total_samples = db.query(func.count(models.Sample.id)).scalar() or 0
    if total_samples < 5:
        raise ValueError(
            f"LD analysis requires at least 5 samples. "
            f"Currently only {total_samples} sample(s) available."
        )

    data_hash = _compute_data_hash(db, ref.id)

    cached = _check_ld_matrix_cache(db, reference_name, data_hash, start, end, min_maf, include_failed_qc)
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
            original_snp_count=cached.snp_count,
            filtered_snp_count=0,
            snp_positions=cached.snp_positions,
            ld_pairs=ld_pairs,
            start=start,
            end=end,
            min_maf=min_maf if min_maf > 0 else None,
            is_stale=cached.is_stale,
            created_at=cached.created_at,
        )

    sample_ids, snp_positions, genotypes_by_sample = _get_snp_data_for_reference(db, ref.id, include_failed_qc=include_failed_qc)
    original_snp_count = len(snp_positions)

    filtered_snps, filtered_count = _filter_snps(
        sample_ids, snp_positions, genotypes_by_sample, start, end, min_maf
    )

    if len(filtered_snps) < 2:
        raise ValueError(
            f"LD analysis requires at least 2 SNP loci after filtering. "
            f"Found only {len(filtered_snps)} SNP(s) on reference '{reference_name}' "
            f"after applying filters (filtered out {filtered_count} SNPs)."
        )

    ld_pairs_raw = compute_ld_matrix(sample_ids, filtered_snps, genotypes_by_sample)

    _save_ld_matrix_result(
        db, ref.id, reference_name, len(sample_ids), filtered_snps,
        ld_pairs_raw, data_hash, start, end, min_maf, include_failed_qc
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
        snp_count=len(filtered_snps),
        original_snp_count=original_snp_count,
        filtered_snp_count=filtered_count,
        snp_positions=filtered_snps,
        ld_pairs=ld_pairs,
        start=start,
        end=end,
        min_maf=min_maf if min_maf > 0 else None,
        is_stale=False,
        created_at=datetime.utcnow(),
    )


def _build_haplotype_cache_key(
    reference_name: str,
    r2_threshold: float,
    start: Optional[int] = None,
    end: Optional[int] = None,
    min_maf: float = 0.0,
    include_failed_qc: bool = False
) -> str:
    key_parts = [
        "haplotype_blocks",
        reference_name,
        f"r2={r2_threshold}",
        f"start={start}",
        f"end={end}",
        f"min_maf={min_maf}",
        f"include_failed_qc={include_failed_qc}"
    ]
    return hashlib.md5(":".join(key_parts).encode()).hexdigest()


def _check_haplotype_block_cache(
    db: Session,
    reference_name: str,
    r2_threshold: float,
    data_hash: str,
    start: Optional[int] = None,
    end: Optional[int] = None,
    min_maf: float = 0.0,
    include_failed_qc: bool = False
) -> Optional[models.HaplotypeBlockResult]:
    cache_key = _build_haplotype_cache_key(
        reference_name, r2_threshold, start, end, min_maf, include_failed_qc
    )

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
    data_hash: str,
    start: Optional[int] = None,
    end: Optional[int] = None,
    min_maf: float = 0.0,
    include_failed_qc: bool = False
) -> models.HaplotypeBlockResult:
    cache_key = _build_haplotype_cache_key(
        reference_name, r2_threshold, start, end, min_maf, include_failed_qc
    )

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
    reference_name: str,
    start: Optional[int] = None,
    end: Optional[int] = None,
    min_maf: float = 0.0,
    include_failed_qc: bool = False
) -> Tuple[List[int], List[Dict], int, int]:
    data_hash = _compute_data_hash(db, reference_id)

    cached = _check_ld_matrix_cache(db, reference_name, data_hash, start, end, min_maf, include_failed_qc)
    if cached:
        return cached.snp_positions, cached.ld_matrix, cached.snp_count, 0

    sample_ids, snp_positions, genotypes_by_sample = _get_snp_data_for_reference(db, reference_id, include_failed_qc=include_failed_qc)
    original_snp_count = len(snp_positions)

    filtered_snps, filtered_count = _filter_snps(
        sample_ids, snp_positions, genotypes_by_sample, start, end, min_maf
    )

    ld_pairs_raw = compute_ld_matrix(sample_ids, filtered_snps, genotypes_by_sample)

    _save_ld_matrix_result(
        db, reference_id, reference_name, len(sample_ids), filtered_snps,
        ld_pairs_raw, data_hash, start, end, min_maf, include_failed_qc
    )

    return filtered_snps, ld_pairs_raw, original_snp_count, filtered_count


def compute_haplotype_blocks_endpoint(
    db: Session,
    reference_name: str,
    r2_threshold: float = 0.8,
    start: Optional[int] = None,
    end: Optional[int] = None,
    min_maf: float = 0.0,
    include_failed_qc: bool = False
) -> schemas.HaplotypeBlockOut:
    if r2_threshold < 0 or r2_threshold > 1:
        raise ValueError("r² threshold must be between 0 and 1")
    if min_maf < 0 or min_maf > 0.5:
        raise ValueError("min_maf must be between 0 and 0.5")

    ref = db.query(models.ReferenceSequence).filter(
        models.ReferenceSequence.name == reference_name
    ).first()
    if not ref:
        raise ValueError(f"Reference sequence '{reference_name}' not found")

    if start is not None and start < 1:
        raise ValueError("start must be >= 1")
    if end is not None and end > ref.length:
        raise ValueError(f"end must be <= reference sequence length ({ref.length})")
    if start is not None and end is not None and start > end:
        raise ValueError("start must be <= end")

    total_samples = db.query(func.count(models.Sample.id)).scalar() or 0
    if total_samples < 5:
        raise ValueError(
            f"LD analysis requires at least 5 samples. "
            f"Currently only {total_samples} sample(s) available."
        )

    data_hash = _compute_data_hash(db, ref.id)

    cached = _check_haplotype_block_cache(
        db, reference_name, r2_threshold, data_hash, start, end, min_maf, include_failed_qc
    )
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
            original_snp_count=cached.snp_count,
            filtered_snp_count=0,
            block_count=cached.block_count,
            blocks=blocks,
            start=start,
            end=end,
            min_maf=min_maf if min_maf > 0 else None,
            is_stale=cached.is_stale,
            created_at=cached.created_at,
        )

    snp_positions, ld_pairs_raw, original_snp_count, filtered_count = _get_ld_matrix_data(
        db, ref.id, reference_name, start, end, min_maf, include_failed_qc
    )

    if len(snp_positions) < 2:
        raise ValueError(
            f"Haplotype block analysis requires at least 2 SNP loci after filtering. "
            f"Found only {len(snp_positions)} SNP(s) on reference '{reference_name}' "
            f"after applying filters (filtered out {filtered_count} SNPs)."
        )

    blocks_raw = find_haplotype_blocks(snp_positions, ld_pairs_raw, r2_threshold)

    total_samples_actual = db.query(func.count(models.Sample.id)).scalar() or 0

    _save_haplotype_block_result(
        db, ref.id, reference_name, r2_threshold,
        total_samples_actual, len(snp_positions), blocks_raw, data_hash,
        start, end, min_maf, include_failed_qc
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
        original_snp_count=original_snp_count,
        filtered_snp_count=filtered_count,
        block_count=len(blocks),
        blocks=blocks,
        start=start,
        end=end,
        min_maf=min_maf if min_maf > 0 else None,
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
