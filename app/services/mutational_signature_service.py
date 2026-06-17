from typing import List, Optional, Dict, Tuple
from sqlalchemy.orm import Session
from sqlalchemy import func
from datetime import datetime

from app import models, schemas
from app.mutation_signature.trinucleotide import (
    ALL_96_MUTATION_TYPES,
    get_trinucleotide_context,
    format_mutation_type,
    get_mutation_type_index,
    create_empty_count_matrix,
    compute_data_hash,
    compute_cache_key,
)
from app.mutation_signature.nmf import (
    nmf_with_multiple_initializations,
    find_optimal_k,
)
from app.mutation_signature.matrix_ops import (
    cosine_similarity,
    normalize_vector,
)
from app.mutation_signature.reference_signatures import (
    seed_reference_signatures,
)


def _get_reference_or_error(db: Session, reference_name: str) -> models.ReferenceSequence:
    ref = db.query(models.ReferenceSequence).filter(
        models.ReferenceSequence.name == reference_name
    ).first()
    if not ref:
        raise ValueError(f"Reference sequence '{reference_name}' not found")
    if not ref.sequence:
        raise ValueError(f"Reference sequence '{reference_name}' has no sequence data")
    return ref


def _get_samples_or_error(db: Session, sample_ids: List[int]) -> List[models.Sample]:
    sorted_ids = sorted(sample_ids)
    samples = []
    for sid in sorted_ids:
        sample = db.query(models.Sample).filter(models.Sample.id == sid).first()
        if not sample:
            raise ValueError(f"Sample {sid} not found")
        samples.append(sample)
    return samples


def _build_count_matrix(
    db: Session,
    samples: List[models.Sample],
    reference: models.ReferenceSequence,
) -> Tuple[List[List[int]], List[str], List[str]]:
    n_samples = len(samples)
    count_matrix = create_empty_count_matrix(n_samples)
    ref_sequence = reference.sequence
    ref_len = len(ref_sequence)

    for sample_idx, sample in enumerate(samples):
        spectrum_rows = db.query(models.SampleVariantSpectrum).filter(
            models.SampleVariantSpectrum.sample_id == sample.id,
            models.SampleVariantSpectrum.reference_id == reference.id,
            models.SampleVariantSpectrum.variant_type == "SNP",
        ).all()

        for row in spectrum_rows:
            ref_pos = row.ref_pos
            if ref_pos <= 0 or ref_pos >= ref_len - 1:
                continue

            try:
                left, ref_base, alt_base, right = get_trinucleotide_context(
                    ref_sequence,
                    ref_pos,
                    row.ref_base,
                    row.alt_base,
                )
                mutation_type = format_mutation_type(left, ref_base, alt_base, right)
                type_idx = get_mutation_type_index(mutation_type)
                if type_idx >= 0:
                    count_matrix[type_idx][sample_idx] += 1
            except ValueError:
                continue

    sample_names = [s.name for s in samples]
    return count_matrix, sample_names, ALL_96_MUTATION_TYPES


def _count_matrix_to_float(
    count_matrix: List[List[int]]
) -> List[List[float]]:
    return [[float(x) for x in row] for row in count_matrix]


def build_mutational_spectrum(
    db: Session,
    sample_ids: List[int],
    reference_name: str,
) -> schemas.MutationalSpectrumOut:
    reference = _get_reference_or_error(db, reference_name)
    samples = _get_samples_or_error(db, sample_ids)
    sorted_ids = sorted(sample_ids)

    cache_key = compute_cache_key("spectrum", sorted_ids, reference_name)
    data_hash = compute_data_hash(sorted_ids, reference_name)

    cached = db.query(models.MutationalSpectrumCache).filter(
        models.MutationalSpectrumCache.cache_key == cache_key
    ).first()

    if cached and not cached.is_stale and cached.data_hash == data_hash:
        cached.last_checked_at = func.now()
        db.commit()
        return schemas.MutationalSpectrumOut(
            reference_name=reference_name,
            sample_ids=cached.sample_ids,
            sample_names=cached.sample_names,
            mutation_types=cached.mutation_types,
            count_matrix=cached.count_matrix,
            cached=True,
            is_stale=cached.is_stale,
            created_at=cached.created_at,
        )

    count_matrix, sample_names, mutation_types = _build_count_matrix(
        db, samples, reference
    )

    if cached:
        cached.count_matrix = count_matrix
        cached.sample_names = sample_names
        cached.data_hash = data_hash
        cached.is_stale = False
        cached.last_checked_at = func.now()
        cached.created_at = func.now()
    else:
        new_cache = models.MutationalSpectrumCache(
            cache_key=cache_key,
            sample_ids=sorted_ids,
            reference_name=reference_name,
            reference_id=reference.id,
            sample_names=sample_names,
            mutation_types=mutation_types,
            count_matrix=count_matrix,
            data_hash=data_hash,
        )
        db.add(new_cache)
    db.commit()

    return schemas.MutationalSpectrumOut(
        reference_name=reference_name,
        sample_ids=sorted_ids,
        sample_names=sample_names,
        mutation_types=mutation_types,
        count_matrix=count_matrix,
        cached=False,
        is_stale=False,
        created_at=datetime.now(),
    )


def extract_signatures(
    db: Session,
    sample_ids: List[int],
    reference_name: str,
    k_value: int,
) -> schemas.NMFSignatureOut:
    reference = _get_reference_or_error(db, reference_name)
    samples = _get_samples_or_error(db, sample_ids)
    sorted_ids = sorted(sample_ids)

    cache_key = compute_cache_key("nmf", sorted_ids, reference_name, k_value)
    data_hash = compute_data_hash(sorted_ids, reference_name, k_value)

    cached = db.query(models.NMFSignatureResult).filter(
        models.NMFSignatureResult.cache_key == cache_key
    ).first()

    if cached and not cached.is_stale and cached.data_hash == data_hash:
        cached.last_checked_at = func.now()
        db.commit()
        return _nmf_cache_to_out(cached, reference_name)

    spectrum_out = build_mutational_spectrum(db, sorted_ids, reference_name)
    V = _count_matrix_to_float(spectrum_out.count_matrix)

    W, H, error, iterations = nmf_with_multiple_initializations(
        V,
        k_value,
        n_init=10,
        max_iter=5000,
        tol=1e-6,
        patience=20,
    )

    signature_matrix = [[float(x) for x in row] for row in W]
    exposure_matrix = [[float(x) for x in row] for row in H]

    if cached:
        cached.signature_matrix = signature_matrix
        cached.exposure_matrix = exposure_matrix
        cached.reconstruction_error = error
        cached.iterations = iterations
        cached.sample_names = spectrum_out.sample_names
        cached.mutation_types = spectrum_out.mutation_types
        cached.data_hash = data_hash
        cached.is_stale = False
        cached.last_checked_at = func.now()
        cached.created_at = func.now()
    else:
        new_cache = models.NMFSignatureResult(
            cache_key=cache_key,
            sample_ids=sorted_ids,
            reference_name=reference_name,
            k_value=k_value,
            signature_matrix=signature_matrix,
            exposure_matrix=exposure_matrix,
            mutation_types=spectrum_out.mutation_types,
            sample_names=spectrum_out.sample_names,
            reconstruction_error=error,
            iterations=iterations,
            data_hash=data_hash,
        )
        db.add(new_cache)
    db.commit()

    return _build_nmf_out(
        reference_name,
        k_value,
        sorted_ids,
        spectrum_out.sample_names,
        signature_matrix,
        exposure_matrix,
        spectrum_out.mutation_types,
        error,
        iterations,
        False,
    )


def _nmf_cache_to_out(
    cached: models.NMFSignatureResult,
    reference_name: str,
) -> schemas.NMFSignatureOut:
    return _build_nmf_out(
        reference_name,
        cached.k_value,
        cached.sample_ids,
        cached.sample_names,
        cached.signature_matrix,
        cached.exposure_matrix,
        cached.mutation_types,
        cached.reconstruction_error,
        cached.iterations,
        True,
        is_stale=cached.is_stale,
        created_at=cached.created_at,
    )


def _build_nmf_out(
    reference_name: str,
    k_value: int,
    sample_ids: List[int],
    sample_names: List[str],
    signature_matrix: List[List[float]],
    exposure_matrix: List[List[float]],
    mutation_types: List[str],
    reconstruction_error: float,
    iterations: int,
    cached: bool,
    is_stale: bool = False,
    created_at: Optional[datetime] = None,
) -> schemas.NMFSignatureOut:
    signatures = []
    for i in range(k_value):
        sig_probs = [signature_matrix[j][i] for j in range(96)]
        sig_probs = normalize_vector(sig_probs)
        signatures.append(schemas.SignatureOut(
            signature_index=i,
            mutation_types=mutation_types,
            probabilities=sig_probs,
        ))

    exposures = []
    for j in range(len(sample_ids)):
        exp_vals = [exposure_matrix[i][j] for i in range(k_value)]
        exposures.append(schemas.SampleExposureOut(
            sample_id=sample_ids[j],
            sample_name=sample_names[j],
            signature_exposures=exp_vals,
        ))

    return schemas.NMFSignatureOut(
        reference_name=reference_name,
        k_value=k_value,
        sample_ids=sample_ids,
        sample_names=sample_names,
        signatures=signatures,
        exposures=exposures,
        reconstruction_error=round(reconstruction_error, 6),
        iterations=iterations,
        cached=cached,
        is_stale=is_stale,
        created_at=created_at,
    )


def find_optimal_k_value(
    db: Session,
    sample_ids: List[int],
    reference_name: str,
) -> schemas.OptimalKOut:
    reference = _get_reference_or_error(db, reference_name)
    samples = _get_samples_or_error(db, sample_ids)
    sorted_ids = sorted(sample_ids)

    cache_key = compute_cache_key("optimal_k", sorted_ids, reference_name)
    data_hash = compute_data_hash(sorted_ids, reference_name)

    cached = db.query(models.OptimalKResult).filter(
        models.OptimalKResult.cache_key == cache_key
    ).first()

    if cached and not cached.is_stale and cached.data_hash == data_hash:
        cached.last_checked_at = func.now()
        db.commit()
        return _optimal_k_cache_to_out(cached, reference_name)

    spectrum_out = build_mutational_spectrum(db, sorted_ids, reference_name)
    V = _count_matrix_to_float(spectrum_out.count_matrix)

    recon_errors, cophenetic_corrs, recommended_k = find_optimal_k(
        V,
        k_min=2,
        k_max=8,
        n_runs_per_k=10,
        max_iter=5000,
        tol=1e-6,
        patience=20,
    )

    k_values = sorted(recon_errors.keys())
    reconstruction_errors = [round(recon_errors[k], 6) for k in k_values]
    cophenetic_correlations = [round(cophenetic_corrs[k], 6) for k in k_values]

    if cached:
        cached.k_values = k_values
        cached.reconstruction_errors = reconstruction_errors
        cached.cophenetic_correlations = cophenetic_correlations
        cached.recommended_k = recommended_k
        cached.data_hash = data_hash
        cached.is_stale = False
        cached.last_checked_at = func.now()
        cached.created_at = func.now()
    else:
        new_cache = models.OptimalKResult(
            cache_key=cache_key,
            sample_ids=sorted_ids,
            reference_name=reference_name,
            k_values=k_values,
            reconstruction_errors=reconstruction_errors,
            cophenetic_correlations=cophenetic_correlations,
            recommended_k=recommended_k,
            data_hash=data_hash,
        )
        db.add(new_cache)
    db.commit()

    metrics = [
        schemas.KMetricEntry(
            k_value=k,
            reconstruction_error=recon_errors[k],
            cophenetic_correlation=cophenetic_corrs[k],
        )
        for k in k_values
    ]

    return schemas.OptimalKOut(
        reference_name=reference_name,
        sample_ids=sorted_ids,
        sample_names=spectrum_out.sample_names,
        metrics=metrics,
        recommended_k=recommended_k,
        cached=False,
        is_stale=False,
        created_at=datetime.now(),
    )


def _optimal_k_cache_to_out(
    cached: models.OptimalKResult,
    reference_name: str,
) -> schemas.OptimalKOut:
    metrics = [
        schemas.KMetricEntry(
            k_value=k,
            reconstruction_error=err,
            cophenetic_correlation=corr,
        )
        for k, err, corr in zip(
            cached.k_values,
            cached.reconstruction_errors,
            cached.cophenetic_correlations,
        )
    ]

    return schemas.OptimalKOut(
        reference_name=reference_name,
        sample_ids=cached.sample_ids,
        sample_names=[],
        metrics=metrics,
        recommended_k=cached.recommended_k,
        cached=True,
        is_stale=cached.is_stale,
        created_at=cached.created_at,
    )


def match_signature_to_reference(
    db: Session,
    query_signature: List[float],
    top_n: int = 3,
    similarity_threshold: float = 0.8,
) -> schemas.SignatureComparisonOut:
    seed_reference_signatures(db, models.ReferenceSignature)

    normalized_query = normalize_vector(query_signature)

    ref_signatures = db.query(models.ReferenceSignature).all()

    matches = []
    for ref_sig in ref_signatures:
        ref_probs = ref_sig.probabilities
        if len(ref_probs) != 96:
            continue
        similarity = cosine_similarity(normalized_query, ref_probs)
        matches.append((ref_sig, similarity))

    matches.sort(key=lambda x: x[1], reverse=True)
    top_matches = matches[:top_n]

    result_matches = []
    for ref_sig, similarity in top_matches:
        result_matches.append(schemas.SignatureMatchOut(
            reference_signature=schemas.ReferenceSignatureOut(
                signature_id=ref_sig.signature_id,
                name=ref_sig.name,
                description=ref_sig.description,
                etiology=ref_sig.etiology,
            ),
            cosine_similarity=round(similarity, 6),
            is_reliable=similarity >= similarity_threshold,
        ))

    return schemas.SignatureComparisonOut(
        query_signature=normalized_query,
        top_matches=result_matches,
    )


def list_reference_signatures(db: Session) -> schemas.ReferenceSignatureListOut:
    seed_reference_signatures(db, models.ReferenceSignature)

    ref_signatures = db.query(models.ReferenceSignature).order_by(
        models.ReferenceSignature.signature_id
    ).all()

    signatures = [
        schemas.ReferenceSignatureOut(
            signature_id=sig.signature_id,
            name=sig.name,
            description=sig.description,
            etiology=sig.etiology,
        )
        for sig in ref_signatures
    ]

    return schemas.ReferenceSignatureListOut(
        total=len(signatures),
        signatures=signatures,
    )


def invalidate_signature_results_for_samples(db: Session, sample_ids: List[int]) -> None:
    sample_id_set = set(sample_ids)

    spectrum_caches = db.query(models.MutationalSpectrumCache).all()
    for cache in spectrum_caches:
        if any(sid in sample_id_set for sid in cache.sample_ids):
            cache.is_stale = True

    nmf_caches = db.query(models.NMFSignatureResult).all()
    for cache in nmf_caches:
        if any(sid in sample_id_set for sid in cache.sample_ids):
            cache.is_stale = True

    optimal_k_caches = db.query(models.OptimalKResult).all()
    for cache in optimal_k_caches:
        if any(sid in sample_id_set for sid in cache.sample_ids):
            cache.is_stale = True

    db.commit()
