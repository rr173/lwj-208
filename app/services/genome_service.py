from typing import List, Optional, Tuple
from sqlalchemy.orm import Session
from app import models, schemas
from app.alignment.smith_waterman import optimized_smith_waterman, sequence_hash
from app.alignment.variant_caller import extract_variants
from app.annotation.gene_annotation import annotate_variants


def parse_fasta(fasta_text: str) -> List[Tuple[str, str]]:
    """Parse FASTA format text into list of (name, sequence) tuples."""
    sequences = []
    current_name = None
    current_seq_parts = []

    for line in fasta_text.strip().split("\n"):
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if current_name is not None:
                sequences.append((current_name, "".join(current_seq_parts)))
            current_name = line[1:].strip().split()[0]
            current_seq_parts = []
        else:
            current_seq_parts.append(line.upper())

    if current_name is not None:
        sequences.append((current_name, "".join(current_seq_parts)))

    return sequences


def parse_gff(gff_text: str) -> List[dict]:
    """Parse simplified GFF format into list of annotation dicts."""
    annotations = []
    for line in gff_text.strip().split("\n"):
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) >= 5:
            annotations.append({
                "seq_name": parts[0],
                "start": int(parts[1]),
                "end": int(parts[2]),
                "feature_type": parts[3],
                "gene_name": parts[4],
            })
    return annotations


def create_reference_sequence(db: Session, name: str, sequence: str) -> models.ReferenceSequence:
    """Create or update a reference sequence."""
    existing = db.query(models.ReferenceSequence).filter(
        models.ReferenceSequence.name == name
    ).first()

    if existing:
        existing.sequence = sequence.upper()
        existing.length = len(sequence)
        db.commit()
        db.refresh(existing)
        return existing

    ref = models.ReferenceSequence(
        name=name,
        sequence=sequence.upper(),
        length=len(sequence),
    )
    db.add(ref)
    db.commit()
    db.refresh(ref)
    return ref


def add_gene_annotations(db: Session, reference_id: int, annotations: List[dict]) -> List[models.GeneAnnotation]:
    """Add gene annotations to a reference sequence."""
    existing = db.query(models.GeneAnnotation).filter(
        models.GeneAnnotation.reference_id == reference_id
    ).all()
    for ann in existing:
        db.delete(ann)
    db.commit()

    created = []
    for ann in annotations:
        gff = models.GeneAnnotation(
            reference_id=reference_id,
            start=ann["start"],
            end=ann["end"],
            feature_type=ann["feature_type"],
            gene_name=ann["gene_name"],
        )
        db.add(gff)
        created.append(gff)
    db.commit()
    for g in created:
        db.refresh(g)
    return created


def get_reference_by_name(db: Session, name: str) -> Optional[models.ReferenceSequence]:
    """Get reference sequence by name."""
    return db.query(models.ReferenceSequence).filter(
        models.ReferenceSequence.name == name
    ).first()


def get_all_references(db: Session) -> List[models.ReferenceSequence]:
    """Get all reference sequences."""
    return db.query(models.ReferenceSequence).all()


def get_cached_alignment(db: Session, reference_id: int, query_hash: str) -> Optional[models.AlignmentResult]:
    """Get cached alignment result for a query on a reference."""
    return db.query(models.AlignmentResult).filter(
        models.AlignmentResult.reference_id == reference_id,
        models.AlignmentResult.query_hash == query_hash,
    ).first()


def perform_alignment(db: Session, reference: models.ReferenceSequence, query_sequence: str) -> models.AlignmentResult:
    """
    Perform alignment of query sequence against reference.
    Uses cache if available, otherwise computes and stores result.
    """
    query_hash = sequence_hash(query_sequence)
    cached = get_cached_alignment(db, reference.id, query_hash)
    if cached:
        return cached

    result = optimized_smith_waterman(query_sequence.upper(), reference.sequence)

    alignment = models.AlignmentResult(
        reference_id=reference.id,
        query_sequence=query_sequence.upper(),
        query_hash=query_hash,
        score=result["score"],
        ref_start=result["ref_start"],
        ref_end=result["ref_end"],
        query_start=result["query_start"],
        query_end=result["query_end"],
        cigar=result["cigar"],
        alignment_query=result["alignment_query"],
        alignment_match=result["alignment_match"],
        alignment_ref=result["alignment_ref"],
    )
    db.add(alignment)
    db.flush()

    variants = extract_variants(result, reference.sequence)

    annotations = db.query(models.GeneAnnotation).filter(
        models.GeneAnnotation.reference_id == reference.id
    ).all()
    ann_dicts = [
        {"start": a.start, "end": a.end, "feature_type": a.feature_type, "gene_name": a.gene_name, "strand": a.strand}
        for a in annotations
    ]

    annotated = annotate_variants(variants, ann_dicts, reference.sequence)

    for v in annotated:
        variant = models.Variant(
            alignment_id=alignment.id,
            reference_id=reference.id,
            variant_type=v["variant_type"],
            ref_pos=v["ref_pos"],
            ref_base=v["ref_base"],
            alt_base=v["alt_base"],
            feature_type=v.get("feature_type"),
            gene_name=v.get("gene_name"),
            impact=v.get("impact"),
            consequence=v.get("consequence"),
            codon_ref=v.get("codon_ref"),
            codon_alt=v.get("codon_alt"),
            aa_ref=v.get("aa_ref"),
            aa_alt=v.get("aa_alt"),
        )
        db.add(variant)

    db.commit()
    db.refresh(alignment)
    return alignment


def align_query_all_references(db: Session, query_sequence: str, threshold_ratio: float = 0.8) -> List[models.AlignmentResult]:
    """
    Align query against all reference sequences.
    Returns results with score > threshold_ratio * best_score, sorted by score descending.
    """
    references = get_all_references(db)
    if not references:
        return []

    all_results = []
    best_score = 0

    for ref in references:
        alignment = perform_alignment(db, ref, query_sequence)
        if alignment.score > best_score:
            best_score = alignment.score
        all_results.append(alignment)

    if best_score == 0:
        return []

    threshold = best_score * threshold_ratio
    filtered = [a for a in all_results if a.score >= threshold and a.score > 0]
    filtered.sort(key=lambda x: x.score, reverse=True)

    return filtered


def get_alignment_by_id(db: Session, alignment_id: int) -> Optional[models.AlignmentResult]:
    """Get alignment result by ID."""
    return db.query(models.AlignmentResult).filter(
        models.AlignmentResult.id == alignment_id
    ).first()


def get_variant_stats(db: Session, reference_id: int, bucket_size: int = 1000) -> dict:
    """
    Get variant statistics for a reference sequence.
    - Density by position (buckets of bucket_size)
    - Count by impact
    - Count by variant type
    """
    reference = db.query(models.ReferenceSequence).filter(
        models.ReferenceSequence.id == reference_id
    ).first()
    if not reference:
        return {"density_by_position": [], "count_by_impact": {}, "count_by_type": {}}

    variants = db.query(models.Variant).filter(
        models.Variant.reference_id == reference_id
    ).all()

    num_buckets = (reference.length + bucket_size - 1) // bucket_size
    density = []
    for i in range(num_buckets):
        start = i * bucket_size
        end = min((i + 1) * bucket_size - 1, reference.length - 1)
        count = sum(1 for v in variants if start <= v.ref_pos <= end)
        density.append({
            "bucket_start": start,
            "bucket_end": end,
            "count": count,
        })

    impact_counts = {}
    for v in variants:
        impact = v.impact or "UNKNOWN"
        impact_counts[impact] = impact_counts.get(impact, 0) + 1

    type_counts = {}
    for v in variants:
        vtype = v.variant_type
        type_counts[vtype] = type_counts.get(vtype, 0) + 1

    return {
        "density_by_position": density,
        "count_by_impact": impact_counts,
        "count_by_type": type_counts,
    }
