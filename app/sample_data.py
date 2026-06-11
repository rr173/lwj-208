import os
from app.database import SessionLocal
from app.services import genome_service
from app import models


def build_sample_reference() -> str:
    """
    Build a 3000bp sample reference sequence with specific regions.
    Coordinates (0-based):
      - 0-99: intergenic (100bp)
      - 100-299: 5' UTR (200bp)
      - 300-799: Exon 1 (500bp)
      - 800-1499: Intron (700bp)
      - 1500-1999: Exon 2 (500bp)
      - 2000-2999: 3' intergenic (1000bp)
    """
    parts = []

    intergenic1 = "A" * 50 + "G" * 50
    parts.append(intergenic1)

    utr = "T" * 100 + "C" * 100
    parts.append(utr)

    exon1_unique_start = "ATGGCCATTGACGTAGCTAGCATCGATCGATCGATCGATCG"
    exon1_mid = "AGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGC"
    exon1_unique_end = "TACGATCGATCGATCGATCGATCGATCGATCGATCGATCGAT"
    exon1 = exon1_unique_start + exon1_mid * 5 + exon1_unique_end
    exon1 = exon1[:500]
    parts.append(exon1)

    intron = "GT" + "ATATATATAT" * 30 + "AG"
    intron = intron[:700]
    if len(intron) < 700:
        intron += "A" * (700 - len(intron))
    parts.append(intron)

    exon2_unique_start = "ATGCGTACGATCGATCGATCGATCGATCGATCGATCGATCGA"
    exon2_mid = "GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCT"
    exon2_unique_end = "CTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCT"
    exon2 = exon2_unique_start + exon2_mid * 5 + exon2_unique_end
    exon2 = exon2[:500]
    parts.append(exon2)

    intergenic2 = "G" * 500 + "C" * 500
    intergenic2 = intergenic2[:1000]
    parts.append(intergenic2)

    full = "".join(parts)
    if len(full) < 3000:
        full += "N" * (3000 - len(full))
    return full[:3000].upper()


SAMPLE_REFERENCE_NAME = "sample_chr1"

SAMPLE_ANNOTATIONS = [
    {"seq_name": SAMPLE_REFERENCE_NAME, "start": 100, "end": 299, "feature_type": "UTR", "gene_name": "SAMPLE_GENE"},
    {"seq_name": SAMPLE_REFERENCE_NAME, "start": 300, "end": 799, "feature_type": "exon", "gene_name": "SAMPLE_GENE"},
    {"seq_name": SAMPLE_REFERENCE_NAME, "start": 800, "end": 1499, "feature_type": "intron", "gene_name": "SAMPLE_GENE"},
    {"seq_name": SAMPLE_REFERENCE_NAME, "start": 1500, "end": 1999, "feature_type": "exon", "gene_name": "SAMPLE_GENE"},
]


def get_sample_queries(reference_seq: str) -> list:
    """
    Generate 5 sample query sequences with known variants.
    Ensures: at least 2 SNPs, 1 insertion, 1 deletion.
    """
    queries = []

    q1_seq = reference_seq[320:470]
    queries.append({
        "name": "query_perfect_match",
        "sequence": q1_seq,
        "description": "Perfect match in exon 1"
    })

    q2_list = list(reference_seq[350:500])
    snp_pos = 30
    original = q2_list[snp_pos]
    new_base = "G" if original in "AT" else "A"
    q2_list[snp_pos] = new_base
    queries.append({
        "name": "query_snp_missense",
        "sequence": "".join(q2_list),
        "description": "SNP (missense) in exon 1"
    })

    q3_list = list(reference_seq[400:550])
    snp_pos2 = 45
    original2 = q3_list[snp_pos2]
    ref_pos = 400 + snp_pos2
    codon_offset = (ref_pos - 300) % 3
    if codon_offset == 2 and original2 == "C":
        new_base2 = "T"
    elif original2 == "T":
        new_base2 = "C"
    elif original2 == "A":
        new_base2 = "G"
    else:
        new_base2 = "T"
    q3_list[snp_pos2] = new_base2
    queries.append({
        "name": "query_snp_synonymous",
        "sequence": "".join(q3_list),
        "description": "SNP (synonymous) in exon 1"
    })

    q4_list = list(reference_seq[380:530])
    ins_pos = 50
    q4_list[ins_pos:ins_pos] = list("GAT")
    queries.append({
        "name": "query_insertion",
        "sequence": "".join(q4_list),
        "description": "3bp insertion in exon 1"
    })

    q5_list = list(reference_seq[1550:1700])
    del_pos = 40
    del_len = 3
    del q5_list[del_pos:del_pos + del_len]
    queries.append({
        "name": "query_deletion",
        "sequence": "".join(q5_list),
        "description": "3bp deletion in exon 2 (inframe)"
    })

    return queries


def init_sample_data():
    """Initialize sample data on startup if enabled."""
    db = SessionLocal()
    try:
        existing = genome_service.get_reference_by_name(db, SAMPLE_REFERENCE_NAME)
        if existing:
            return

        ref_seq = build_sample_reference()
        ref = genome_service.create_reference_sequence(db, SAMPLE_REFERENCE_NAME, ref_seq)

        genome_service.add_gene_annotations(db, ref.id, SAMPLE_ANNOTATIONS)

        queries = get_sample_queries(ref_seq)
        for q in queries:
            genome_service.align_query_all_references(db, q["sequence"])

        db.commit()
    finally:
        db.close()
