from typing import List, Dict, Tuple
from .trinucleotide import ALL_96_MUTATION_TYPES
from .matrix_ops import normalize_vector


SIGNATURE_DESCRIPTIONS = [
    {
        "id": "SBS1",
        "name": "Signature 1",
        "description": "Age-related signature due to spontaneous deamination of 5-methylcytosine",
        "etiology": "Endogenous - spontaneous deamination of 5-methylcytosine",
        "peak_types": ["N[C>T]G", "N[C>T]T"],
    },
    {
        "id": "SBS2",
        "name": "Signature 2",
        "description": "APOBEC cytidine deaminase activity signature",
        "etiology": "APOBEC family cytidine deaminases",
        "peak_types": ["T[C>T]T", "T[C>G]T"],
    },
    {
        "id": "SBS3",
        "name": "Signature 3",
        "description": "Defective DNA double-strand break repair by homologous recombination",
        "etiology": "BRCA1/BRCA2 deficiency, HRD",
        "peak_types": ["C[C>A]C", "T[T>C]T"],
    },
    {
        "id": "SBS4",
        "name": "Signature 4",
        "description": "Tobacco smoking-related signature",
        "etiology": "Tobacco smoke mutagens (polycyclic aromatic hydrocarbons)",
        "peak_types": ["C[A>G]C", "G[T>C]T"],
    },
    {
        "id": "SBS5",
        "name": "Signature 5",
        "description": "Clock-like signature of unknown etiology",
        "etiology": "Unknown - possibly endogenous processes",
        "peak_types": ["T[C>T]C", "C[T>C]T"],
    },
    {
        "id": "SBS6",
        "name": "Signature 6",
        "description": "Defective DNA mismatch repair",
        "etiology": "MMR deficiency (MSI)",
        "peak_types": ["C[T>C]T", "G[A>G]A"],
    },
    {
        "id": "SBS7a",
        "name": "Signature 7a",
        "description": "Ultraviolet light exposure signature",
        "etiology": "UV radiation (sunlight)",
        "peak_types": ["T[C>T]T", "C[C>T]T"],
    },
    {
        "id": "SBS7b",
        "name": "Signature 7b",
        "description": "Ultraviolet light exposure signature subtype",
        "etiology": "UV radiation (sunlight)",
        "peak_types": ["T[C>T]C", "T[C>T]T"],
    },
    {
        "id": "SBS7c",
        "name": "Signature 7c",
        "description": "Ultraviolet light exposure signature subtype",
        "etiology": "UV radiation (sunlight)",
        "peak_types": ["G[C>T]A", "T[C>T]G"],
    },
    {
        "id": "SBS7d",
        "name": "Signature 7d",
        "description": "Ultraviolet light exposure signature subtype",
        "etiology": "UV radiation (sunlight)",
        "peak_types": ["T[C>T]A", "A[C>T]T"],
    },
    {
        "id": "SBS8",
        "name": "Signature 8",
        "description": "Signature of unknown etiology common in many cancers",
        "etiology": "Unknown",
        "peak_types": ["T[A>C]A", "A[T>C]T"],
    },
    {
        "id": "SBS9",
        "name": "Signature 9",
        "description": "Polymerase eta somatic hypermutation",
        "etiology": "Pol eta activity during somatic hypermutation",
        "peak_types": ["C[T>A]T", "G[A>T]A"],
    },
    {
        "id": "SBS10a",
        "name": "Signature 10a",
        "description": "POLE exonuclease domain mutation signature",
        "etiology": "POLE proofreading domain mutation",
        "peak_types": ["T[C>T]T", "C[T>G]T"],
    },
    {
        "id": "SBS10b",
        "name": "Signature 10b",
        "description": "POLE exonuclease domain mutation signature subtype",
        "etiology": "POLE proofreading domain mutation",
        "peak_types": ["G[T>G]T", "A[C>G]C"],
    },
    {
        "id": "SBS11",
        "name": "Signature 11",
        "description": "Temozolomide treatment signature",
        "etiology": "Alkylating agent (temozolomide) treatment",
        "peak_types": ["N[C>T]N"],
    },
    {
        "id": "SBS12",
        "name": "Signature 12",
        "description": "Signature of unknown etiology in liver cancer",
        "etiology": "Unknown - possibly aristolochic acid",
        "peak_types": ["C[T>A]G", "G[A>T]C"],
    },
    {
        "id": "SBS13",
        "name": "Signature 13",
        "description": "APOBEC cytidine deaminase activity signature",
        "etiology": "APOBEC family cytidine deaminases",
        "peak_types": ["T[C>G]T", "T[C>G]C"],
    },
    {
        "id": "SBS14",
        "name": "Signature 14",
        "description": "Defective NER and transcription-coupled NER",
        "etiology": "NER deficiency (e.g. ERCC2 mutations)",
        "peak_types": ["G[T>A]T", "A[C>A]T"],
    },
    {
        "id": "SBS15",
        "name": "Signature 15",
        "description": "Defective DNA mismatch repair with C>T at CpG",
        "etiology": "MMR deficiency with CpG methylation",
        "peak_types": ["C[C>T]G", "G[G>A]C"],
    },
    {
        "id": "SBS16",
        "name": "Signature 16",
        "description": "Signature of unknown etiology in liver cancer",
        "etiology": "Unknown - possibly alcohol metabolism",
        "peak_types": ["T[C>G]A", "A[G>C]T"],
    },
    {
        "id": "SBS17a",
        "name": "Signature 17a",
        "description": "Signature of unknown etiology in esophagus and stomach cancers",
        "etiology": "Unknown - possibly oxidative damage",
        "peak_types": ["T[T>C]A", "A[A>G]T"],
    },
    {
        "id": "SBS17b",
        "name": "Signature 17b",
        "description": "Signature of unknown etiology in esophagus and stomach cancers",
        "etiology": "Unknown - possibly oxidative damage",
        "peak_types": ["T[T>G]A", "A[A>C]T"],
    },
    {
        "id": "SBS18",
        "name": "Signature 18",
        "description": "Reactive oxygen species damage signature",
        "etiology": "Oxidative damage (ROS)",
        "peak_types": ["C[A>C]A", "T[G>T]T"],
    },
    {
        "id": "SBS19",
        "name": "Signature 19",
        "description": "Pilot study artifact signature",
        "etiology": "Sequencing artifact",
        "peak_types": ["C[G>T]C", "G[C>A]G"],
    },
    {
        "id": "SBS20",
        "name": "Signature 20",
        "description": "Concurrent POLD1 and MMR deficiency",
        "etiology": "POLD1 proofreading + MMR deficiency",
        "peak_types": ["T[T>G]T", "A[A>C]A"],
    },
    {
        "id": "SBS21",
        "name": "Signature 21",
        "description": "Defective DNA mismatch repair signature subtype",
        "etiology": "MMR deficiency (MSI) subtype",
        "peak_types": ["G[T>C]G", "C[A>G]C"],
    },
    {
        "id": "SBS22",
        "name": "Signature 22",
        "description": "Aristolochic acid exposure signature",
        "etiology": "Aristolochic acid exposure",
        "peak_types": ["A[T>A]T", "T[A>T]A"],
    },
    {
        "id": "SBS23",
        "name": "Signature 23",
        "description": "Signature of unknown etiology",
        "etiology": "Unknown",
        "peak_types": ["T[C>A]T", "A[G>T]A"],
    },
    {
        "id": "SBS24",
        "name": "Signature 24",
        "description": "Aflatoxin exposure signature",
        "etiology": "Aflatoxin B1 exposure",
        "peak_types": ["C[A>T]C", "G[T>A]G"],
    },
    {
        "id": "SBS30",
        "name": "Signature 30",
        "description": "Defective base excision repair (NTHL1)",
        "etiology": "NTHL1 deficiency - BER pathway defect",
        "peak_types": ["C[G>T]T", "A[C>A]G"],
    },
]


def generate_signature_probabilities(signature_idx: int) -> List[float]:
    import random
    rng = random.Random(signature_idx * 12345)

    probs = []
    for i, mt in enumerate(ALL_96_MUTATION_TYPES):
        desc = SIGNATURE_DESCRIPTIONS[signature_idx % len(SIGNATURE_DESCRIPTIONS)]
        peak_types = desc.get("peak_types", [])

        is_peak = False
        for peak_pattern in peak_types:
            left, middle, right = peak_pattern[0], peak_pattern[1:6], peak_pattern[6]
            if left == 'N' or mt[0] == left:
                if middle in mt:
                    if right == 'N' or mt[-1] == right:
                        is_peak = True
                        break

        if is_peak:
            prob = rng.uniform(0.02, 0.05)
        else:
            prob = rng.uniform(0.001, 0.01)
        probs.append(prob)

    return normalize_vector(probs)


def get_all_reference_signatures() -> List[Dict]:
    signatures = []
    for i, desc in enumerate(SIGNATURE_DESCRIPTIONS):
        probs = generate_signature_probabilities(i)
        signatures.append({
            "signature_id": desc["id"],
            "name": desc["name"],
            "description": desc["description"],
            "etiology": desc["etiology"],
            "mutation_types": ALL_96_MUTATION_TYPES,
            "probabilities": probs,
        })
    return signatures


def seed_reference_signatures(db, model_class) -> None:
    from sqlalchemy import func
    existing_count = db.query(func.count(model_class.id)).scalar()
    if existing_count > 0:
        return

    signatures = get_all_reference_signatures()
    for sig_data in signatures:
        sig = model_class(
            signature_id=sig_data["signature_id"],
            name=sig_data["name"],
            description=sig_data["description"],
            etiology=sig_data["etiology"],
            mutation_types=sig_data["mutation_types"],
            probabilities=sig_data["probabilities"],
        )
        db.add(sig)
    db.commit()
