from typing import List, Dict, Tuple
import hashlib
import json

REF_BASES = ['A', 'C', 'G', 'T']
ALT_BASES = ['C', 'G', 'T', 'A']
CONTEXT_BASES = ['A', 'C', 'G', 'T']

SNP_SUBSTITUTIONS = [
    ('C', 'A'), ('C', 'G'), ('C', 'T'),
    ('T', 'A'), ('T', 'C'), ('T', 'G'),
]


def generate_all_mutation_types() -> List[str]:
    mutation_types = []
    for ref, alt in SNP_SUBSTITUTIONS:
        for left in CONTEXT_BASES:
            for right in CONTEXT_BASES:
                mutation_type = f"{left}[{ref}>{alt}]{right}"
                mutation_types.append(mutation_type)
    return mutation_types


ALL_96_MUTATION_TYPES = generate_all_mutation_types()
MUTATION_TYPE_INDEX = {mt: i for i, mt in enumerate(ALL_96_MUTATION_TYPES)}


def get_trinucleotide_context(
    reference_sequence: str,
    ref_pos: int,
    ref_base: str,
    alt_base: str,
) -> Tuple[str, str, str, str]:
    if ref_pos <= 0 or ref_pos >= len(reference_sequence) - 1:
        raise ValueError("Position at edge of reference sequence")

    left_base = reference_sequence[ref_pos - 1]
    right_base = reference_sequence[ref_pos + 1]

    actual_ref = reference_sequence[ref_pos]
    if actual_ref != ref_base:
        ref_base = actual_ref

    if ref_base in ['C', 'T']:
        return left_base, ref_base, alt_base, right_base
    else:
        comp_map = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
        left_base = comp_map[left_base]
        right_base = comp_map[right_base]
        ref_base = comp_map[ref_base]
        alt_base = comp_map[alt_base]
        return right_base, ref_base, alt_base, left_base


def format_mutation_type(left: str, ref: str, alt: str, right: str) -> str:
    return f"{left}[{ref}>{alt}]{right}"


def get_mutation_type_index(mutation_type: str) -> int:
    return MUTATION_TYPE_INDEX.get(mutation_type, -1)


def create_empty_count_matrix(n_samples: int) -> List[List[int]]:
    return [[0] * n_samples for _ in range(96)]


def compute_data_hash(
    sample_ids: List[int],
    reference_name: str,
    k_value: int = None,
) -> str:
    sorted_ids = sorted(sample_ids)
    data = {
        'sample_ids': sorted_ids,
        'reference_name': reference_name,
    }
    if k_value is not None:
        data['k_value'] = k_value
    json_str = json.dumps(data, sort_keys=True)
    return hashlib.sha256(json_str.encode()).hexdigest()


def compute_cache_key(
    prefix: str,
    sample_ids: List[int],
    reference_name: str,
    k_value: int = None,
) -> str:
    data_hash = compute_data_hash(sample_ids, reference_name, k_value)
    if k_value is not None:
        return f"{prefix}:{reference_name}:{data_hash}:k{k_value}"
    return f"{prefix}:{reference_name}:{data_hash}"
