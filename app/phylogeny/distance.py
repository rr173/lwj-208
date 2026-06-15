import math
from typing import List, Tuple, Dict, Set, Optional


IUPAC_AMBIGUITY = {
    'A': {'A'}, 'T': {'T'}, 'C': {'C'}, 'G': {'G'},
    'R': {'A', 'G'}, 'Y': {'C', 'T'}, 'S': {'G', 'C'}, 'W': {'A', 'T'},
    'K': {'G', 'T'}, 'M': {'A', 'C'}, 'B': {'C', 'G', 'T'},
    'D': {'A', 'G', 'T'}, 'H': {'A', 'C', 'T'}, 'V': {'A', 'C', 'G'},
    'N': {'A', 'T', 'C', 'G'}, '-': set(),
}

TRANSITIONS = {('A', 'G'), ('G', 'A'), ('C', 'T'), ('T', 'C')}


def _is_transition(base1: str, base2: str) -> bool:
    return (base1, base2) in TRANSITIONS


def _is_transversion(base1: str, base2: str) -> bool:
    b1 = base1.upper()
    b2 = base2.upper()
    if b1 == b2:
        return False
    if _is_transition(b1, b2):
        return False
    s1 = IUPAC_AMBIGUITY.get(b1, set())
    s2 = IUPAC_AMBIGUITY.get(b2, set())
    if not s1 or not s2:
        return False
    for a in s1:
        for b in s2:
            if a != b and not _is_transition(a, b):
                return True
    return False


def _count_sites(seq1: List[str], seq2: List[str]) -> Tuple[int, int, int]:
    """
    Count (same_sites, transitions, transversions) between two sequences.
    Only considers SNP variants (single-base).
    """
    same = 0
    transitions = 0
    transversions = 0
    
    for b1, b2 in zip(seq1, seq2):
        b1u = b1.upper()
        b2u = b2.upper()
        
        s1 = IUPAC_AMBIGUITY.get(b1u, set())
        s2 = IUPAC_AMBIGUITY.get(b2u, set())
        
        if not s1 or not s2:
            continue
        
        if s1 == s2 and len(s1) == 1 and len(s2) == 1:
            base1 = next(iter(s1))
            base2 = next(iter(s2))
            if base1 == base2:
                same += 1
            elif _is_transition(base1, base2):
                transitions += 1
            else:
                transversions += 1
        else:
            is_same = False
            has_transition = False
            has_transversion = False
            for a in s1:
                for b in s2:
                    if a == b:
                        is_same = True
                    elif _is_transition(a, b):
                        has_transition = True
                    else:
                        has_transversion = True
            
            if is_same and len(s1) == 1 and len(s2) == 1:
                same += 1
            elif has_transition and not has_transversion:
                transitions += 1
            elif has_transversion and not has_transition:
                transversions += 1
            elif has_transition and has_transversion:
                transversions += 1
    
    return same, transitions, transversions


def kimura_two_parameter_distance(
    seq1: List[str], 
    seq2: List[str]
) -> Tuple[Optional[float], Optional[str]]:
    """
    Calculate Kimura 2-parameter distance between two sequences.
    
    K2P distance formula:
        d = -0.5 * ln[(1 - 2p - q) * sqrt(1 - 2q)]
    
    where:
        p = proportion of transitional differences
        q = proportion of transversional differences
    
    Returns (distance, warning_message). If there are no common sites,
    returns (NaN, warning).
    """
    same, ti, tv = _count_sites(seq1, seq2)
    total = same + ti + tv
    
    if total == 0:
        return float('nan'), "No common variant sites covered by both samples"
    
    p = ti / total
    q = tv / total
    
    term1 = 1 - 2 * p - q
    term2 = 1 - 2 * q
    
    if term1 <= 0 or term2 <= 0:
        if q > 0:
            correction = (p + q) * (1 + 0.5 * (p + q)) + (p * p) / (2 * q * q)
        else:
            correction = -math.log(1 - p - q) if (p + q) < 1 else float('inf')
        return correction, None
    
    try:
        d = -0.5 * math.log(term1 * math.sqrt(term2))
        return d, None
    except (ValueError, ZeroDivisionError):
        return float('nan'), "Distance calculation failed due to numerical issues"


def build_sequence_matrix(
    variants_by_sample: Dict[int, Dict[int, str]],
    all_positions: List[int],
    sample_ids: List[int],
    ref_base_by_pos: Dict[int, str]
) -> Dict[int, List[str]]:
    """
    Build the multi-sequence alignment matrix for variant sites.
    
    Args:
        variants_by_sample: {sample_id: {pos: alt_base}}
        all_positions: sorted list of all variant positions
        sample_ids: list of sample IDs in order
        ref_base_by_pos: {pos: reference_base}
    
    Returns:
        {sample_id: [base_at_pos_1, base_at_pos_2, ...]}
    """
    matrix = {}
    for sid in sample_ids:
        seq = []
        sample_vars = variants_by_sample.get(sid, {})
        for pos in all_positions:
            if pos in sample_vars:
                seq.append(sample_vars[pos])
            elif pos in ref_base_by_pos:
                seq.append(ref_base_by_pos[pos])
            else:
                seq.append('N')
        matrix[sid] = seq
    return matrix


def compute_distance_matrix(
    sequence_matrix: Dict[int, List[str]],
    sample_ids: List[int]
) -> Tuple[List[List[float]], List[Dict], List[str]]:
    """
    Compute pairwise distance matrix using Kimura 2-parameter model.
    
    Args:
        sequence_matrix: {sample_id: [base1, base2, ...]}
        sample_ids: ordered list of sample IDs
    
    Returns:
        (distance_matrix, warnings, sample_names_in_order)
    """
    n = len(sample_ids)
    dist_matrix = [[0.0] * n for _ in range(n)]
    warnings = []
    
    for i in range(n):
        for j in range(i + 1, n):
            seq_i = sequence_matrix[sample_ids[i]]
            seq_j = sequence_matrix[sample_ids[j]]
            d, warn = kimura_two_parameter_distance(seq_i, seq_j)
            dist_matrix[i][j] = d
            dist_matrix[j][i] = d
            if warn:
                warnings.append({
                    "sample_i_id": sample_ids[i],
                    "sample_j_id": sample_ids[j],
                    "row": i,
                    "col": j,
                    "reason": warn
                })
    
    sample_names = [f"sample_{sid}" for sid in sample_ids]
    return dist_matrix, warnings, sample_names
