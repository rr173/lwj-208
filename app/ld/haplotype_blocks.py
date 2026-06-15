from typing import List, Dict, Tuple


def build_r2_matrix(
    snp_positions: List[int],
    ld_pairs: List[Dict]
) -> Dict[Tuple[int, int], float]:
    """
    Build a lookup dictionary for r² values between SNP pairs.
    
    Args:
        snp_positions: Sorted list of SNP positions
        ld_pairs: List of LD pair dicts with pos_i, pos_j, r_squared
    
    Returns:
        Dict mapping (pos_i, pos_j) -> r² where pos_i < pos_j
    """
    r2_map = {}
    for pair in ld_pairs:
        key = (pair["pos_i"], pair["pos_j"])
        r2_map[key] = pair["r_squared"]
    return r2_map


def _all_pairs_above_threshold(
    snp_indices: List[int],
    snp_positions: List[int],
    r2_map: Dict[Tuple[int, int], float],
    threshold: float
) -> bool:
    """
    Check if all pairs of SNPs in the given index list have r² >= threshold.
    """
    n = len(snp_indices)
    for i in range(n):
        for j in range(i + 1, n):
            pos_i = snp_positions[snp_indices[i]]
            pos_j = snp_positions[snp_indices[j]]
            key = (pos_i, pos_j)
            r2 = r2_map.get(key, 0.0)
            if r2 < threshold:
                return False
    return True


def _block_avg_r2(
    snp_indices: List[int],
    snp_positions: List[int],
    r2_map: Dict[Tuple[int, int], float]
) -> float:
    """
    Compute average r² across all pairs in a block.
    """
    n = len(snp_indices)
    if n < 2:
        return 1.0
    
    total = 0.0
    count = 0
    for i in range(n):
        for j in range(i + 1, n):
            pos_i = snp_positions[snp_indices[i]]
            pos_j = snp_positions[snp_indices[j]]
            key = (pos_i, pos_j)
            total += r2_map.get(key, 0.0)
            count += 1
    
    return total / count if count > 0 else 0.0


def find_haplotype_blocks(
    snp_positions: List[int],
    ld_pairs: List[Dict],
    r2_threshold: float = 0.8
) -> List[Dict]:
    """
    Identify haplotype blocks where all pairs within a block have r² >= threshold.
    
    Uses a greedy algorithm:
    - Start with each SNP as a potential block start
    - Extend the block to the right as long as all pairwise r² values remain >= threshold
    - Move to the next SNP after the end of the current block
    
    Args:
        snp_positions: Sorted list of SNP positions
        ld_pairs: List of LD pair dicts
        r2_threshold: Minimum r² threshold for block membership
    
    Returns:
        List of block dicts with keys: block_index, start_pos, end_pos,
        snp_count, snp_positions, avg_r_squared
    """
    n_snps = len(snp_positions)
    if n_snps < 2:
        return []

    r2_map = build_r2_matrix(snp_positions, ld_pairs)

    blocks = []
    block_index = 0
    i = 0

    while i < n_snps:
        block_end = i
        
        for j in range(i + 1, n_snps):
            candidate_indices = list(range(i, j + 1))
            if _all_pairs_above_threshold(candidate_indices, snp_positions, r2_map, r2_threshold):
                block_end = j
            else:
                break
        
        if block_end > i:
            block_indices = list(range(i, block_end + 1))
            block_snps = [snp_positions[idx] for idx in block_indices]
            avg_r2 = _block_avg_r2(block_indices, snp_positions, r2_map)
            
            blocks.append({
                "block_index": block_index,
                "start_pos": block_snps[0],
                "end_pos": block_snps[-1],
                "snp_count": len(block_snps),
                "snp_positions": block_snps,
                "avg_r_squared": avg_r2,
            })
            block_index += 1
            i = block_end + 1
        else:
            i += 1

    return blocks
