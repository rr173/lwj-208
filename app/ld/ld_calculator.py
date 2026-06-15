from typing import List, Dict, Tuple


def compute_ld_pair(
    genotypes_a: List[bool],
    genotypes_b: List[bool]
) -> Tuple[float, float, float]:
    """
    Compute linkage disequilibrium statistics for a pair of SNP loci.
    
    Args:
        genotypes_a: List of boolean values indicating variant status at locus A
                     (True = variant allele, False = reference allele)
        genotypes_b: List of boolean values indicating variant status at locus B
    
    Returns:
        Tuple of (D, D', r²)
    
    Note:
        Assumes haploid data (one allele per sample). For diploid data, this would
        need to be adjusted to account for phase information.
    """
    n = len(genotypes_a)
    if n == 0:
        return 0.0, 0.0, 0.0

    count_rr = 0
    count_ra = 0
    count_ar = 0
    count_aa = 0

    for ga, gb in zip(genotypes_a, genotypes_b):
        if not ga and not gb:
            count_rr += 1
        elif not ga and gb:
            count_ra += 1
        elif ga and not gb:
            count_ar += 1
        else:
            count_aa += 1

    p_a = (count_ar + count_aa) / n
    p_b = (count_ra + count_aa) / n
    p_ab = count_aa / n

    D = p_ab - p_a * p_b

    if D == 0:
        D_prime = 0.0
    elif D > 0:
        D_max = min(p_a * (1 - p_b), (1 - p_a) * p_b)
        D_prime = D / D_max if D_max != 0 else 0.0
    else:
        D_max = min(p_a * p_b, (1 - p_a) * (1 - p_b))
        D_prime = D / D_max if D_max != 0 else 0.0

    D_prime = abs(D_prime)

    denom = p_a * (1 - p_a) * p_b * (1 - p_b)
    if denom == 0:
        r_squared = 0.0
    else:
        r_squared = (D * D) / denom

    return D, D_prime, r_squared


def compute_ld_matrix(
    sample_ids: List[int],
    snp_positions: List[int],
    genotypes_by_sample: Dict[int, Dict[int, bool]]
) -> List[Dict]:
    """
    Compute pairwise LD matrix for all SNP positions.
    
    Args:
        sample_ids: List of sample IDs
        snp_positions: Sorted list of SNP positions
        genotypes_by_sample: Dict mapping sample_id -> {pos: is_variant}
    
    Returns:
        List of dicts with keys: pos_i, pos_j, d, d_prime, r_squared
        Only includes pairs where i < j (upper triangle)
    """
    n_snps = len(snp_positions)
    n_samples = len(sample_ids)

    genotype_matrix: Dict[int, List[bool]] = {}
    for pos in snp_positions:
        genotypes = []
        for sid in sample_ids:
            sample_geno = genotypes_by_sample.get(sid, {})
            genotypes.append(sample_geno.get(pos, False))
        genotype_matrix[pos] = genotypes

    ld_pairs = []
    for i in range(n_snps):
        for j in range(i + 1, n_snps):
            pos_i = snp_positions[i]
            pos_j = snp_positions[j]
            D, D_prime, r_squared = compute_ld_pair(
                genotype_matrix[pos_i],
                genotype_matrix[pos_j]
            )
            ld_pairs.append({
                "pos_i": pos_i,
                "pos_j": pos_j,
                "d": D,
                "d_prime": D_prime,
                "r_squared": r_squared,
            })

    return ld_pairs
