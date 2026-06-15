from app.ld.ld_calculator import compute_ld_pair, compute_ld_matrix
from app.ld.haplotype_blocks import find_haplotype_blocks
from app.services.ld_service import _compute_snp_maf, _filter_snps


def test_ld_calculation():
    print("=== Testing LD Calculation ===")
    
    # Test case 1: Perfect LD (all variants co-occur)
    geno_a = [True, True, True, False, False, False, False, False, False, False]
    geno_b = [True, True, True, False, False, False, False, False, False, False]
    D, Dp, r2 = compute_ld_pair(geno_a, geno_b)
    print(f"Test 1 - Perfect LD: D={D:.4f}, D'={Dp:.4f}, r²={r2:.4f}")
    assert r2 > 0.99, "Perfect LD should have r² ≈ 1"
    assert abs(Dp) > 0.99, "Perfect LD should have D' ≈ 1"
    
    # Test case 2: No LD (independent)
    geno_a = [True, True, True, True, True, False, False, False, False, False]
    geno_b = [True, True, True, True, False, True, True, True, False, False]
    D, Dp, r2 = compute_ld_pair(geno_a, geno_b)
    print(f"Test 2 - Partial LD: D={D:.4f}, D'={Dp:.4f}, r²={r2:.4f}")
    
    # Test case 3: Complete equilibrium
    geno_a = [True, True, False, False]
    geno_b = [True, False, True, False]
    D, Dp, r2 = compute_ld_pair(geno_a, geno_b)
    print(f"Test 3 - Equilibrium: D={D:.4f}, D'={Dp:.4f}, r²={r2:.4f}")
    assert abs(D) < 0.001, "Equilibrium should have D ≈ 0"
    assert r2 < 0.01, "Equilibrium should have r² ≈ 0"
    
    # Test case 4: Negative D
    geno_a = [True, True, False, False, False, False]
    geno_b = [False, False, True, True, False, False]
    D, Dp, r2 = compute_ld_pair(geno_a, geno_b)
    print(f"Test 4 - Negative D: D={D:.4f}, D'={Dp:.4f}, r²={r2:.4f}")
    assert D < 0, "Should have negative D"
    
    # Test case 5: All reference
    geno_a = [False, False, False, False]
    geno_b = [False, False, False, False]
    D, Dp, r2 = compute_ld_pair(geno_a, geno_b)
    print(f"Test 5 - All reference: D={D:.4f}, D'={Dp:.4f}, r²={r2:.4f}")
    assert r2 == 0.0, "All reference should have r² = 0"
    
    print("\nAll LD calculation tests passed!")


def test_ld_matrix():
    print("\n=== Testing LD Matrix ===")
    
    sample_ids = [1, 2, 3, 4, 5, 6]
    snp_positions = [100, 200, 300]
    
    genotypes_by_sample = {
        1: {100: True, 200: True, 300: True},
        2: {100: True, 200: True, 300: True},
        3: {100: True, 200: True, 300: False},
        4: {100: False, 200: False, 300: False},
        5: {100: False, 200: False, 300: False},
        6: {100: False, 200: False, 300: True},
    }
    
    ld_pairs = compute_ld_matrix(sample_ids, snp_positions, genotypes_by_sample)
    
    print(f"Number of pairs: {len(ld_pairs)}")
    expected_pairs = len(snp_positions) * (len(snp_positions) - 1) // 2
    assert len(ld_pairs) == expected_pairs, f"Expected {expected_pairs} pairs"
    
    for pair in ld_pairs:
        print(f"  {pair['pos_i']} - {pair['pos_j']}: r²={pair['r_squared']:.4f}, D'={pair['d_prime']:.4f}")
        assert "d" in pair
        assert "d_prime" in pair
        assert "r_squared" in pair
        assert pair["pos_i"] < pair["pos_j"]
    
    print("\nLD matrix tests passed!")


def test_haplotype_blocks():
    print("\n=== Testing Haplotype Block Detection ===")
    
    snp_positions = [100, 200, 300, 400, 500, 600]
    
    ld_pairs = [
        {"pos_i": 100, "pos_j": 200, "r_squared": 0.95},
        {"pos_i": 100, "pos_j": 300, "r_squared": 0.90},
        {"pos_i": 100, "pos_j": 400, "r_squared": 0.30},
        {"pos_i": 100, "pos_j": 500, "r_squared": 0.20},
        {"pos_i": 100, "pos_j": 600, "r_squared": 0.10},
        {"pos_i": 200, "pos_j": 300, "r_squared": 0.92},
        {"pos_i": 200, "pos_j": 400, "r_squared": 0.25},
        {"pos_i": 200, "pos_j": 500, "r_squared": 0.15},
        {"pos_i": 200, "pos_j": 600, "r_squared": 0.05},
        {"pos_i": 300, "pos_j": 400, "r_squared": 0.20},
        {"pos_i": 300, "pos_j": 500, "r_squared": 0.10},
        {"pos_i": 300, "pos_j": 600, "r_squared": 0.05},
        {"pos_i": 400, "pos_j": 500, "r_squared": 0.97},
        {"pos_i": 400, "pos_j": 600, "r_squared": 0.93},
        {"pos_i": 500, "pos_j": 600, "r_squared": 0.98},
    ]
    
    blocks = find_haplotype_blocks(snp_positions, ld_pairs, r2_threshold=0.8)
    
    print(f"Found {len(blocks)} blocks")
    
    for block in blocks:
        print(f"  Block {block['block_index']}: {block['start_pos']}-{block['end_pos']} "
              f"({block['snp_count']} SNPs), avg r²={block['avg_r_squared']:.4f}")
        print(f"    SNPs: {block['snp_positions']}")
    
    assert len(blocks) == 2, f"Expected 2 blocks, got {len(blocks)}"
    
    assert blocks[0]["snp_count"] == 3
    assert blocks[0]["start_pos"] == 100
    assert blocks[0]["end_pos"] == 300
    
    assert blocks[1]["snp_count"] == 3
    assert blocks[1]["start_pos"] == 400
    assert blocks[1]["end_pos"] == 600
    
    print("\nHaplotype block tests passed!")


def test_haplotype_blocks_all_pairs_check():
    print("\n=== Testing All-Pairs Check in Haplotype Blocks ===")
    
    snp_positions = [100, 200, 300, 400]
    
    ld_pairs = [
        {"pos_i": 100, "pos_j": 200, "r_squared": 0.95},
        {"pos_i": 100, "pos_j": 300, "r_squared": 0.70},
        {"pos_i": 100, "pos_j": 400, "r_squared": 0.50},
        {"pos_i": 200, "pos_j": 300, "r_squared": 0.95},
        {"pos_i": 200, "pos_j": 400, "r_squared": 0.70},
        {"pos_i": 300, "pos_j": 400, "r_squared": 0.95},
    ]
    
    blocks = find_haplotype_blocks(snp_positions, ld_pairs, r2_threshold=0.8)
    
    print(f"Found {len(blocks)} blocks")
    for block in blocks:
        print(f"  Block {block['block_index']}: {block['start_pos']}-{block['end_pos']} "
              f"({block['snp_count']} SNPs)")
    
    assert len(blocks) == 2, f"Expected 2 blocks (100-200 and 300-400), got {len(blocks)}"
    
    assert blocks[0]["snp_count"] == 2
    assert blocks[0]["start_pos"] == 100
    assert blocks[0]["end_pos"] == 200
    
    assert blocks[1]["snp_count"] == 2
    assert blocks[1]["start_pos"] == 300
    assert blocks[1]["end_pos"] == 400
    
    print("\nAll-pairs check test passed!")


def test_maf_calculation():
    print("\n=== Testing MAF Calculation ===")
    
    sample_ids = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    
    genotypes_by_sample = {
        1: {100: True, 200: True, 300: True},
        2: {100: True, 200: True, 300: True},
        3: {100: True, 200: False, 300: False},
        4: {100: True, 200: False, 300: False},
        5: {100: False, 200: False, 300: False},
        6: {100: False, 200: False, 300: False},
        7: {100: False, 200: False, 300: False},
        8: {100: False, 200: False, 300: False},
        9: {100: False, 200: False, 300: False},
        10: {100: False, 200: False, 300: False},
    }
    
    maf_100 = _compute_snp_maf(100, sample_ids, genotypes_by_sample)
    print(f"SNP at 100 (4/10 variant): MAF = {maf_100:.4f} (expected 0.4000)")
    assert abs(maf_100 - 0.4) < 0.001, f"Expected MAF 0.4, got {maf_100}"
    
    maf_200 = _compute_snp_maf(200, sample_ids, genotypes_by_sample)
    print(f"SNP at 200 (2/10 variant): MAF = {maf_200:.4f} (expected 0.2000)")
    assert abs(maf_200 - 0.2) < 0.001, f"Expected MAF 0.2, got {maf_200}"
    
    maf_300 = _compute_snp_maf(300, sample_ids, genotypes_by_sample)
    print(f"SNP at 300 (2/10 variant): MAF = {maf_300:.4f} (expected 0.2000)")
    assert abs(maf_300 - 0.2) < 0.001, f"Expected MAF 0.2, got {maf_300}"
    
    maf_rare = _compute_snp_maf(400, sample_ids, genotypes_by_sample)
    print(f"SNP at 400 (0/10 variant): MAF = {maf_rare:.4f} (expected 0.0000)")
    assert abs(maf_rare - 0.0) < 0.001, f"Expected MAF 0.0, got {maf_rare}"
    
    sample_ids_2 = [1, 2, 3, 4, 5]
    genotypes_by_sample_2 = {
        1: {500: True},
        2: {500: True},
        3: {500: True},
        4: {500: True},
        5: {500: True},
    }
    maf_all_variant = _compute_snp_maf(500, sample_ids_2, genotypes_by_sample_2)
    print(f"SNP at 500 (5/5 variant): MAF = {maf_all_variant:.4f} (expected 0.0000)")
    assert abs(maf_all_variant - 0.0) < 0.001, f"Expected MAF 0.0 (all variant), got {maf_all_variant}"
    
    print("\nMAF calculation tests passed!")


def test_snp_filtering():
    print("\n=== Testing SNP Filtering ===")
    
    sample_ids = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    snp_positions = [100, 200, 300, 400, 500, 600]
    
    genotypes_by_sample = {}
    for sid in sample_ids:
        genotypes_by_sample[sid] = {}
    
    genotypes_by_sample[1][100] = True
    genotypes_by_sample[2][100] = True
    genotypes_by_sample[1][200] = True
    genotypes_by_sample[1][300] = True
    genotypes_by_sample[2][300] = True
    genotypes_by_sample[3][300] = True
    genotypes_by_sample[4][300] = True
    for sid in sample_ids[:5]:
        genotypes_by_sample[sid][400] = True
    for sid in sample_ids:
        genotypes_by_sample[sid][500] = True
    
    print("\nTest 1: No filtering")
    filtered, count = _filter_snps(sample_ids, snp_positions, genotypes_by_sample)
    print(f"  Original: {len(snp_positions)}, Filtered: {len(filtered)}, Removed: {count}")
    assert count == 0, f"Expected 0 removed, got {count}"
    assert len(filtered) == len(snp_positions), f"Expected {len(snp_positions)} SNPs, got {len(filtered)}"
    
    print("\nTest 2: Region filtering (start=250)")
    filtered, count = _filter_snps(sample_ids, snp_positions, genotypes_by_sample, start=250)
    print(f"  Original: {len(snp_positions)}, Filtered: {len(filtered)}, Removed: {count}")
    print(f"  Filtered positions: {filtered}")
    assert count == 2, f"Expected 2 removed, got {count}"
    assert 100 not in filtered
    assert 200 not in filtered
    assert 300 in filtered
    
    print("\nTest 3: Region filtering (end=450)")
    filtered, count = _filter_snps(sample_ids, snp_positions, genotypes_by_sample, end=450)
    print(f"  Original: {len(snp_positions)}, Filtered: {len(filtered)}, Removed: {count}")
    print(f"  Filtered positions: {filtered}")
    assert count == 2, f"Expected 2 removed, got {count}"
    assert 500 not in filtered
    assert 600 not in filtered
    assert 400 in filtered
    
    print("\nTest 4: Region filtering (start=250, end=450)")
    filtered, count = _filter_snps(sample_ids, snp_positions, genotypes_by_sample, start=250, end=450)
    print(f"  Original: {len(snp_positions)}, Filtered: {len(filtered)}, Removed: {count}")
    print(f"  Filtered positions: {filtered}")
    assert count == 4, f"Expected 4 removed, got {count}"
    assert filtered == [300, 400], f"Expected [300, 400], got {filtered}"
    
    print("\nTest 5: MAF filtering (min_maf=0.3)")
    filtered, count = _filter_snps(sample_ids, snp_positions, genotypes_by_sample, min_maf=0.3)
    print(f"  Original: {len(snp_positions)}, Filtered: {len(filtered)}, Removed: {count}")
    print(f"  Filtered positions: {filtered}")
    assert 100 not in filtered, "SNP 100 (MAF 0.2) should be filtered out"
    assert 200 not in filtered, "SNP 200 (MAF 0.1) should be filtered out"
    assert 300 in filtered, "SNP 300 (MAF 0.4) should be kept"
    assert 400 in filtered, "SNP 400 (MAF 0.5) should be kept"
    assert 500 not in filtered, "SNP 500 (MAF 0.0) should be filtered out"
    
    print("\nTest 6: Combined region and MAF filtering")
    filtered, count = _filter_snps(
        sample_ids, snp_positions, genotypes_by_sample,
        start=250, end=450, min_maf=0.3
    )
    print(f"  Original: {len(snp_positions)}, Filtered: {len(filtered)}, Removed: {count}")
    print(f"  Filtered positions: {filtered}")
    assert 300 in filtered, "SNP 300 should pass both filters"
    assert 400 in filtered, "SNP 400 should pass both filters"
    assert len(filtered) == 2, f"Expected 2 SNPs, got {len(filtered)}"
    
    print("\nSNP filtering tests passed!")


if __name__ == "__main__":
    test_ld_calculation()
    test_ld_matrix()
    test_haplotype_blocks()
    test_haplotype_blocks_all_pairs_check()
    test_maf_calculation()
    test_snp_filtering()
    print("\n" + "=" * 50)
    print("ALL TESTS PASSED!")
    print("=" * 50)
