import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from app.mutation_signature.matrix_ops import (
    matmul, transpose, cosine_similarity, normalize_vector,
    create_matrix, create_random_matrix, frobenius_norm, subtract,
    elementwise_multiply, elementwise_divide
)
from app.mutation_signature.trinucleotide import (
    ALL_96_MUTATION_TYPES, MUTATION_TYPE_INDEX,
    get_trinucleotide_context, format_mutation_type,
    get_mutation_type_index, compute_data_hash, compute_cache_key
)
from app.mutation_signature.nmf import (
    nmf_multiplicative_update, nmf_with_multiple_initializations,
    find_optimal_k, cophenetic_correlation, sample_connectivity_matrix
)
from app.mutation_signature.reference_signatures import (
    SIGNATURE_DESCRIPTIONS, get_all_reference_signatures
)
import random


def test_matrix_ops():
    print("=== Testing Matrix Operations ===")
    A = [[1.0, 2.0], [3.0, 4.0]]
    B = [[5.0, 6.0], [7.0, 8.0]]
    print(f"A = {A}")
    print(f"B = {B}")
    result = matmul(A, B)
    print(f"A * B = {result}")
    assert result == [[19.0, 22.0], [43.0, 50.0]], f"Matrix multiply failed: {result}"
    
    At = transpose(A)
    print(f"A^T = {At}")
    assert At == [[1.0, 3.0], [2.0, 4.0]], f"Transpose failed: {At}"
    
    sim1 = cosine_similarity([1.0, 0.0], [0.0, 1.0])
    print(f"cosine_similarity([1,0], [0,1]) = {sim1}")
    assert abs(sim1 - 0.0) < 1e-10, f"Cosine similarity orthogonal failed: {sim1}"
    
    sim2 = cosine_similarity([1.0, 0.0], [1.0, 0.0])
    print(f"cosine_similarity([1,0], [1,0]) = {sim2}")
    assert abs(sim2 - 1.0) < 1e-10, f"Cosine similarity identical failed: {sim2}"
    
    nv = normalize_vector([1.0, 2.0, 3.0])
    print(f"normalize_vector([1,2,3]) = {nv}")
    assert abs(sum(nv) - 1.0) < 1e-10, f"Normalize sum failed: {sum(nv)}"
    
    print("✓ Matrix operations passed\n")


def test_trinucleotide():
    print("=== Testing Trinucleotide Context ===")
    print(f"Total mutation types: {len(ALL_96_MUTATION_TYPES)}")
    assert len(ALL_96_MUTATION_TYPES) == 96, f"Expected 96, got {len(ALL_96_MUTATION_TYPES)}"
    print(f"First 5 types: {ALL_96_MUTATION_TYPES[:5]}")
    print(f"Last 5 types: {ALL_96_MUTATION_TYPES[-5:]}")
    
    ref_seq = "ATCGATCGATCG"
    left, ref, alt, right = get_trinucleotide_context(ref_seq, 3, "G", "A")
    mt = format_mutation_type(left, ref, alt, right)
    print(f"Position 3 in '{ref_seq}': {mt}")
    assert mt == "T[C>T]G", f"Expected T[C>T]G (complemented from G>A), got {mt}"
    
    idx = get_mutation_type_index("A[C>A]A")
    print(f"Index of 'A[C>A]A': {idx}")
    assert idx == 0, f"Expected 0, got {idx}"
    
    ref_seq2 = "AACGTT"
    left2, ref2, alt2, right2 = get_trinucleotide_context(ref_seq2, 2, "C", "T")
    mt2 = format_mutation_type(left2, ref2, alt2, right2)
    print(f"Position 2 in '{ref_seq2}', C>T: {mt2}")
    assert mt2 == "A[C>T]G", f"Expected A[C>T]G, got {mt2}"
    
    data_hash = compute_data_hash([1, 2, 3], "test_ref")
    print(f"Data hash: {data_hash[:16]}...")
    assert len(data_hash) == 64, f"SHA256 should be 64 chars"
    
    cache_key = compute_cache_key("test", [1, 2, 3], "test_ref", k_value=5)
    print(f"Cache key: {cache_key}")
    assert "k5" in cache_key
    
    print("✓ Trinucleotide operations passed\n")


def test_nmf():
    print("=== Testing NMF Algorithm ===")
    random.seed(42)
    
    m, n, k = 96, 10, 3
    V = [[random.random() for _ in range(n)] for _ in range(m)]
    print(f"Input V: {m}x{n} matrix")
    
    W, H, error, iterations = nmf_multiplicative_update(
        V, k, max_iter=200, tol=1e-4, patience=10, seed=42
    )
    
    print(f"W: {len(W)}x{len(W[0])} matrix")
    print(f"H: {len(H)}x{len(H[0])} matrix")
    print(f"Reconstruction error: {error:.6f}")
    print(f"Iterations: {iterations}")
    
    assert len(W) == m and len(W[0]) == k, f"W shape wrong: {len(W)}x{len(W[0])}"
    assert len(H) == k and len(H[0]) == n, f"H shape wrong: {len(H)}x{len(H[0])}"
    assert all(x >= 0 for row in W for x in row), "W has negative values"
    assert all(x >= 0 for row in H for x in row), "H has negative values"
    
    W_cols = [[W[i][j] for i in range(len(W))] for j in range(len(W[0]))]
    for j, col in enumerate(W_cols):
        s = sum(col)
        print(f"Column {j} sum: {s:.6f}")
        assert abs(s - 1.0) < 1e-6, f"Column {j} not normalized: {s}"
    
    print("\nTesting connectivity matrix...")
    cm = sample_connectivity_matrix(H)
    print(f"Connectivity matrix shape: {len(cm)}x{len(cm[0])}")
    assert len(cm) == n and len(cm[0]) == n
    
    print("\nTesting optimal K...")
    recon_errors, cophenetic_corrs, recommended_k = find_optimal_k(
        V, k_min=2, k_max=4, n_runs_per_k=3, max_iter=100, tol=1e-3, patience=5
    )
    print(f"Reconstruction errors: {recon_errors}")
    print(f"Cophenetic correlations: {cophenetic_corrs}")
    print(f"Recommended K: {recommended_k}")
    assert recommended_k >= 2 and recommended_k <= 4
    
    print("✓ NMF operations passed\n")


def test_reference_signatures():
    print("=== Testing Reference Signatures ===")
    print(f"Total signatures: {len(SIGNATURE_DESCRIPTIONS)}")
    assert len(SIGNATURE_DESCRIPTIONS) == 30, f"Expected 30, got {len(SIGNATURE_DESCRIPTIONS)}"
    
    sigs = get_all_reference_signatures()
    assert len(sigs) == 30
    
    first_sig = sigs[0]
    print(f"First signature: {first_sig['signature_id']} - {first_sig['name']}")
    print(f"  Etiology: {first_sig['etiology']}")
    print(f"  Probabilities length: {len(first_sig['probabilities'])}")
    print(f"  Probabilities sum: {sum(first_sig['probabilities']):.6f}")
    
    assert len(first_sig['probabilities']) == 96
    assert abs(sum(first_sig['probabilities']) - 1.0) < 1e-6
    
    for sig in sigs:
        assert len(sig['probabilities']) == 96
        assert abs(sum(sig['probabilities']) - 1.0) < 1e-6
    
    print("✓ Reference signatures passed\n")


def test_cophenetic():
    print("=== Testing Cophenetic Correlation ===")
    cm1 = [[1.0, 1.0, 0.0], [1.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
    cm2 = [[1.0, 1.0, 0.0], [1.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
    corr = cophenetic_correlation([cm1, cm2])
    print(f"Identical matrices correlation: {corr:.6f}")
    assert abs(corr - 1.0) < 1e-6
    
    print("✓ Cophenetic correlation passed\n")


if __name__ == "__main__":
    print("Running mutational signature module tests...\n")
    try:
        test_matrix_ops()
        test_trinucleotide()
        test_nmf()
        test_reference_signatures()
        test_cophenetic()
        print("=" * 50)
        print("ALL TESTS PASSED! ✓")
        print("=" * 50)
    except AssertionError as e:
        print(f"\nTEST FAILED: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
    except Exception as e:
        print(f"\nUNEXPECTED ERROR: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
