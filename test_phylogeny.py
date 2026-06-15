import math
from app.phylogeny.distance import (
    kimura_two_parameter_distance,
    build_sequence_matrix,
    compute_distance_matrix,
)
from app.phylogeny.tree_building import (
    build_tree,
    upgma,
    neighbor_joining,
    get_tree_splits,
    robinson_foulds_distance,
    to_newick_format,
)
from app.phylogeny.molecular_clock import (
    linear_regression,
    compute_root_to_tip_distances,
    date_to_year,
)
from datetime import datetime


def test_kimura_distance():
    print("Testing Kimura 2-parameter distance...")
    
    seq1 = ['A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A']
    seq2 = ['A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A']
    d, warn = kimura_two_parameter_distance(seq1, seq2)
    assert d == 0.0, f"Identical sequences should have distance 0, got {d}"
    assert warn is None
    print("  ✓ Identical sequences: distance = 0")
    
    seq1 = ['A', 'A', 'A', 'A']
    seq2 = ['G', 'G', 'G', 'G']
    d, warn = kimura_two_parameter_distance(seq1, seq2)
    assert d > 0, f"All transitions should have positive distance, got {d}"
    print(f"  ✓ All transitions (A→G): distance = {d:.4f}")
    
    seq1 = ['A', 'A', 'A', 'A']
    seq2 = ['C', 'C', 'C', 'C']
    d, warn = kimura_two_parameter_distance(seq1, seq2)
    assert d > 0, f"All transversions should have positive distance, got {d}"
    print(f"  ✓ All transversions (A→C): distance = {d:.4f}")
    
    seq1 = ['N', 'N', 'N']
    seq2 = ['A', 'C', 'G']
    d, warn = kimura_two_parameter_distance(seq1, seq2)
    print(f"  ✓ Ambiguous bases handled: distance = {d:.4f}")
    
    seq1 = []
    seq2 = []
    d, warn = kimura_two_parameter_distance(seq1, seq2)
    assert math.isnan(d), "No common sites should return NaN"
    assert warn is not None
    print("  ✓ No common sites: returns NaN with warning")
    
    print("✓ Kimura distance tests passed\n")


def test_distance_matrix():
    print("Testing distance matrix computation...")
    
    variants_by_sample = {
        1: {100: 'A', 200: 'G', 300: 'C'},
        2: {100: 'A', 200: 'A', 300: 'C'},
        3: {100: 'T', 200: 'G', 300: 'G'},
    }
    all_positions = [100, 200, 300]
    sample_ids = [1, 2, 3]
    ref_base_by_pos = {100: 'A', 200: 'G', 300: 'C'}
    
    seq_matrix = build_sequence_matrix(
        variants_by_sample, all_positions, sample_ids, ref_base_by_pos
    )
    assert len(seq_matrix) == 3
    assert len(seq_matrix[1]) == 3
    print("  ✓ Sequence matrix built correctly")
    
    dist_matrix, warnings, names = compute_distance_matrix(seq_matrix, sample_ids)
    assert len(dist_matrix) == 3
    assert len(dist_matrix[0]) == 3
    assert dist_matrix[0][0] == 0.0
    assert dist_matrix[1][1] == 0.0
    assert dist_matrix[2][2] == 0.0
    assert dist_matrix[0][1] == dist_matrix[1][0]
    assert dist_matrix[0][2] == dist_matrix[2][0]
    assert dist_matrix[1][2] == dist_matrix[2][1]
    print(f"  ✓ Distance matrix is symmetric, diagonal = 0")
    print(f"  ✓ Distances: [{dist_matrix[0][1]:.4f}, {dist_matrix[0][2]:.4f}, {dist_matrix[1][2]:.4f}]")
    
    print("✓ Distance matrix tests passed\n")


def test_upgma():
    print("Testing UPGMA algorithm...")
    
    dist_matrix = [
        [0.0, 0.1, 0.4, 0.5],
        [0.1, 0.0, 0.4, 0.5],
        [0.4, 0.4, 0.0, 0.3],
        [0.5, 0.5, 0.3, 0.0],
    ]
    labels = ['A', 'B', 'C', 'D']
    
    root, newick = build_tree(dist_matrix, labels, method="upgma")
    
    assert root is not None
    assert len(root.children) == 2
    print(f"  ✓ UPGMA tree built, root has 2 children")
    
    leaf_names = root.get_leaf_names()
    assert leaf_names == {'A', 'B', 'C', 'D'}
    print("  ✓ All leaves present in tree")
    
    assert newick.endswith(';')
    assert 'A' in newick and 'B' in newick and 'C' in newick and 'D' in newick
    print(f"  ✓ Newick format: {newick[:60]}...")
    
    print("✓ UPGMA tests passed\n")


def test_neighbor_joining():
    print("Testing Neighbor-Joining algorithm...")
    
    dist_matrix = [
        [0.0, 0.05, 0.10, 0.15],
        [0.05, 0.0, 0.10, 0.15],
        [0.10, 0.10, 0.0, 0.10],
        [0.15, 0.15, 0.10, 0.0],
    ]
    labels = ['W', 'X', 'Y', 'Z']
    
    root, newick = build_tree(dist_matrix, labels, method="nj")
    
    assert root is not None
    assert len(root.children) == 2
    print(f"  ✓ NJ tree built, root has 2 children")
    
    leaf_names = root.get_leaf_names()
    assert leaf_names == {'W', 'X', 'Y', 'Z'}, f"Expected {{'W', 'X', 'Y', 'Z'}}, got {leaf_names}"
    print("  ✓ All leaves present in tree")
    
    for child in root.children:
        assert child.branch_length is not None
        assert child.branch_length >= 0
    print("  ✓ All branch lengths are non-negative")
    
    print(f"  ✓ Newick format: {newick[:60]}...")
    
    print("✓ Neighbor-Joining tests passed\n")


def test_tree_splits_and_rf():
    print("Testing tree splits and Robinson-Foulds distance...")
    
    dist_matrix = [
        [0.0, 0.1, 0.4, 0.5],
        [0.1, 0.0, 0.4, 0.5],
        [0.4, 0.4, 0.0, 0.3],
        [0.5, 0.5, 0.3, 0.0],
    ]
    labels = ['A', 'B', 'C', 'D']
    
    tree1, _ = build_tree(dist_matrix, labels, method="upgma")
    tree2, _ = build_tree(dist_matrix, labels, method="upgma")
    
    splits1 = get_tree_splits(tree1)
    splits2 = get_tree_splits(tree2)
    
    print(f"  ✓ Tree 1 splits: {[sorted(s) for s in splits1]}")
    print(f"  ✓ Tree 2 splits: {[sorted(s) for s in splits2]}")
    
    rf, normalized, matching, total = robinson_foulds_distance(
        splits1, splits2, {'A', 'B', 'C', 'D'}
    )
    assert rf == 0, f"Identical trees should have RF distance 0, got {rf}"
    assert normalized == 0.0
    print("  ✓ Identical trees: RF distance = 0")
    
    dist_matrix2 = [
        [0.0, 0.5, 0.1, 0.4],
        [0.5, 0.0, 0.5, 0.1],
        [0.1, 0.5, 0.0, 0.5],
        [0.4, 0.1, 0.5, 0.0],
    ]
    tree3, _ = build_tree(dist_matrix2, labels, method="upgma")
    splits3 = get_tree_splits(tree3)
    
    rf2, normalized2, matching2, total2 = robinson_foulds_distance(
        splits1, splits3, {'A', 'B', 'C', 'D'}
    )
    print(f"  ✓ Different trees: RF distance = {rf2}, normalized = {normalized2:.2f}")
    
    print("✓ Tree splits and RF distance tests passed\n")


def test_linear_regression():
    print("Testing linear regression...")
    
    x = [2020.0, 2021.0, 2022.0, 2023.0, 2024.0]
    y = [0.0, 0.01, 0.02, 0.03, 0.04]
    
    slope, intercept, r_squared, residuals = linear_regression(x, y)
    
    assert abs(slope - 0.01) < 0.001
    assert r_squared > 0.99
    print(f"  ✓ Slope = {slope:.4f}, Intercept = {intercept:.4f}, R² = {r_squared:.4f}")
    
    x = [1.0, 2.0, 3.0, 4.0, 5.0]
    y = [2.0, 4.0, 6.0, 8.0, 10.0]
    slope, intercept, r_squared, residuals = linear_regression(x, y)
    assert abs(slope - 2.0) < 0.001
    assert abs(intercept) < 0.001
    print(f"  ✓ Perfect linear: y = 2x, slope = {slope:.4f}")
    
    print("✓ Linear regression tests passed\n")


def test_date_to_year():
    print("Testing date to year conversion...")
    
    d1 = datetime(2024, 1, 1)
    y1 = date_to_year(d1)
    assert abs(y1 - 2024.0) < 0.01
    print(f"  ✓ 2024-01-01 → {y1:.4f}")
    
    d2 = datetime(2024, 6, 30)
    y2 = date_to_year(d2)
    assert y2 > 2024.4 and y2 < 2024.6
    print(f"  ✓ 2024-06-30 → {y2:.4f}")
    
    d3 = datetime(2024, 12, 31)
    y3 = date_to_year(d3)
    assert y3 > 2024.9 and y3 < 2025.0
    print(f"  ✓ 2024-12-31 → {y3:.4f}")
    
    print("✓ Date conversion tests passed\n")


def test_root_to_tip_distances():
    print("Testing root-to-tip distances...")
    
    dist_matrix = [
        [0.0, 0.1, 0.2, 0.3],
        [0.1, 0.0, 0.2, 0.3],
        [0.2, 0.2, 0.0, 0.1],
        [0.3, 0.3, 0.1, 0.0],
    ]
    labels = ['S1', 'S2', 'S3', 'S4']
    
    root, _ = build_tree(dist_matrix, labels, method="nj")
    distances = compute_root_to_tip_distances(root)
    
    assert len(distances) == 4
    for name, dist in distances.items():
        assert dist >= 0
        print(f"  ✓ {name}: {dist:.4f}")
    
    print("✓ Root-to-tip distance tests passed\n")


def test_performance_200_samples():
    print("Testing performance with 200 samples...")
    import time
    import random
    
    n = 200
    n_sites = 100
    
    sample_ids = list(range(1, n + 1))
    variants_by_sample = {}
    ref_base_by_pos = {}
    all_positions = list(range(1, n_sites + 1))
    
    bases = ['A', 'T', 'C', 'G']
    for pos in all_positions:
        ref_base_by_pos[pos] = random.choice(bases)
    
    for sid in sample_ids:
        variants_by_sample[sid] = {}
        for pos in all_positions:
            if random.random() < 0.1:
                alt_bases = [b for b in bases if b != ref_base_by_pos[pos]]
                variants_by_sample[sid][pos] = random.choice(alt_bases)
    
    start = time.time()
    seq_matrix = build_sequence_matrix(
        variants_by_sample, all_positions, sample_ids, ref_base_by_pos
    )
    dist_matrix, warnings, names = compute_distance_matrix(seq_matrix, sample_ids)
    elapsed = time.time() - start
    
    assert len(dist_matrix) == n
    assert len(dist_matrix[0]) == n
    assert elapsed < 30, f"Distance matrix for {n} samples took {elapsed:.2f}s, should be < 30s"
    print(f"  ✓ {n} samples, {n_sites} sites: {elapsed:.2f}s (< 30s ✓)")
    
    print("✓ Performance tests passed\n")


if __name__ == "__main__":
    print("=" * 60)
    print("Running Phylogeny Module Tests")
    print("=" * 60 + "\n")
    
    test_kimura_distance()
    test_distance_matrix()
    test_upgma()
    test_neighbor_joining()
    test_tree_splits_and_rf()
    test_linear_regression()
    test_date_to_year()
    test_root_to_tip_distances()
    test_performance_200_samples()
    
    print("=" * 60)
    print("All tests passed! ✓")
    print("=" * 60)
