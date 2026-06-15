import random
from typing import List, Dict, Set, Callable, Optional
from .distance import kimura_two_parameter_distance
from .tree_building import TreeNode, build_tree, get_tree_splits


def bootstrap_resample_sites(
    sequence_matrix: Dict[int, List[str]],
    sample_ids: List[int]
) -> Dict[int, List[str]]:
    """
    Perform bootstrap resampling of variant sites.
    Randomly sample sites with replacement to create a new sequence matrix.
    """
    n_sites = len(sequence_matrix[sample_ids[0]]) if sample_ids else 0
    if n_sites == 0:
        return sequence_matrix
    
    sampled_indices = [random.randint(0, n_sites - 1) for _ in range(n_sites)]
    
    resampled = {}
    for sid in sample_ids:
        seq = sequence_matrix[sid]
        resampled[sid] = [seq[i] for i in sampled_indices]
    
    return resampled


def compute_distance_matrix_from_matrix(
    sequence_matrix: Dict[int, List[str]],
    sample_ids: List[int]
) -> List[List[float]]:
    """Compute distance matrix from a sequence matrix."""
    n = len(sample_ids)
    dist_matrix = [[0.0] * n for _ in range(n)]
    
    for i in range(n):
        for j in range(i + 1, n):
            seq_i = sequence_matrix[sample_ids[i]]
            seq_j = sequence_matrix[sample_ids[j]]
            d, _ = kimura_two_parameter_distance(seq_i, seq_j)
            dist_matrix[i][j] = d
            dist_matrix[j][i] = d
    
    return dist_matrix


def compute_bootstrap_support(
    original_tree: TreeNode,
    sequence_matrix: Dict[int, List[str]],
    sample_ids: List[int],
    labels: List[str],
    n_replicates: int = 100,
    method: str = "nj",
    progress_callback: Optional[Callable[[int], None]] = None
) -> TreeNode:
    """
    Compute bootstrap support values for each internal node.
    
    Args:
        original_tree: The original tree built from full data
        sequence_matrix: Original sequence matrix
        sample_ids: List of sample IDs in order
        labels: Sample names in order
        n_replicates: Number of bootstrap replicates (default 100)
        method: Tree building method ("nj" or "upgma")
        progress_callback: Optional callback for progress updates
    
    Returns:
        Original tree with bootstrap_support values set on internal nodes
    """
    if n_replicates <= 0:
        return original_tree
    
    original_splits = get_tree_splits(original_tree)
    split_counts: Dict[frozenset, int] = {s: 0 for s in original_splits}
    
    for rep in range(n_replicates):
        resampled_matrix = bootstrap_resample_sites(sequence_matrix, sample_ids)
        dist_matrix = compute_distance_matrix_from_matrix(resampled_matrix, sample_ids)
        bootstrap_tree, _ = build_tree(dist_matrix, labels, method)
        bootstrap_splits = get_tree_splits(bootstrap_tree)
        
        for s in original_splits:
            if s in bootstrap_splits:
                split_counts[s] += 1
        
        if progress_callback and (rep + 1) % 10 == 0:
            progress_callback(rep + 1)
    
    def add_support(node: TreeNode):
        if node.is_leaf:
            return
        
        node_leaves = node.get_leaf_names()
        all_leaves = original_tree.get_leaf_names()
        complement = all_leaves - node_leaves
        
        if len(node_leaves) > 1 and len(complement) > 1:
            if frozenset(node_leaves) < frozenset(complement):
                key = frozenset(node_leaves)
            else:
                key = frozenset(complement)
            
            if key in split_counts:
                node.bootstrap_support = (split_counts[key] / n_replicates) * 100
            else:
                for s in original_splits:
                    if s == key or s == frozenset(complement):
                        node.bootstrap_support = (split_counts[s] / n_replicates) * 100
                        break
        
        for child in node.children:
            add_support(child)
    
    add_support(original_tree)
    return original_tree
