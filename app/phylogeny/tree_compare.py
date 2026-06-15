from typing import List, Dict, Tuple, Set
from .tree_building import TreeNode, get_tree_splits


def collect_branch_lengths(node: TreeNode) -> Dict[frozenset, float]:
    """
    Collect branch lengths keyed by the split they represent.
    
    Returns:
        {split: branch_length}
    """
    branch_lengths = {}
    all_leaves = node.get_leaf_names()
    
    def traverse(n: TreeNode):
        if n.is_leaf:
            return {n.name} if n.name else set()
        
        leaves = set()
        for child in n.children:
            child_leaves = traverse(child)
            leaves.update(child_leaves)
            if not child.is_leaf:
                child_all = child.get_leaf_names()
                split = frozenset(child_all)
                complement = all_leaves - child_all
                if len(split) > 1 and len(complement) > 1:
                    if split < complement:
                        key = split
                    else:
                        key = frozenset(complement)
                    if key not in branch_lengths or (child.branch_length is not None and child.branch_length > branch_lengths.get(key, 0)):
                        if child.branch_length is not None:
                            branch_lengths[key] = child.branch_length
        
        return leaves
    
    traverse(node)
    return branch_lengths


def weighted_branch_length_distance(
    branch_lengths_a: Dict[frozenset, float],
    branch_lengths_b: Dict[frozenset, float]
) -> float:
    """
    Calculate the weighted branch length distance between two trees.
    
    This is the sum of absolute differences of branch lengths for matching splits
    plus the sum of branch lengths for splits that exist in only one tree.
    """
    all_splits = set(branch_lengths_a.keys()) | set(branch_lengths_b.keys())
    
    total = 0.0
    for s in all_splits:
        len_a = branch_lengths_a.get(s, 0.0)
        len_b = branch_lengths_b.get(s, 0.0)
        total += abs(len_a - len_b)
    
    return total


def compare_trees(
    tree_a: TreeNode,
    tree_b: TreeNode,
    sample_names_a: List[str],
    sample_names_b: List[str]
) -> Dict:
    """
    Compare two phylogenetic trees.
    
    Args:
        tree_a: First tree
        tree_b: Second tree
        sample_names_a: Sample names in tree_a
        sample_names_b: Sample names in tree_b
    
    Returns:
        Dictionary with comparison results
    """
    leaves_a = set(sample_names_a)
    leaves_b = set(sample_names_b)
    common_leaves = leaves_a & leaves_b
    
    if len(common_leaves) < 3:
        return {
            "rf_distance": -1,
            "normalized_rf_distance": 1.0,
            "matching_splits": 0,
            "total_splits": 0,
            "weight_dist": 0.0,
            "inconsistent_branches": [],
            "warning": f"Trees share only {len(common_leaves)} samples. Need at least 3 for meaningful comparison.",
        }
    
    splits_a = get_tree_splits(tree_a)
    splits_b = get_tree_splits(tree_b)
    
    splits_a_common = set()
    for s in splits_a:
        common = s & common_leaves
        complement = common_leaves - common
        if len(common) > 1 and len(complement) > 1:
            if common < complement:
                splits_a_common.add(frozenset(common))
            else:
                splits_a_common.add(frozenset(complement))
    
    splits_b_common = set()
    for s in splits_b:
        common = s & common_leaves
        complement = common_leaves - common
        if len(common) > 1 and len(complement) > 1:
            if common < complement:
                splits_b_common.add(frozenset(common))
            else:
                splits_b_common.add(frozenset(complement))
    
    matching = splits_a_common & splits_b_common
    only_a = splits_a_common - splits_b_common
    only_b = splits_b_common - splits_a_common
    
    rf = len(only_a) + len(only_b)
    total = len(splits_a_common) + len(splits_b_common)
    normalized = rf / total if total > 0 else 0.0
    
    bl_a = collect_branch_lengths(tree_a)
    bl_b = collect_branch_lengths(tree_b)
    weight_dist = weighted_branch_length_distance(bl_a, bl_b)
    
    inconsistent = []
    for s in only_a:
        comp = common_leaves - s
        inconsistent.append({
            "split_a": sorted(list(s)),
            "split_b": sorted(list(comp))
        })
    for s in only_b:
        comp = common_leaves - s
        inconsistent.append({
            "split_a": sorted(list(comp)),
            "split_b": sorted(list(s))
        })
    
    return {
        "rf_distance": rf,
        "normalized_rf_distance": normalized,
        "matching_splits": len(matching),
        "total_splits": len(splits_a_common),
        "weight_dist": weight_dist,
        "inconsistent_branches": inconsistent,
    }
