from typing import List, Dict, Tuple, Optional, Set
import copy
import math


class TreeNode:
    def __init__(self, name: Optional[str] = None):
        self.name = name
        self.children: List['TreeNode'] = []
        self.branch_length: Optional[float] = None
        self.bootstrap_support: Optional[float] = None
        self.divergence_time: Optional[float] = None
        self.is_leaf = False
        self.id = None
        self.parent: Optional['TreeNode'] = None
        self.leaf_names: Set[str] = set()

    def add_child(self, child: 'TreeNode', branch_length: float):
        child.branch_length = branch_length
        child.parent = self
        self.children.append(child)

    def get_leaf_names(self) -> Set[str]:
        if self.is_leaf:
            return {self.name} if self.name else set()
        if not self.leaf_names:
            names = set()
            for child in self.children:
                names.update(child.get_leaf_names())
            self.leaf_names = names
        return self.leaf_names

    def to_dict(self) -> Dict:
        return {
            "name": self.name,
            "children": [c.to_dict() for c in self.children],
            "branch_length": self.branch_length,
            "bootstrap_support": self.bootstrap_support,
            "divergence_time": self.divergence_time,
            "is_leaf": self.is_leaf,
        }


def _to_newick(node: TreeNode, parent_branch_length: Optional[float] = None) -> str:
    bl = node.branch_length if node.branch_length is not None else parent_branch_length
    
    if node.is_leaf:
        name = node.name if node.name else ""
        if bl is not None:
            return f"{name}:{bl:.6f}"
        return name
    
    parts = []
    for child in node.children:
        parts.append(_to_newick(child, bl))
    
    inner = ",".join(parts)
    result = f"({inner})"
    
    if not node.is_leaf and node.bootstrap_support is not None:
        result += f"{node.bootstrap_support:.0f}"
    
    if bl is not None and node.parent is not None:
        result += f":{bl:.6f}"
    
    return result


def to_newick_format(root: TreeNode) -> str:
    return _to_newick(root) + ";"


def upgma(dist_matrix: List[List[float]], labels: List[str]) -> TreeNode:
    """
    Unweighted Pair Group Method with Arithmetic Mean (UPGMA) algorithm.
    Builds an ultrametric tree.
    """
    n = len(labels)
    D = [row[:] for row in dist_matrix]
    node_list = [TreeNode(labels[i]) for i in range(n)]
    for node in node_list:
        node.is_leaf = True
    
    cluster_sizes = [1] * n
    heights = [0.0] * n
    
    for _ in range(n - 1):
        m = len(node_list)
        if m < 2:
            break
            
        min_dist = float('inf')
        i_min, j_min = -1, -1
        
        for i in range(m):
            for j in range(i + 1, m):
                if not math.isnan(D[i][j]) and D[i][j] < min_dist:
                    min_dist = D[i][j]
                    i_min, j_min = i, j
        
        if i_min == -1:
            for i in range(m):
                for j in range(i + 1, m):
                    if not math.isnan(D[i][j]):
                        min_dist = D[i][j]
                        i_min, j_min = i, j
                        break
                if i_min != -1:
                    break
        
        if i_min == -1:
            break
        
        new_node = TreeNode()
        
        target_height = min_dist / 2.0
        bl_i = max(0.0, target_height - heights[i_min])
        bl_j = max(0.0, target_height - heights[j_min])
        
        new_node.add_child(node_list[i_min], bl_i)
        new_node.add_child(node_list[j_min], bl_j)
        
        new_size = cluster_sizes[i_min] + cluster_sizes[j_min]
        new_height = target_height
        
        new_D = []
        for p in range(m):
            if p == i_min or p == j_min:
                continue
            new_row = []
            for q in range(m):
                if q == i_min or q == j_min:
                    continue
                new_row.append(D[p][q])
            d_p_new = (cluster_sizes[i_min] * D[p][i_min] + cluster_sizes[j_min] * D[p][j_min]) / new_size
            new_row.append(d_p_new)
            new_D.append(new_row)
        
        if new_D:
            last_row = [row[-1] for row in new_D]
            last_row.append(0.0)
            new_D.append(last_row)
        else:
            new_D = [[0.0]]
        
        D = new_D
        
        new_node_list = []
        new_sizes = []
        new_heights = []
        for p in range(m):
            if p != i_min and p != j_min:
                new_node_list.append(node_list[p])
                new_sizes.append(cluster_sizes[p])
                new_heights.append(heights[p])
        new_node_list.append(new_node)
        new_sizes.append(new_size)
        new_heights.append(new_height)
        
        node_list = new_node_list
        cluster_sizes = new_sizes
        heights = new_heights
    
    return node_list[0] if node_list else TreeNode()


def _fix_upgma_branch_lengths(node: TreeNode, total_dist: float):
    """Fix branch lengths for UPGMA to ensure ultrametric property."""
    if len(node.children) != 2:
        return
    
    child1, child2 = node.children
    
    def get_height(n: TreeNode) -> float:
        if n.is_leaf:
            return 0.0
        if n.branch_length is not None:
            return n.branch_length + get_height(n.children[0])
        return 0.0
    
    h1 = get_height(child1)
    h2 = get_height(child2)
    
    target = total_dist / 2.0
    bl1 = max(0.0, target - h1)
    bl2 = max(0.0, target - h2)
    
    child1.branch_length = bl1
    child2.branch_length = bl2


def neighbor_joining(dist_matrix: List[List[float]], labels: List[str]) -> TreeNode:
    """
    Neighbor-Joining algorithm for phylogenetic tree reconstruction.
    Does not assume a molecular clock.
    """
    n = len(labels)
    D = [row[:] for row in dist_matrix]
    node_list = [TreeNode(labels[i]) for i in range(n)]
    for node in node_list:
        node.is_leaf = True
    
    while len(node_list) > 2:
        m = len(node_list)
        
        r = [0.0] * m
        for i in range(m):
            valid_sum = 0.0
            count = 0
            for j in range(m):
                if j != i and not math.isnan(D[i][j]):
                    valid_sum += D[i][j]
                    count += 1
            if count > 0:
                r[i] = valid_sum / (m - 2) if (m - 2) > 0 else 0.0
        
        min_Q = float('inf')
        i_min, j_min = -1, -1
        
        for i in range(m):
            for j in range(i + 1, m):
                if math.isnan(D[i][j]):
                    continue
                Q = D[i][j] - r[i] - r[j]
                if Q < min_Q:
                    min_Q = Q
                    i_min, j_min = i, j
        
        if i_min == -1:
            for i in range(m):
                for j in range(i + 1, m):
                    if not math.isnan(D[i][j]):
                        i_min, j_min = i, j
                        break
                if i_min != -1:
                    break
        
        if i_min == -1:
            break
        
        new_node = TreeNode()
        
        d_ij = D[i_min][j_min]
        bl_i = (d_ij + r[i_min] - r[j_min]) / 2.0
        bl_j = (d_ij + r[j_min] - r[i_min]) / 2.0
        
        bl_i = max(0.0, bl_i)
        bl_j = max(0.0, bl_j)
        
        new_node.add_child(node_list[i_min], bl_i)
        new_node.add_child(node_list[j_min], bl_j)
        
        new_D = []
        for p in range(m):
            if p == i_min or p == j_min:
                continue
            new_row = []
            for q in range(m):
                if q == i_min or q == j_min:
                    continue
                new_row.append(D[p][q])
            d_p_new = (D[p][i_min] + D[p][j_min] - d_ij) / 2.0
            new_row.append(d_p_new)
            new_D.append(new_row)
        
        last_row = [row[-1] for row in new_D]
        last_row.append(0.0)
        new_D.append(last_row)
        
        D = new_D
        
        new_node_list = []
        for p in range(m):
            if p != i_min and p != j_min:
                new_node_list.append(node_list[p])
        new_node_list.append(new_node)
        node_list = new_node_list
    
    if len(node_list) == 2:
        root = TreeNode()
        d = D[0][1] if not math.isnan(D[0][1]) else 0.0
        bl = max(0.0, d / 2.0)
        root.add_child(node_list[0], bl)
        root.add_child(node_list[1], bl)
        return root
    
    return node_list[0] if node_list else TreeNode()


def build_tree(
    dist_matrix: List[List[float]],
    labels: List[str],
    method: str = "nj"
) -> Tuple[TreeNode, str]:
    """
    Build a phylogenetic tree using the specified method.
    
    Args:
        dist_matrix: N x N distance matrix
        labels: list of sample names
        method: "nj" or "upgma"
    
    Returns:
        (root_node, newick_string)
    """
    if method == "upgma":
        root = upgma(dist_matrix, labels)
    else:
        root = neighbor_joining(dist_matrix, labels)
    
    newick = to_newick_format(root)
    return root, newick


def get_tree_splits(root: TreeNode) -> Set[frozenset]:
    """
    Get all splits (bipartitions) from a rooted or unrooted tree.
    Each split is represented as a frozenset of leaf names.
    """
    splits = set()
    
    def traverse(node: TreeNode):
        if node.is_leaf:
            return {node.name}
        
        leaves = set()
        for child in node.children:
            child_leaves = traverse(child)
            leaves.update(child_leaves)
            if not child.is_leaf:
                child_all = child.get_leaf_names()
                splits.add(frozenset(child_all))
        
        return leaves
    
    all_leaves = traverse(root)
    
    nontrivial_splits = set()
    for s in splits:
        complement = all_leaves - s
        if len(s) > 1 and len(complement) > 1:
            if frozenset(s) < frozenset(complement):
                nontrivial_splits.add(frozenset(s))
            else:
                nontrivial_splits.add(frozenset(complement))
    
    return nontrivial_splits


def count_matching_splits(
    splits_a: Set[frozenset],
    splits_b: Set[frozenset]
) -> int:
    """Count the number of matching splits between two sets."""
    return len(splits_a & splits_b)


def robinson_foulds_distance(
    splits_a: Set[frozenset],
    splits_b: Set[frozenset],
    all_leaves: Set[str]
) -> Tuple[int, float, int, int]:
    """
    Calculate Robinson-Foulds distance between two trees.
    
    Returns:
        (rf_distance, normalized_rf, matching_splits, total_splits)
    """
    a_only = splits_a - splits_b
    b_only = splits_b - splits_a
    rf = len(a_only) + len(b_only)
    
    total = len(splits_a) + len(splits_b)
    normalized = rf / total if total > 0 else 0.0
    
    matching = count_matching_splits(splits_a, splits_b)
    
    return rf, normalized, matching, len(splits_a)
