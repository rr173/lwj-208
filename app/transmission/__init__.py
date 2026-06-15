from typing import List, Dict, Set, Tuple, TYPE_CHECKING

if TYPE_CHECKING:
    from app.phylogeny.tree_building import TreeNode


def jaccard_similarity(set_a: Set, set_b: Set) -> float:
    if not set_a and not set_b:
        return 1.0
    intersection = set_a & set_b
    union = set_a | set_b
    return len(intersection) / len(union) if union else 0.0


def compute_jaccard_matrix(
    variant_sets: Dict[int, Set[Tuple]],
    sample_ids: List[int],
) -> List[List[float]]:
    n = len(sample_ids)
    matrix = [[0.0] * n for _ in range(n)]
    for i in range(n):
        matrix[i][i] = 1.0
        for j in range(i + 1, n):
            sim = jaccard_similarity(
                variant_sets.get(sample_ids[i], set()),
                variant_sets.get(sample_ids[j], set()),
            )
            matrix[i][j] = sim
            matrix[j][i] = sim
    return matrix


def jaccard_to_distance(jaccard: float) -> float:
    return 1.0 - jaccard


def compute_distance_matrix_from_jaccard(
    jaccard_matrix: List[List[float]],
) -> List[List[float]]:
    n = len(jaccard_matrix)
    dist = [[0.0] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            dist[i][j] = jaccard_to_distance(jaccard_matrix[i][j])
    return dist


def upgma_cluster(
    dist_matrix: List[List[float]],
    labels: List[str],
    threshold: float = 0.7,
) -> List[List[int]]:
    from app.phylogeny.tree_building import upgma

    n = len(labels)
    if n == 0:
        return []
    if n == 1:
        return [[0]]

    root = upgma(dist_matrix, labels)

    def get_max_height(node) -> float:
        if node.is_leaf:
            return 0.0
        h = 0.0
        for child in node.children:
            child_h = get_max_height(child)
            bl = child.branch_length if child.branch_length is not None else 0.0
            h = max(h, child_h + bl)
        return h

    def cut_at_height(node, cut_height: float) -> List[List[int]]:
        if node.is_leaf:
            if node.name is not None:
                try:
                    idx = labels.index(node.name)
                    return [[idx]]
                except ValueError:
                    return []
            return []

        node_h = get_max_height(node)
        if node_h <= cut_height:
            leaves = _collect_leaf_indices(node, labels)
            if leaves:
                return [leaves]

        result = []
        for child in node.children:
            result.extend(cut_at_height(child, cut_height))
        return result

    max_h = get_max_height(root)
    dist_threshold = jaccard_to_distance(threshold)
    cut_h = dist_threshold / 2.0
    if cut_h > max_h:
        all_leaves = list(range(n))
        return [all_leaves]

    clusters = cut_at_height(root, cut_h)

    clusters = [c for c in clusters if len(c) > 0]

    covered = set()
    for c in clusters:
        covered.update(c)
    for i in range(n):
        if i not in covered:
            clusters.append([i])

    return clusters


def _collect_leaf_indices(node, labels: List[str]) -> List[int]:
    if node.is_leaf:
        if node.name is not None:
            try:
                idx = labels.index(node.name)
                return [idx]
            except ValueError:
                return []
        return []
    result = []
    for child in node.children:
        result.extend(_collect_leaf_indices(child, labels))
    return result


def transitive_reduction(
    edges: List[Tuple[int, int]],
    nodes: Set[int],
) -> List[Tuple[int, int]]:
    if not edges:
        return []

    adj: Dict[int, Set[int]] = {n: set() for n in nodes}
    for src, tgt in edges:
        adj[src].add(tgt)

    descendants: Dict[int, Set[int]] = {}

    def get_descendants(node: int, visited: Set[int]) -> Set[int]:
        if node in descendants:
            return descendants[node]
        if node in visited:
            return set()
        visited.add(node)
        desc = set()
        for child in adj[node]:
            desc.add(child)
            desc.update(get_descendants(child, visited))
        descendants[node] = desc
        return desc

    all_nodes = set(nodes)
    for n in all_nodes:
        get_descendants(n, set())

    reduced = []
    for src, tgt in edges:
        indirect = False
        for intermediate in adj[src]:
            if intermediate != tgt and tgt in descendants.get(intermediate, set()):
                indirect = True
                break
        if not indirect:
            reduced.append((src, tgt))

    return reduced
