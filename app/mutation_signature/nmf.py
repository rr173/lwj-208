from typing import List, Tuple, Dict
from .matrix_ops import (
    matmul, transpose, elementwise_multiply, elementwise_divide,
    frobenius_norm, subtract, create_random_matrix, normalize_columns, normalize_rows,
    dot_product
)


def nmf_multiplicative_update(
    V: List[List[float]],
    k: int,
    max_iter: int = 5000,
    tol: float = 1e-6,
    patience: int = 20,
    seed: int = None,
) -> Tuple[List[List[float]], List[List[float]], float, int]:
    m = len(V)
    n = len(V[0])

    W = create_random_matrix(m, k, seed=seed)
    H = create_random_matrix(k, n, seed=seed + 1 if seed is not None else None)

    prev_error = float('inf')
    stable_count = 0
    iterations = 0

    for iteration in range(max_iter):
        iterations = iteration + 1

        Wt = transpose(W)
        Ht = transpose(H)

        WtV = matmul(Wt, V)
        WtWH = matmul(matmul(Wt, W), H)
        H_update = elementwise_multiply(H, elementwise_divide(WtV, WtWH))
        H = H_update

        VHT = matmul(V, Ht)
        WHHT = matmul(matmul(W, H), Ht)
        W_update = elementwise_multiply(W, elementwise_divide(VHT, WHHT))
        W = W_update

        WH = matmul(W, H)
        diff = subtract(V, WH)
        current_error = frobenius_norm(diff)

        if prev_error > 0:
            relative_change = abs(prev_error - current_error) / prev_error
            if relative_change < tol:
                stable_count += 1
                if stable_count >= patience:
                    break
            else:
                stable_count = 0

        prev_error = current_error

    W, col_sums_W = normalize_columns(W)
    for j in range(n):
        for i in range(k):
            H[i][j] *= col_sums_W[i]

    return W, H, prev_error, iterations


def nmf_with_multiple_initializations(
    V: List[List[float]],
    k: int,
    n_init: int = 10,
    max_iter: int = 5000,
    tol: float = 1e-6,
    patience: int = 20,
    base_seed: int = 42,
) -> Tuple[List[List[float]], List[List[float]], float, int]:
    best_W = None
    best_H = None
    best_error = float('inf')
    best_iterations = 0

    for init_idx in range(n_init):
        seed = base_seed + init_idx * 1000
        W, H, error, iterations = nmf_multiplicative_update(
            V, k, max_iter=max_iter, tol=tol, patience=patience, seed=seed
        )
        if error < best_error:
            best_error = error
            best_W = W
            best_H = H
            best_iterations = iterations

    return best_W, best_H, best_error, best_iterations


def compute_reconstruction_error(
    V: List[List[float]],
    W: List[List[float]],
    H: List[List[float]],
) -> float:
    WH = matmul(W, H)
    diff = subtract(V, WH)
    return frobenius_norm(diff)


def _upgma_cluster(dist_matrix: List[List[float]]) -> List[Tuple[int, int, float]]:
    n = len(dist_matrix)
    if n <= 1:
        return []

    active = [True] * n
    sizes = [1] * n
    cluster_map = {i: i for i in range(n)}
    next_id = n
    merges = []
    dist = [[dist_matrix[i][j] for j in range(n)] for i in range(n)]
    max_id = n
    row_cache = {}

    def get_dist(a: int, b: int) -> float:
        key = (min(a, b), max(a, b))
        if key in row_cache:
            return row_cache[key]
        return dist[a][b] if a < max_id and b < max_id else 0.0

    remaining = list(range(n))

    for step in range(n - 1):
        min_dist = float('inf')
        merge_i = -1
        merge_j = -1
        for idx_a in range(len(remaining)):
            for idx_b in range(idx_a + 1, len(remaining)):
                a = remaining[idx_a]
                b = remaining[idx_b]
                d = get_dist(a, b)
                if d < min_dist:
                    min_dist = d
                    merge_i = a
                    merge_j = b

        new_id = next_id
        next_id += 1
        sizes.append(sizes[merge_i] + sizes[merge_j])
        active.append(True)

        new_dists = {}
        for other in remaining:
            if other == merge_i or other == merge_j:
                continue
            d_i = get_dist(min(merge_i, other), max(merge_i, other))
            d_j = get_dist(min(merge_j, other), max(merge_j, other))
            new_d = (d_i * sizes[merge_i] + d_j * sizes[merge_j]) / (sizes[merge_i] + sizes[merge_j])
            key = (min(new_id, other), max(new_id, other))
            row_cache[key] = new_d

        merges.append((merge_i, merge_j, min_dist))

        remaining = [r for r in remaining if r != merge_i and r != merge_j]
        remaining.append(new_id)

    return merges


def _cophenetic_distances(
    merges: List[Tuple[int, int, float]],
    n_items: int,
) -> List[List[float]]:
    coph = [[0.0] * n_items for _ in range(n_items)]
    cluster_members = {i: [i] for i in range(n_items)}
    next_id = n_items

    for a, b, height in merges:
        members_a = cluster_members.get(a, [a])
        members_b = cluster_members.get(b, [b])
        for i in members_a:
            for j in members_b:
                coph[i][j] = height
                coph[j][i] = height
        combined = members_a + members_b
        cluster_members[next_id] = combined
        next_id += 1

    return coph


def _pearson_correlation(x: List[float], y: List[float]) -> float:
    n = len(x)
    if n < 2:
        return 0.0

    mean_x = sum(x) / n
    mean_y = sum(y) / n

    num = 0.0
    denom_x = 0.0
    denom_y = 0.0
    for i in range(n):
        dx = x[i] - mean_x
        dy = y[i] - mean_y
        num += dx * dy
        denom_x += dx * dx
        denom_y += dy * dy

    if denom_x == 0.0 or denom_y == 0.0:
        return 0.0

    return num / ((denom_x ** 0.5) * (denom_y ** 0.5))


def cophenetic_correlation(
    connectivity_matrices: List[List[List[float]]]
) -> float:
    n_runs = len(connectivity_matrices)
    if n_runs < 2:
        return 1.0

    n_samples = len(connectivity_matrices[0])

    consensus = [[0.0 for _ in range(n_samples)] for _ in range(n_samples)]
    for i in range(n_samples):
        for j in range(n_samples):
            for cm in connectivity_matrices:
                consensus[i][j] += cm[i][j]
            consensus[i][j] /= n_runs

    dist = [[0.0 for _ in range(n_samples)] for _ in range(n_samples)]
    for i in range(n_samples):
        for j in range(n_samples):
            dist[i][j] = 1.0 - consensus[i][j]

    merges = _upgma_cluster(dist)
    coph_dist = _cophenetic_distances(merges, n_samples)

    original_vals = []
    coph_vals = []
    for i in range(n_samples):
        for j in range(i + 1, n_samples):
            original_vals.append(dist[i][j])
            coph_vals.append(coph_dist[i][j])

    return _pearson_correlation(original_vals, coph_vals)


def sample_connectivity_matrix(
    H: List[List[float]]
) -> List[List[float]]:
    k = len(H)
    n_samples = len(H[0])

    assignments = []
    for j in range(n_samples):
        max_val = -1.0
        max_idx = 0
        for i in range(k):
            if H[i][j] > max_val:
                max_val = H[i][j]
                max_idx = i
        assignments.append(max_idx)

    n = n_samples
    connectivity = [[0.0 for _ in range(n)] for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if assignments[i] == assignments[j]:
                connectivity[i][j] = 1.0

    return connectivity


def find_optimal_k(
    V: List[List[float]],
    k_min: int = 2,
    k_max: int = 8,
    n_runs_per_k: int = 10,
    max_iter: int = 5000,
    tol: float = 1e-6,
    patience: int = 20,
) -> Tuple[Dict[int, float], Dict[int, float], int]:
    reconstruction_errors = {}
    cophenetic_correlations = {}

    for k in range(k_min, k_max + 1):
        connectivity_matrices = []
        best_error_for_k = float('inf')

        for run in range(n_runs_per_k):
            seed = k * 10000 + run
            W, H, error, iterations = nmf_multiplicative_update(
                V, k, max_iter=max_iter, tol=tol, patience=patience, seed=seed
            )

            if error < best_error_for_k:
                best_error_for_k = error

            cm = sample_connectivity_matrix(H)
            connectivity_matrices.append(cm)

        reconstruction_errors[k] = best_error_for_k
        cophenetic_correlations[k] = cophenetic_correlation(connectivity_matrices)

    recommended_k = k_min
    for k in range(k_min + 1, k_max + 1):
        if cophenetic_correlations[k] >= 0.8:
            recommended_k = k
        else:
            if cophenetic_correlations[k] < cophenetic_correlations[k - 1] - 0.1:
                break
            recommended_k = k

    return reconstruction_errors, cophenetic_correlations, recommended_k
