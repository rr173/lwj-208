from typing import List, Tuple, Dict
from .matrix_ops import (
    matmul, transpose, elementwise_multiply, elementwise_divide,
    frobenius_norm, subtract, create_random_matrix, normalize_columns, normalize_rows
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

    dispersion = 0.0
    total_pairs = 0
    for i in range(n_samples):
        for j in range(i + 1, n_samples):
            c = consensus[i][j]
            dispersion += c * (1 - c)
            total_pairs += 1

    if total_pairs == 0:
        return 1.0

    avg_dispersion = dispersion / total_pairs
    stability = 1.0 - 4.0 * avg_dispersion
    return max(0.0, min(1.0, stability))


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
