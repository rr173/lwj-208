from typing import List, Tuple
import random


def create_matrix(rows: int, cols: int, fill_value: float = 0.0) -> List[List[float]]:
    return [[fill_value for _ in range(cols)] for _ in range(rows)]


def create_random_matrix(rows: int, cols: int, seed: int = None) -> List[List[float]]:
    if seed is not None:
        rng = random.Random(seed)
    else:
        rng = random.Random()
    return [[rng.random() for _ in range(cols)] for _ in range(rows)]


def matmul(A: List[List[float]], B: List[List[float]]) -> List[List[float]]:
    m = len(A)
    n = len(A[0])
    p = len(B[0])
    if len(B) != n:
        raise ValueError(f"Matrix dimensions incompatible: {m}x{n} * {len(B)}x{p}")
    result = create_matrix(m, p)
    for i in range(m):
        for j in range(p):
            s = 0.0
            for k in range(n):
                s += A[i][k] * B[k][j]
            result[i][j] = s
    return result


def transpose(A: List[List[float]]) -> List[List[float]]:
    m = len(A)
    n = len(A[0])
    result = create_matrix(n, m)
    for i in range(m):
        for j in range(n):
            result[j][i] = A[i][j]
    return result


def elementwise_multiply(A: List[List[float]], B: List[List[float]]) -> List[List[float]]:
    m = len(A)
    n = len(A[0])
    result = create_matrix(m, n)
    for i in range(m):
        for j in range(n):
            result[i][j] = A[i][j] * B[i][j]
    return result


def elementwise_divide(A: List[List[float]], B: List[List[float]], eps: float = 1e-10) -> List[List[float]]:
    m = len(A)
    n = len(A[0])
    result = create_matrix(m, n)
    for i in range(m):
        for j in range(n):
            result[i][j] = A[i][j] / (B[i][j] + eps)
    return result


def frobenius_norm(A: List[List[float]]) -> float:
    s = 0.0
    for row in A:
        for val in row:
            s += val * val
    return s ** 0.5


def subtract(A: List[List[float]], B: List[List[float]]) -> List[List[float]]:
    m = len(A)
    n = len(A[0])
    result = create_matrix(m, n)
    for i in range(m):
        for j in range(n):
            result[i][j] = A[i][j] - B[i][j]
    return result


def sum_columns(A: List[List[float]]) -> List[float]:
    n = len(A[0])
    result = [0.0] * n
    for row in A:
        for j in range(n):
            result[j] += row[j]
    return result


def normalize_columns(A: List[List[float]]) -> Tuple[List[List[float]], List[float]]:
    col_sums = sum_columns(A)
    m = len(A)
    n = len(A[0])
    result = create_matrix(m, n)
    for i in range(m):
        for j in range(n):
            result[i][j] = A[i][j] / (col_sums[j] + 1e-10)
    return result, col_sums


def normalize_rows(A: List[List[float]]) -> List[List[float]]:
    m = len(A)
    n = len(A[0])
    result = create_matrix(m, n)
    for i in range(m):
        row_sum = sum(A[i])
        for j in range(n):
            result[i][j] = A[i][j] / (row_sum + 1e-10)
    return result


def normalize_vector(v: List[float]) -> List[float]:
    s = sum(v)
    if s == 0:
        return v.copy()
    return [x / s for x in v]


def dot_product(a: List[float], b: List[float]) -> float:
    if len(a) != len(b):
        raise ValueError("Vectors must have same length")
    s = 0.0
    for i in range(len(a)):
        s += a[i] * b[i]
    return s


def vector_norm(v: List[float]) -> float:
    s = 0.0
    for x in v:
        s += x * x
    return s ** 0.5


def cosine_similarity(a: List[float], b: List[float]) -> float:
    norm_a = vector_norm(a)
    norm_b = vector_norm(b)
    if norm_a == 0 or norm_b == 0:
        return 0.0
    return dot_product(a, b) / (norm_a * norm_b)


def hstack(A: List[List[float]], B: List[List[float]]) -> List[List[float]]:
    if len(A) != len(B):
        raise ValueError("Matrices must have same number of rows")
    m = len(A)
    n1 = len(A[0])
    n2 = len(B[0])
    result = create_matrix(m, n1 + n2)
    for i in range(m):
        for j in range(n1):
            result[i][j] = A[i][j]
        for j in range(n2):
            result[i][n1 + j] = B[i][j]
    return result
