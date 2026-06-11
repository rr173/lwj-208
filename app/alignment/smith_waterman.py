import hashlib
import time
from typing import Tuple, List, Dict, Optional


MATCH_SCORE = 2
MISMATCH_SCORE = -1
GAP_SCORE = -2

ALIGNMENT_TIMEOUT_SECONDS = 60


class AlignmentTimeoutError(Exception):
    """Raised when alignment takes longer than the timeout."""
    pass


def sequence_hash(seq: str) -> str:
    return hashlib.md5(seq.encode()).hexdigest()


def smith_waterman(query: str, reference: str, timeout: int = ALIGNMENT_TIMEOUT_SECONDS) -> Dict:
    """
    Smith-Waterman local alignment algorithm (optimized).
    Returns alignment result with score, positions, cigar, and alignment strings.
    Raises AlignmentTimeoutError if computation exceeds timeout seconds.
    """
    n = len(query)
    m = len(reference)

    if n == 0 or m == 0:
        return {
            "score": 0,
            "query_start": 0,
            "query_end": 0,
            "ref_start": 0,
            "ref_end": 0,
            "cigar": "",
            "alignment_query": "",
            "alignment_match": "",
            "alignment_ref": "",
        }

    query_bytes = query.encode('ascii')
    ref_bytes = reference.encode('ascii')

    score_prev = [0] * (m + 1)
    score_curr = [0] * (m + 1)
    traceback = [[0] * (m + 1) for _ in range(n + 1)]

    max_score = 0
    max_i, max_j = 0, 0

    match_score = MATCH_SCORE
    mismatch_score = MISMATCH_SCORE
    gap_score = GAP_SCORE

    start_time = time.time()
    check_interval = max(1, n // 10)

    for i in range(1, n + 1):
        if i % check_interval == 0:
            if time.time() - start_time > timeout:
                raise AlignmentTimeoutError(
                    f"Alignment timed out after {timeout} seconds "
                    f"({n} x {m} = {n*m:,} cells)"
                )

        qi = query_bytes[i - 1]
        score_curr[0] = 0
        diag = 0

        for j in range(1, m + 1):
            if qi == ref_bytes[j - 1]:
                match_val = diag + match_score
            else:
                match_val = diag + mismatch_score

            delete_val = score_prev[j] + gap_score
            insert_val = score_curr[j - 1] + gap_score

            current_max = 0
            tb = 0

            if match_val > current_max:
                current_max = match_val
                tb = 1
            if delete_val > current_max:
                current_max = delete_val
                tb = 2
            if insert_val > current_max:
                current_max = insert_val
                tb = 3

            diag = score_prev[j]
            score_curr[j] = current_max
            traceback[i][j] = tb

            if current_max > max_score:
                max_score = current_max
                max_i, max_j = i, j

        score_prev, score_curr = score_curr, score_prev

    if max_score == 0:
        return {
            "score": 0,
            "query_start": 0,
            "query_end": 0,
            "ref_start": 0,
            "ref_end": 0,
            "cigar": "",
            "alignment_query": "",
            "alignment_match": "",
            "alignment_ref": "",
        }

    align_query = []
    align_match = []
    align_ref = []

    i, j = max_i, max_j

    while i > 0 and j > 0 and traceback[i][j] != 0:
        tb = traceback[i][j]
        if tb == 1:
            qb = query_bytes[i - 1]
            rb = ref_bytes[j - 1]
            align_query.append(chr(qb))
            align_ref.append(chr(rb))
            align_match.append("|" if qb == rb else "*")
            i -= 1
            j -= 1
        elif tb == 2:
            align_query.append(chr(query_bytes[i - 1]))
            align_ref.append("-")
            align_match.append(" ")
            i -= 1
        elif tb == 3:
            align_query.append("-")
            align_ref.append(chr(ref_bytes[j - 1]))
            align_match.append(" ")
            j -= 1
        else:
            break

    align_query.reverse()
    align_match.reverse()
    align_ref.reverse()

    cigar = _build_cigar(align_query, align_ref)

    return {
        "score": max_score,
        "query_start": i,
        "query_end": max_i - 1,
        "ref_start": j,
        "ref_end": max_j - 1,
        "cigar": cigar,
        "alignment_query": "".join(align_query),
        "alignment_match": "".join(align_match),
        "alignment_ref": "".join(align_ref),
    }


def _build_cigar(align_query: List[str], align_ref: List[str]) -> str:
    """Build CIGAR string from alignment."""
    cigar_parts = []
    current_op = None
    count = 0

    for q, r in zip(align_query, align_ref):
        if q == "-":
            op = "D"
        elif r == "-":
            op = "I"
        else:
            op = "M"

        if op == current_op:
            count += 1
        else:
            if current_op is not None:
                cigar_parts.append(f"{count}{current_op}")
            current_op = op
            count = 1

    if current_op is not None:
        cigar_parts.append(f"{count}{current_op}")

    return "".join(cigar_parts)


def find_seed_positions(query: str, reference: str, seed_size: int = 15, max_seeds: int = 10) -> List[int]:
    """Find seed positions in reference using k-mer matching for optimization."""
    n = len(query)
    m = len(reference)

    if n < seed_size:
        return [0]

    step = max(1, (n - seed_size + 1) // max_seeds)
    query_kmers = []
    seen = set()
    for i in range(0, n - seed_size + 1, step):
        kmer = query[i:i + seed_size]
        if kmer not in seen:
            seen.add(kmer)
            query_kmers.append((i, kmer))
        if len(query_kmers) >= max_seeds:
            break

    seed_hits = []
    for query_offset, kmer in query_kmers:
        start = 0
        while True:
            pos = reference.find(kmer, start)
            if pos == -1:
                break
            seed_hits.append(pos - query_offset)
            start = pos + 1

    if not seed_hits:
        return []

    seed_hits.sort()

    clusters = []
    current_cluster_start = seed_hits[0]
    current_cluster_end = seed_hits[0]
    cluster_range = n * 2

    for pos in seed_hits[1:]:
        if pos - current_cluster_end <= cluster_range:
            current_cluster_end = pos
        else:
            clusters.append((current_cluster_start, current_cluster_end))
            current_cluster_start = pos
            current_cluster_end = pos
    clusters.append((current_cluster_start, current_cluster_end))

    result = []
    for start, end in clusters:
        result.append((start + end) // 2)

    return result


def optimized_smith_waterman(query: str, reference: str, timeout: int = ALIGNMENT_TIMEOUT_SECONDS) -> Dict:
    """
    Optimized Smith-Waterman using seed finding + banded DP around seeds.
    Falls back to full DP for short sequences.
    """
    n = len(query)
    m = len(reference)

    if n == 0 or m == 0:
        return smith_waterman(query, reference, timeout)

    if n * m <= 5000000:
        return smith_waterman(query, reference, timeout)

    start_time = time.time()

    seed_size = min(20, max(11, n // 20))
    seed_positions = find_seed_positions(query, reference, seed_size)

    if not seed_positions:
        seed_size = max(7, seed_size // 2)
        seed_positions = find_seed_positions(query, reference, seed_size)

    if not seed_positions:
        remaining = max(1, timeout - int(time.time() - start_time))
        return smith_waterman(query, reference, remaining)

    best_result = None
    best_score = -1

    band_extension = max(n // 2, 50)

    for seed_pos in seed_positions:
        elapsed = time.time() - start_time
        if elapsed > timeout:
            raise AlignmentTimeoutError(
                f"Alignment timed out after {timeout} seconds "
                f"({n} x {m} = {n*m:,} cells, seed-based optimization)"
            )

        region_start = max(0, seed_pos - band_extension)
        region_end = min(m, seed_pos + n + band_extension)
        region_len = region_end - region_start

        ref_subset = reference[region_start:region_end]
        remaining = max(1, timeout - int(time.time() - start_time))
        result = smith_waterman(query, ref_subset, remaining)

        if result["score"] > best_score:
            best_score = result["score"]
            best_result = dict(result)
            best_result["ref_start"] += region_start
            best_result["ref_end"] += region_start

    if best_result is None or best_score < MATCH_SCORE * 5:
        remaining = max(1, timeout - int(time.time() - start_time))
        return smith_waterman(query, reference, remaining)

    return best_result
