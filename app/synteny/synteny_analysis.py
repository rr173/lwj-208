import time
import hashlib
from dataclasses import dataclass, field
from typing import List, Tuple, Optional, Dict
from app.alignment.smith_waterman import optimized_smith_waterman, AlignmentTimeoutError


ANCHOR_LENGTH = 50
SCORE_THRESHOLD_RATIO = 1.5
SYNTENY_TIMEOUT_SECONDS = 60
MAX_GAP_RATIO = 3.0


@dataclass
class SyntenyAnchor:
    a_start: int
    a_end: int
    b_start: int
    b_end: int
    score: int
    direction: str = "+"


@dataclass
class SyntenyBlock:
    a_start: int
    a_end: int
    b_start: int
    b_end: int
    anchor_count: int
    direction: str
    anchors: List[SyntenyAnchor] = field(default_factory=list)


@dataclass
class RearrangementEvent:
    event_type: str
    a_start: int
    a_end: int
    b_start: int
    b_end: int
    anchor_count: int


_COMP = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}


def _reverse_complement(seq: str) -> str:
    return ''.join(_COMP.get(c, 'N') for c in reversed(seq))


def _generate_anchors(seq_a: str, anchor_length: int = ANCHOR_LENGTH) -> List[Tuple[int, str]]:
    anchors = []
    n = len(seq_a)
    for i in range(0, n - anchor_length + 1, anchor_length):
        anchors.append((i, seq_a[i:i + anchor_length]))
    return anchors


def _align_anchors(
    anchors: List[Tuple[int, str]],
    seq_b: str,
    score_threshold_ratio: float = SCORE_THRESHOLD_RATIO,
    timeout: int = SYNTENY_TIMEOUT_SECONDS,
) -> List[SyntenyAnchor]:
    aligned_anchors = []
    start_time = time.time()

    for a_start, anchor_seq in anchors:
        if time.time() - start_time > timeout:
            raise AlignmentTimeoutError(
                f"Synteny analysis timed out after {timeout} seconds"
            )

        anchor_len = len(anchor_seq)
        threshold = anchor_len * score_threshold_ratio

        rc_anchor = _reverse_complement(anchor_seq)

        forward_result = optimized_smith_waterman(anchor_seq, seq_b)
        rc_result = optimized_smith_waterman(rc_anchor, seq_b)

        best_result = forward_result
        direction = "+"

        if rc_result["score"] > forward_result["score"]:
            best_result = rc_result
            direction = "-"

        if best_result["score"] >= threshold and best_result["score"] > 0:
            aligned_anchors.append(SyntenyAnchor(
                a_start=a_start,
                a_end=a_start + anchor_len - 1,
                b_start=best_result["ref_start"],
                b_end=best_result["ref_end"],
                score=best_result["score"],
                direction=direction,
            ))

    return aligned_anchors


def find_synteny_blocks(
    seq_a: str,
    seq_b: str,
    anchor_length: int = ANCHOR_LENGTH,
    score_threshold_ratio: float = SCORE_THRESHOLD_RATIO,
    timeout: int = SYNTENY_TIMEOUT_SECONDS,
    max_gap_ratio: float = MAX_GAP_RATIO,
) -> List[SyntenyBlock]:
    start_time = time.time()

    anchors = _generate_anchors(seq_a, anchor_length)

    remaining_time = max(1, timeout - int(time.time() - start_time))
    aligned_anchors = _align_anchors(anchors, seq_b, score_threshold_ratio, remaining_time)

    if not aligned_anchors:
        return []

    aligned_anchors.sort(key=lambda x: x.a_start)

    max_gap = anchor_length * max_gap_ratio

    blocks = []
    current_anchors = [aligned_anchors[0]]
    current_direction = aligned_anchors[0].direction

    for anchor in aligned_anchors[1:]:
        if time.time() - start_time > timeout:
            raise AlignmentTimeoutError(
                f"Synteny analysis timed out after {timeout} seconds"
            )

        is_monotonic = False
        if current_direction == "+":
            is_monotonic = anchor.b_start >= current_anchors[-1].b_start
        else:
            is_monotonic = anchor.b_start <= current_anchors[-1].b_start

        same_direction = anchor.direction == current_direction

        gap_on_b = abs(anchor.b_start - current_anchors[-1].b_end)
        within_gap_limit = gap_on_b <= max_gap

        if is_monotonic and same_direction and within_gap_limit:
            current_anchors.append(anchor)
        else:
            block = _build_block(current_anchors, current_direction)
            blocks.append(block)
            current_anchors = [anchor]
            current_direction = anchor.direction

    block = _build_block(current_anchors, current_direction)
    blocks.append(block)

    return blocks


def _build_block(anchors: List[SyntenyAnchor], direction: str) -> SyntenyBlock:
    a_start = anchors[0].a_start
    a_end = anchors[-1].a_end

    if direction == "+":
        b_start = anchors[0].b_start
        b_end = anchors[-1].b_end
    else:
        b_start = anchors[0].b_end
        b_end = anchors[-1].b_start

    return SyntenyBlock(
        a_start=a_start,
        a_end=a_end,
        b_start=min(b_start, b_end),
        b_end=max(b_start, b_end),
        anchor_count=len(anchors),
        direction=direction,
        anchors=anchors,
    )


def detect_rearrangements(
    blocks: List[SyntenyBlock],
    seq_a_len: int,
    seq_b_len: int,
) -> List[RearrangementEvent]:
    events = []

    for block in blocks:
        if block.direction == "-":
            events.append(RearrangementEvent(
                event_type="inversion",
                a_start=block.a_start,
                a_end=block.a_end,
                b_start=block.b_start,
                b_end=block.b_end,
                anchor_count=block.anchor_count,
            ))

    if len(blocks) < 2:
        _detect_duplications(blocks, events)
        return events

    for i in range(len(blocks) - 1):
        curr_block = blocks[i]
        next_block = blocks[i + 1]

        a_gap_start = curr_block.a_end + 1
        a_gap_end = next_block.a_start - 1

        if a_gap_start > a_gap_end:
            continue

        if curr_block.direction != next_block.direction:
            events.append(RearrangementEvent(
                event_type="inversion",
                a_start=a_gap_start,
                a_end=a_gap_end,
                b_start=min(curr_block.b_end, next_block.b_start),
                b_end=max(curr_block.b_end, next_block.b_start),
                anchor_count=0,
            ))

    for i in range(len(blocks) - 1):
        curr = blocks[i]
        nxt = blocks[i + 1]
        if curr.direction == nxt.direction and curr.b_start > nxt.b_start:
            a_gap_start = curr.a_end + 1
            a_gap_end = nxt.a_start - 1

            if a_gap_start > a_gap_end:
                a_gap_start = min(curr.a_start, nxt.a_start)
                a_gap_end = max(curr.a_end, nxt.a_end)

            b_gap_start = min(curr.b_end, nxt.b_start)
            b_gap_end = max(curr.b_end, nxt.b_start)

            events.append(RearrangementEvent(
                event_type="translocation",
                a_start=a_gap_start,
                a_end=a_gap_end,
                b_start=b_gap_start,
                b_end=b_gap_end,
                anchor_count=0,
            ))
            break

    _detect_duplications(blocks, events)

    return events


def _detect_duplications(blocks: List[SyntenyBlock], events: List[RearrangementEvent]):
    a_intervals = []
    for block in blocks:
        a_intervals.append((block.a_start, block.a_end, block))

    for i in range(len(a_intervals)):
        for j in range(i + 1, len(a_intervals)):
            a_start_i, a_end_i, block_i = a_intervals[i]
            a_start_j, a_end_j, block_j = a_intervals[j]

            overlap_start = max(a_start_i, a_start_j)
            overlap_end = min(a_end_i, a_end_j)

            if overlap_start <= overlap_end:
                b_distinct = (
                    block_i.b_end < block_j.b_start
                    or block_j.b_end < block_i.b_start
                )

                if b_distinct:
                    events.append(RearrangementEvent(
                        event_type="duplication",
                        a_start=overlap_start,
                        a_end=overlap_end,
                        b_start=min(block_i.b_start, block_j.b_start),
                        b_end=max(block_i.b_end, block_j.b_end),
                        anchor_count=block_i.anchor_count + block_j.anchor_count,
                    ))


def analyze_synteny(
    seq_a: str,
    seq_b: str,
    anchor_length: int = ANCHOR_LENGTH,
    score_threshold_ratio: float = SCORE_THRESHOLD_RATIO,
    timeout: int = SYNTENY_TIMEOUT_SECONDS,
    max_gap_ratio: float = MAX_GAP_RATIO,
) -> Dict:
    start_time = time.time()

    blocks = find_synteny_blocks(
        seq_a, seq_b, anchor_length, score_threshold_ratio, timeout, max_gap_ratio
    )

    rearrangements = detect_rearrangements(blocks, len(seq_a), len(seq_b))

    b_non_monotonic_transitions = []
    for i in range(len(blocks) - 1):
        curr = blocks[i]
        nxt = blocks[i + 1]
        if nxt.b_start < curr.b_end:
            b_non_monotonic_transitions.append(i)

    return {
        "seq_a_length": len(seq_a),
        "seq_b_length": len(seq_b),
        "anchor_length": anchor_length,
        "score_threshold_ratio": score_threshold_ratio,
        "max_gap_ratio": max_gap_ratio,
        "synteny_blocks": blocks,
        "rearrangements": rearrangements,
        "b_non_monotonic_transitions": b_non_monotonic_transitions,
        "total_anchors_aligned": sum(b.anchor_count for b in blocks),
    }


def sequence_pair_hash(seq_a_name: str, seq_b_name: str) -> str:
    sorted_names = sorted([seq_a_name, seq_b_name])
    combined = f"{sorted_names[0]}__{sorted_names[1]}"
    return hashlib.md5(combined.encode()).hexdigest()


@dataclass
class PairwiseSyntenyResult:
    ref_a_name: str
    ref_b_name: str
    blocks: List[SyntenyBlock]
    total_anchors: int
    cached: bool = False


@dataclass
class ConservedCoreInterval:
    start: int
    end: int
    anchor_scores: List[float]
    other_intervals: Dict[str, Tuple[int, int]]


def _interval_overlap(start1: int, end1: int, start2: int, end2: int) -> Tuple[int, int]:
    overlap_start = max(start1, start2)
    overlap_end = min(end1, end2)
    if overlap_start > overlap_end:
        return (-1, -1)
    return (overlap_start, overlap_end)


def _find_blocks_covering_position(blocks: List[SyntenyBlock], pos: int) -> List[SyntenyBlock]:
    covering = []
    for block in blocks:
        if block.a_start <= pos <= block.a_end:
            covering.append(block)
    return covering


def _get_b_interval_for_a_position(block: SyntenyBlock, a_pos: int) -> Optional[int]:
    if block.a_start > a_pos or block.a_end < a_pos:
        return None

    relative_pos = (a_pos - block.a_start) / max(1, block.a_end - block.a_start)

    if block.direction == "+":
        b_pos = int(block.b_start + relative_pos * (block.b_end - block.b_start))
    else:
        b_pos = int(block.b_end - relative_pos * (block.b_end - block.b_start))

    return max(block.b_start, min(block.b_end, b_pos))


def _get_anchor_scores_in_interval(block: SyntenyBlock, start: int, end: int) -> List[int]:
    scores = []
    for anchor in block.anchors:
        if anchor.a_end < start or anchor.a_start > end:
            continue
        overlap_start = max(anchor.a_start, start)
        overlap_end = min(anchor.a_end, end)
        if overlap_start <= overlap_end:
            overlap_len = overlap_end - overlap_start + 1
            anchor_len = anchor.a_end - anchor.a_start + 1
            weighted_score = anchor.score * (overlap_len / anchor_len)
            scores.append(int(weighted_score))
    return scores


def infer_conserved_core_regions(
    pairwise_results: List[PairwiseSyntenyResult],
    reference_names: List[str],
    reference_lengths: Dict[str, int],
    anchor_length: int,
) -> Dict:
    if len(reference_names) < 3:
        raise ValueError("At least 3 reference sequences are required for multi-sequence synteny")

    primary_name = reference_names[0]
    primary_length = reference_lengths.get(primary_name, 0)
    if primary_length == 0:
        raise ValueError(f"Reference sequence '{primary_name}' has zero length")

    pairwise_by_other = {}
    for pr in pairwise_results:
        other_name = pr.ref_b_name if pr.ref_a_name == primary_name else pr.ref_a_name
        pairwise_by_other[other_name] = pr

    coverage = [0] * primary_length
    position_scores = [[] for _ in range(primary_length)]
    position_other_intervals: List[Dict[str, Tuple[int, int]]] = [{} for _ in range(primary_length)]

    for pos in range(primary_length):
        all_covered = True
        for other_name in reference_names[1:]:
            pr = pairwise_by_other.get(other_name)
            if not pr:
                all_covered = False
                break

            blocks = pr.blocks
            covering_blocks = _find_blocks_covering_position(blocks, pos)
            if not covering_blocks:
                all_covered = False
                break

            best_block = max(covering_blocks, key=lambda b: b.anchor_count)
            b_pos = _get_b_interval_for_a_position(best_block, pos)
            if b_pos is None:
                all_covered = False
                break

            if primary_name == pr.ref_a_name:
                other_start, other_end = best_block.b_start, best_block.b_end
            else:
                other_start, other_end = best_block.a_start, best_block.a_end

            position_other_intervals[pos][other_name] = (other_start, other_end)

            scores = _get_anchor_scores_in_interval(best_block, pos, pos + anchor_length - 1)
            if scores:
                position_scores[pos].extend(scores)

        if all_covered:
            coverage[pos] = 1

    core_regions: List[ConservedCoreInterval] = []
    current_start = -1

    for pos in range(primary_length):
        if coverage[pos] == 1 and current_start == -1:
            current_start = pos
        elif coverage[pos] == 0 and current_start != -1:
            core_end = pos - 1
            if core_end >= current_start:
                all_scores = []
                other_intervals = {}
                for other_name in reference_names[1:]:
                    intervals = set()
                    for p in range(current_start, core_end + 1):
                        if other_name in position_other_intervals[p]:
                            intervals.add(position_other_intervals[p][other_name])
                    if intervals:
                        all_other_starts = [i[0] for i in intervals]
                        all_other_ends = [i[1] for i in intervals]
                        other_intervals[other_name] = (min(all_other_starts), max(all_other_ends))
                for p in range(current_start, core_end + 1):
                    all_scores.extend(position_scores[p])
                avg_score = sum(all_scores) / len(all_scores) if all_scores else 0.0
                core_regions.append(ConservedCoreInterval(
                    start=current_start,
                    end=core_end,
                    anchor_scores=all_scores,
                    other_intervals=other_intervals,
                ))
            current_start = -1

    if current_start != -1:
        core_end = primary_length - 1
        all_scores = []
        other_intervals = {}
        for other_name in reference_names[1:]:
            intervals = set()
            for p in range(current_start, core_end + 1):
                if other_name in position_other_intervals[p]:
                    intervals.add(position_other_intervals[p][other_name])
            if intervals:
                all_other_starts = [i[0] for i in intervals]
                all_other_ends = [i[1] for i in intervals]
                other_intervals[other_name] = (min(all_other_starts), max(all_other_ends))
        for p in range(current_start, core_end + 1):
            all_scores.extend(position_scores[p])
        avg_score = sum(all_scores) / len(all_scores) if all_scores else 0.0
        core_regions.append(ConservedCoreInterval(
            start=current_start,
            end=core_end,
            anchor_scores=all_scores,
            other_intervals=other_intervals,
        ))

    total_conserved = sum(c.end - c.start + 1 for c in core_regions)
    conservation_ratio = total_conserved / primary_length if primary_length > 0 else 0.0

    result_regions = []
    for idx, core in enumerate(core_regions):
        intervals = [{
            "sequence_name": primary_name,
            "start": core.start,
            "end": core.end,
        }]
        for other_name, (s, e) in core.other_intervals.items():
            intervals.append({
                "sequence_name": other_name,
                "start": s,
                "end": e,
            })

        avg_score = sum(core.anchor_scores) / len(core.anchor_scores) if core.anchor_scores else 0.0

        result_regions.append({
            "region_index": idx,
            "intervals": intervals,
            "avg_anchor_score": round(avg_score, 2),
            "sequence_count": len(reference_names),
            "total_length": core.end - core.start + 1,
        })

    return {
        "conserved_core_regions": result_regions,
        "total_conserved_length": total_conserved,
        "conservation_ratio": round(conservation_ratio, 6),
    }
