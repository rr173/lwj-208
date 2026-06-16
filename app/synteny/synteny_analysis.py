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
