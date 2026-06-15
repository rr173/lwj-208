import hashlib
import json
import math
from typing import List, Dict, Tuple, Optional

from app.alignment.smith_waterman import smith_waterman
from app.primer.thermodynamics import compute_tm, gc_content, reverse_complement


FLANK_SIZE = 200
MIN_PRIMER_LEN = 18
MAX_PRIMER_LEN = 25
MIN_GC = 40.0
MAX_GC = 60.0
MIN_TM = 55.0
MAX_TM = 65.0
OPTIMAL_TM = 60.0
MAX_MISMATCHES = 2
MAX_NONSPECIFIC_DIST = 3000
TOP_N = 10


def _cache_key(ref_name: str, target_start: int, target_end: int) -> str:
    raw = f"{ref_name}:{target_start}-{target_end}"
    return hashlib.sha256(raw.encode()).hexdigest()


def _find_binding_sites(primer: str, reference: str, max_mismatches: int = MAX_MISMATCHES) -> List[Dict]:
    """
    Use Smith-Waterman to find off-target binding sites for a primer on the full reference.
    Repeatedly run SW, mask the best hit region, and re-run until no significant hit found.
    A hit is considered valid if alignment covers >= (len(primer) - max_mismatches) bases
    with at most max_mismatches mismatches.
    """
    sites = []
    ref_masked = list(reference)
    primer_len = len(primer)
    min_cover = primer_len - max_mismatches

    for _ in range(50):
        ref_str = "".join(ref_masked)
        result = smith_waterman(primer, ref_str, timeout=30)

        if result["score"] == 0:
            break

        aln_len = len(result["alignment_query"])
        mismatches = result["alignment_match"].count("*")
        gaps = result["alignment_match"].count(" ")

        if aln_len < min_cover:
            break

        total_errors = mismatches + gaps
        if total_errors > max_mismatches:
            break

        ref_s = result["ref_start"]
        ref_e = result["ref_end"]
        sites.append({
            "ref_start": ref_s,
            "ref_end": ref_e,
            "mismatches": mismatches,
            "gaps": gaps,
            "score": result["score"],
            "alignment_ref": result["alignment_ref"],
        })

        for i in range(max(0, ref_s - 2), min(len(ref_masked), ref_e + 3)):
            ref_masked[i] = "N"

    return sites


def _check_nonspecific_amplification(
    fwd_sites: List[Dict],
    rev_sites: List[Dict],
    target_start: int,
    target_end: int,
    rev_len: int,
) -> Tuple[bool, List[Dict]]:
    """
    Check if forward and reverse primers have off-target binding sites
    in opposite directions with distance < 3000bp, producing non-specific amplicons.
    """
    nonspecific_pairs = []

    fwd_off = [s for s in fwd_sites if not (s["ref_start"] >= target_start and s["ref_end"] <= target_end)]
    rev_off = [s for s in rev_sites if not (s["ref_start"] >= target_start and s["ref_end"] <= target_end)]

    for fs in fwd_off:
        for rs in rev_off:
            amp_start = fs["ref_start"]
            amp_end = rs["ref_end"]
            if amp_start < amp_end:
                amp_len = amp_end - amp_start + 1
                if amp_len < MAX_NONSPECIFIC_DIST:
                    in_target = (fs["ref_start"] >= target_start and fs["ref_end"] <= target_end and
                                 rs["ref_start"] >= target_start and rs["ref_end"] <= target_end)
                    if not in_target:
                        nonspecific_pairs.append({
                            "fwd_site": fs,
                            "rev_site": rs,
                            "amplicon_start": amp_start,
                            "amplicon_end": amp_end,
                            "amplicon_length": amp_len,
                        })

    return len(nonspecific_pairs) > 0, nonspecific_pairs


def _generate_primer_candidates(
    seq_region: str,
    region_offset: int,
    is_reverse: bool = False,
) -> List[Dict]:
    """
    Generate primer candidates from a sequence region.
    For forward primers: take subsequences directly.
    For reverse primers: take reverse complement of subsequences.
    """
    candidates = []
    for length in range(MIN_PRIMER_LEN, MAX_PRIMER_LEN + 1):
        for i in range(len(seq_region) - length + 1):
            subseq = seq_region[i:i + length]

            if is_reverse:
                primer_seq = reverse_complement(subseq)
                bind_start = region_offset + i
                bind_end = region_offset + i + length - 1
            else:
                primer_seq = subseq
                bind_start = region_offset + i
                bind_end = region_offset + i + length - 1

            gc = gc_content(primer_seq)
            if gc < MIN_GC or gc > MAX_GC:
                continue

            tm = compute_tm(primer_seq)
            if tm < MIN_TM or tm > MAX_TM:
                continue

            candidates.append({
                "sequence": primer_seq,
                "length": length,
                "gc": gc,
                "tm": tm,
                "tm_deviation": abs(tm - OPTIMAL_TM),
                "bind_start": bind_start,
                "bind_end": bind_end,
                "is_reverse": is_reverse,
            })

    return candidates


def design_primers(
    reference_sequence: str,
    ref_name: str,
    target_start: int,
    target_end: int,
) -> Dict:
    """
    Design primer pairs for a target region.
    Returns top 10 primer pairs sorted by Tm deviation from 60°C.
    """
    ref_len = len(reference_sequence)

    if target_start < 0 or target_end >= ref_len or target_start >= target_end:
        return {"error": "Invalid target coordinates", "primer_pairs": []}

    left_start = max(0, target_start - FLANK_SIZE)
    left_end = target_start
    right_start = target_end + 1
    right_end = min(ref_len, target_end + 1 + FLANK_SIZE)

    left_region = reference_sequence[left_start:left_end]
    right_region = reference_sequence[right_start:right_end]

    fwd_candidates = _generate_primer_candidates(left_region, left_start, is_reverse=False)
    rev_candidates = _generate_primer_candidates(right_region, right_start, is_reverse=True)

    if not fwd_candidates or not rev_candidates:
        return {"error": "No valid primer candidates found", "primer_pairs": []}

    fwd_candidates.sort(key=lambda x: (x["tm_deviation"], x["length"]))
    rev_candidates.sort(key=lambda x: (x["tm_deviation"], x["length"]))

    top_fwd = fwd_candidates[:30]
    top_rev = rev_candidates[:30]

    pairs = []
    for fwd in top_fwd:
        for rev in top_rev:
            avg_tm_dev = (fwd["tm_deviation"] + rev["tm_deviation"]) / 2.0
            tm_diff = abs(fwd["tm"] - rev["tm"])
            amplicon_len = rev["bind_end"] - fwd["bind_start"] + 1

            pairs.append({
                "forward": fwd,
                "reverse": rev,
                "avg_tm_deviation": round(avg_tm_dev, 2),
                "tm_difference": round(tm_diff, 2),
                "amplicon_length": amplicon_len,
            })

    pairs.sort(key=lambda x: (x["avg_tm_deviation"], x["tm_difference"]))

    top_pairs = pairs[:TOP_N]

    results = []
    for idx, pair in enumerate(top_pairs):
        fwd = pair["forward"]
        rev = pair["reverse"]

        fwd_sites = _find_binding_sites(fwd["sequence"], reference_sequence, MAX_MISMATCHES)
        rev_sites = _find_binding_sites(rev["sequence"], reference_sequence, MAX_MISMATCHES)

        has_nonspecific, nonspecific_details = _check_nonspecific_amplification(
            fwd_sites, rev_sites,
            target_start, target_end,
            rev["length"],
        )

        fwd_on_target = [
            s for s in fwd_sites
            if s["ref_start"] >= target_start - FLANK_SIZE and s["ref_end"] <= target_start + 10
        ]
        rev_on_target = [
            s for s in rev_sites
            if s["ref_start"] >= target_end - 10 and s["ref_end"] <= target_end + FLANK_SIZE
        ]

        results.append({
            "rank": idx + 1,
            "forward_primer": {
                "sequence": fwd["sequence"],
                "length": fwd["length"],
                "gc_percent": fwd["gc"],
                "tm": fwd["tm"],
                "bind_start": fwd["bind_start"],
                "bind_end": fwd["bind_end"],
                "off_target_sites": [
                    {
                        "ref_start": s["ref_start"],
                        "ref_end": s["ref_end"],
                        "mismatches": s["mismatches"],
                    }
                    for s in fwd_sites
                    if s not in fwd_on_target
                ],
            },
            "reverse_primer": {
                "sequence": rev["sequence"],
                "length": rev["length"],
                "gc_percent": rev["gc"],
                "tm": rev["tm"],
                "bind_start": rev["bind_start"],
                "bind_end": rev["bind_end"],
                "off_target_sites": [
                    {
                        "ref_start": s["ref_start"],
                        "ref_end": s["ref_end"],
                        "mismatches": s["mismatches"],
                    }
                    for s in rev_sites
                    if s not in rev_on_target
                ],
            },
            "amplicon_length": pair["amplicon_length"],
            "avg_tm_deviation": pair["avg_tm_deviation"],
            "tm_difference": pair["tm_difference"],
            "has_nonspecific_risk": has_nonspecific,
            "nonspecific_amplicons": nonspecific_details if has_nonspecific else [],
        })

    return {
        "reference_name": ref_name,
        "target_start": target_start,
        "target_end": target_end,
        "flank_size": FLANK_SIZE,
        "primer_pairs": results,
    }


def get_cached_primer_result(db, ref_name: str, target_start: int, target_end: int) -> Optional[dict]:
    from app.models import PrimerDesignCache
    key = _cache_key(ref_name, target_start, target_end)
    cached = db.query(PrimerDesignCache).filter(
        PrimerDesignCache.cache_key == key
    ).first()
    if cached:
        return json.loads(cached.result_json)
    return None


def save_primer_result(db, ref_name: str, target_start: int, target_end: int, result: dict) -> None:
    from app.models import PrimerDesignCache
    key = _cache_key(ref_name, target_start, target_end)
    existing = db.query(PrimerDesignCache).filter(
        PrimerDesignCache.cache_key == key
    ).first()
    result_json = json.dumps(result, ensure_ascii=False)
    if existing:
        existing.result_json = result_json
    else:
        entry = PrimerDesignCache(
            cache_key=key,
            reference_name=ref_name,
            target_start=target_start,
            target_end=target_end,
            result_json=result_json,
        )
        db.add(entry)
    db.commit()


def _std(values: List[float]) -> float:
    if len(values) < 2:
        return 0.0
    mean = sum(values) / len(values)
    variance = sum((v - mean) ** 2 for v in values) / len(values)
    return round(math.sqrt(variance), 4)


def compute_tm_consistency(best_pairs: List[Dict]) -> Optional[Dict]:
    valid_pairs = [p for p in best_pairs if p is not None]
    if not valid_pairs:
        return None

    all_fwd_tm = [p["forward_primer"]["tm"] for p in valid_pairs]
    all_rev_tm = [p["reverse_primer"]["tm"] for p in valid_pairs]
    combined_tm = all_fwd_tm + all_rev_tm

    return {
        "all_fwd_tm": all_fwd_tm,
        "all_rev_tm": all_rev_tm,
        "fwd_tm_std": _std(all_fwd_tm),
        "rev_tm_std": _std(all_rev_tm),
        "combined_tm_std": _std(combined_tm),
    }


def get_primer_design_history(db, reference_name: str) -> List[Dict]:
    from app.models import PrimerDesignCache
    records = (
        db.query(PrimerDesignCache)
        .filter(PrimerDesignCache.reference_name == reference_name)
        .order_by(PrimerDesignCache.created_at.desc())
        .all()
    )

    result = []
    for record in records:
        result_data = json.loads(record.result_json)
        primer_pairs = result_data.get("primer_pairs", [])
        best_pair = primer_pairs[0] if primer_pairs else None
        result.append({
            "id": record.id,
            "reference_name": record.reference_name,
            "target_start": record.target_start,
            "target_end": record.target_end,
            "best_primer_pair": best_pair,
            "created_at": record.created_at,
        })
    return result
