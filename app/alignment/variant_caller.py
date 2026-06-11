from typing import List, Dict


def extract_variants(alignment: Dict, ref_sequence: str) -> List[Dict]:
    """
    Extract variants (SNP, INS, DEL) from alignment result.
    Returns list of variant records with reference coordinates.
    """
    variants = []
    align_query = alignment["alignment_query"]
    align_ref = alignment["alignment_ref"]
    ref_start = alignment["ref_start"]

    ref_pos = ref_start
    query_pos = 0
    i = 0

    while i < len(align_query):
        q = align_query[i]
        r = align_ref[i]

        if q == "-":
            del_bases = []
            start_pos = ref_pos
            while i < len(align_query) and align_query[i] == "-":
                del_bases.append(align_ref[i])
                ref_pos += 1
                i += 1
            variants.append({
                "variant_type": "DEL",
                "ref_pos": start_pos,
                "ref_base": "".join(del_bases),
                "alt_base": "",
                "ref_pos_end": start_pos + len(del_bases) - 1,
            })
        elif r == "-":
            ins_bases = []
            start_pos = ref_pos
            while i < len(align_query) and align_ref[i] == "-":
                ins_bases.append(align_query[i])
                query_pos += 1
                i += 1
            variants.append({
                "variant_type": "INS",
                "ref_pos": start_pos,
                "ref_base": "",
                "alt_base": "".join(ins_bases),
                "ref_pos_end": start_pos,
            })
        else:
            if q != r:
                variants.append({
                    "variant_type": "SNP",
                    "ref_pos": ref_pos,
                    "ref_base": r,
                    "alt_base": q,
                    "ref_pos_end": ref_pos,
                })
            ref_pos += 1
            query_pos += 1
            i += 1

    return variants


def normalize_variant(variant: Dict, ref_sequence: str) -> Dict:
    """
    Normalize variant representation (left-align, trim common prefix/suffix).
    For insertions and deletions, include flanking base for VCF-style representation.
    """
    vtype = variant["variant_type"]
    ref_pos = variant["ref_pos"]

    if vtype == "SNP":
        return variant

    if vtype == "INS":
        if ref_pos > 0:
            flank_base = ref_sequence[ref_pos - 1]
            return {
                "variant_type": "INS",
                "ref_pos": ref_pos - 1,
                "ref_base": flank_base,
                "alt_base": flank_base + variant["alt_base"],
                "ref_pos_end": ref_pos - 1,
            }
        return variant

    if vtype == "DEL":
        return variant

    return variant
