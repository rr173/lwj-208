from typing import List, Dict, Optional, Tuple


CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

STOP_CODONS = {"TAA", "TAG", "TGA"}


def translate_codon(codon: str) -> str:
    """Translate a codon to amino acid. Returns 'X' for invalid codons."""
    codon = codon.upper()
    if len(codon) != 3:
        return "X"
    return CODON_TABLE.get(codon, "X")


def find_feature_at_position(annotations: List[Dict], pos: int) -> Optional[Dict]:
    """Find the gene feature that contains the given position."""
    for ann in annotations:
        if ann["start"] <= pos <= ann["end"]:
            return ann
    return None


def get_codon_position(exon_start: int, pos: int, strand: str = "+") -> int:
    """
    Calculate the position within the codon (0, 1, or 2)
    given the exon start and the position.
    Assumes coding starts at exon_start (0-based or 1-based? need to be consistent).
    Using 0-based positions.
    """
    if strand == "+":
        offset = pos - exon_start
    else:
        offset = exon_start - pos
    return offset % 3


def annotate_variant(variant: Dict, annotations: List[Dict], ref_sequence: str) -> Dict:
    """
    Annotate a single variant with gene feature and impact.
    variant format: {variant_type, ref_pos, ref_base, alt_base, ref_pos_end}
    """
    result = dict(variant)
    result["feature_type"] = "intergenic"
    result["gene_name"] = None
    result["impact"] = "MODIFIER"
    result["consequence"] = None
    result["codon_ref"] = None
    result["codon_alt"] = None
    result["aa_ref"] = None
    result["aa_alt"] = None

    ref_pos = variant["ref_pos"]
    vtype = variant["variant_type"]

    feature = find_feature_at_position(annotations, ref_pos)

    if feature is None:
        result["feature_type"] = "intergenic"
        result["impact"] = "MODIFIER"
        return result

    result["feature_type"] = feature["feature_type"]
    result["gene_name"] = feature["gene_name"]

    ftype = feature["feature_type"]

    if ftype == "intron":
        result["impact"] = "LOW"
        result["consequence"] = "intron_variant"
        return result

    if ftype == "UTR":
        result["impact"] = "LOW"
        result["consequence"] = "UTR_variant"
        return result

    if ftype == "intergenic":
        result["impact"] = "MODIFIER"
        result["consequence"] = "intergenic_variant"
        return result

    if ftype != "exon":
        result["impact"] = "MODIFIER"
        return result

    return _annotate_exonic_variant(result, variant, feature, ref_sequence)


def _annotate_exonic_variant(result: Dict, variant: Dict, feature: Dict, ref_sequence: str) -> Dict:
    """Annotate a variant that falls within an exon."""
    vtype = variant["variant_type"]
    ref_pos = variant["ref_pos"]
    exon_start = feature["start"]
    strand = feature.get("strand", "+")

    result["consequence"] = "exon_variant"

    if vtype == "SNP":
        codon_offset = (ref_pos - exon_start) % 3
        codon_start = ref_pos - codon_offset

        if codon_start + 2 < len(ref_sequence):
            ref_codon = ref_sequence[codon_start:codon_start + 3].upper()
            ref_aa = translate_codon(ref_codon)

            alt_codon_list = list(ref_codon)
            alt_codon_list[codon_offset] = variant["alt_base"].upper()
            alt_codon = "".join(alt_codon_list)
            alt_aa = translate_codon(alt_codon)

            result["codon_ref"] = ref_codon
            result["codon_alt"] = alt_codon
            result["aa_ref"] = ref_aa
            result["aa_alt"] = alt_aa

            if ref_aa == alt_aa:
                result["consequence"] = "synonymous_variant"
                result["impact"] = "LOW"
            elif alt_aa == "*":
                result["consequence"] = "stop_gained"
                result["impact"] = "HIGH"
            elif ref_aa == "*":
                result["consequence"] = "stop_lost"
                result["impact"] = "HIGH"
            else:
                result["consequence"] = "missense_variant"
                result["impact"] = "MODERATE"
        else:
            result["impact"] = "MODERATE"
            result["consequence"] = "missense_variant"

    elif vtype == "INS":
        alt_bases = variant["alt_base"]
        if len(alt_bases) % 3 != 0:
            result["consequence"] = "frameshift_variant"
            result["impact"] = "HIGH"
        else:
            result["consequence"] = "inframe_insertion"
            result["impact"] = "MODERATE"

    elif vtype == "DEL":
        ref_bases = variant["ref_base"]
        if len(ref_bases) % 3 != 0:
            result["consequence"] = "frameshift_variant"
            result["impact"] = "HIGH"
        else:
            result["consequence"] = "inframe_deletion"
            result["impact"] = "MODERATE"

    return result


def annotate_variants(variants: List[Dict], annotations: List[Dict], ref_sequence: str) -> List[Dict]:
    """Annotate a list of variants."""
    return [annotate_variant(v, annotations, ref_sequence) for v in variants]
