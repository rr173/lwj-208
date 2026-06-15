import math


R = 1.987

OLIGONUCLEOTIDE_CONCENTRATION = 250e-9

NN_DH = {
    "AA": -7.9, "TT": -7.9,
    "AT": -7.2,
    "TA": -7.2,
    "CA": -8.5, "TG": -8.5,
    "GT": -8.4, "AC": -8.4,
    "CT": -7.8, "AG": -7.8,
    "GA": -8.2, "TC": -8.2,
    "CG": -10.6,
    "GC": -9.8,
    "GG": -8.0, "CC": -8.0,
}

NN_DS = {
    "AA": -22.2, "TT": -22.2,
    "AT": -20.4,
    "TA": -21.3,
    "CA": -22.7, "TG": -22.7,
    "GT": -22.4, "AC": -22.4,
    "CT": -21.0, "AG": -21.0,
    "GA": -22.2, "TC": -22.2,
    "CG": -27.2,
    "GC": -24.4,
    "GG": -19.9, "CC": -19.9,
}

INIT_DH = {"G": 0.1, "C": 0.1, "A": 2.3, "T": 2.3}
INIT_DS = {"G": -2.8, "C": -2.8, "A": 4.1, "T": 4.1}


def compute_tm(seq: str, c: float = OLIGONUCLEOTIDE_CONCENTRATION) -> float:
    """
    Nearest-neighbor thermodynamic Tm calculation.
    Tm = ΔH / (ΔS + R * ln(C/4)) - 273.15
    ΔH in cal/mol, ΔS in cal/(mol·K), R = 1.987 cal/(mol·K).
    """
    seq = seq.upper()
    if len(seq) < 2:
        return 0.0

    dh = 0.0
    ds = 0.0

    dh += INIT_DH[seq[0]] * 1000.0
    ds += INIT_DS[seq[0]]
    dh += INIT_DH[seq[-1]] * 1000.0
    ds += INIT_DS[seq[-1]]

    for i in range(len(seq) - 1):
        dinuc = seq[i:i + 2]
        dh += NN_DH[dinuc] * 1000.0
        ds += NN_DS[dinuc]

    tm = dh / (ds + R * math.log(c / 4.0)) - 273.15
    return round(tm, 2)


def gc_content(seq: str) -> float:
    seq = seq.upper()
    if len(seq) == 0:
        return 0.0
    gc = sum(1 for b in seq if b in "GC")
    return round(gc / len(seq) * 100, 2)


def reverse_complement(seq: str) -> str:
    comp = {"A": "T", "T": "A", "G": "C", "C": "G"}
    return "".join(comp[b] for b in reversed(seq.upper()))
