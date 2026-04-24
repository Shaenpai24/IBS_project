import argparse
import math
from collections import Counter, defaultdict
from pathlib import Path

import pandas as pd
from Bio import SeqIO, pairwise2
from Bio.Align import substitution_matrices
from Bio.SeqUtils.ProtParam import ProteinAnalysis


STANDARD_AA = "ACDEFGHIKLMNPQRSTVWY"

HYDROPATHY = {
    "A": 1.8,
    "C": 2.5,
    "D": -3.5,
    "E": -3.5,
    "F": 2.8,
    "G": -0.4,
    "H": -3.2,
    "I": 4.5,
    "K": -3.9,
    "L": 3.8,
    "M": 1.9,
    "N": -3.5,
    "P": -1.6,
    "Q": -3.5,
    "R": -4.5,
    "S": -0.8,
    "T": -0.7,
    "V": 4.2,
    "W": -0.9,
    "Y": -1.3,
}

VOLUME = {
    "A": 88.6,
    "C": 108.5,
    "D": 111.1,
    "E": 138.4,
    "F": 189.9,
    "G": 60.1,
    "H": 153.2,
    "I": 166.7,
    "K": 168.6,
    "L": 166.7,
    "M": 162.9,
    "N": 114.1,
    "P": 112.7,
    "Q": 143.8,
    "R": 173.4,
    "S": 89.0,
    "T": 116.1,
    "V": 140.0,
    "W": 227.8,
    "Y": 193.6,
}

CHARGE = {
    "A": 0,
    "C": 0,
    "D": -1,
    "E": -1,
    "F": 0,
    "G": 0,
    "H": 0,
    "I": 0,
    "K": 1,
    "L": 0,
    "M": 0,
    "N": 0,
    "P": 0,
    "Q": 0,
    "R": 1,
    "S": 0,
    "T": 0,
    "V": 0,
    "W": 0,
    "Y": 0,
}

POLAR = set("STNQYCHKRDE")

GRANTHAM = {
    ("A", "C"): 195, ("A", "D"): 126, ("A", "E"): 107, ("A", "F"): 113, ("A", "G"): 60,
    ("A", "H"): 86, ("A", "I"): 94, ("A", "K"): 106, ("A", "L"): 96, ("A", "M"): 84,
    ("A", "N"): 111, ("A", "P"): 27, ("A", "Q"): 91, ("A", "R"): 112, ("A", "S"): 99,
    ("A", "T"): 58, ("A", "V"): 64, ("A", "W"): 148, ("A", "Y"): 112,
    ("C", "D"): 154, ("C", "E"): 170, ("C", "F"): 205, ("C", "G"): 159, ("C", "H"): 174,
    ("C", "I"): 198, ("C", "K"): 202, ("C", "L"): 198, ("C", "M"): 196, ("C", "N"): 139,
    ("C", "P"): 169, ("C", "Q"): 154, ("C", "R"): 180, ("C", "S"): 112, ("C", "T"): 149,
    ("C", "V"): 192, ("C", "W"): 215, ("C", "Y"): 194,
    ("D", "E"): 45, ("D", "F"): 177, ("D", "G"): 94, ("D", "H"): 81, ("D", "I"): 168,
    ("D", "K"): 101, ("D", "L"): 172, ("D", "M"): 160, ("D", "N"): 23, ("D", "P"): 108,
    ("D", "Q"): 61, ("D", "R"): 96, ("D", "S"): 65, ("D", "T"): 85, ("D", "V"): 152,
    ("D", "W"): 181, ("D", "Y"): 160,
    ("E", "F"): 140, ("E", "G"): 98, ("E", "H"): 40, ("E", "I"): 134, ("E", "K"): 56,
    ("E", "L"): 138, ("E", "M"): 126, ("E", "N"): 42, ("E", "P"): 93, ("E", "Q"): 29,
    ("E", "R"): 54, ("E", "S"): 80, ("E", "T"): 65, ("E", "V"): 121, ("E", "W"): 152,
    ("E", "Y"): 122,
    ("F", "G"): 153, ("F", "H"): 100, ("F", "I"): 21, ("F", "K"): 102, ("F", "L"): 22,
    ("F", "M"): 28, ("F", "N"): 158, ("F", "P"): 114, ("F", "Q"): 116, ("F", "R"): 97,
    ("F", "S"): 155, ("F", "T"): 103, ("F", "V"): 50, ("F", "W"): 40, ("F", "Y"): 22,
    ("G", "H"): 98, ("G", "I"): 135, ("G", "K"): 127, ("G", "L"): 138, ("G", "M"): 127,
    ("G", "N"): 80, ("G", "P"): 42, ("G", "Q"): 87, ("G", "R"): 125, ("G", "S"): 56,
    ("G", "T"): 59, ("G", "V"): 109, ("G", "W"): 184, ("G", "Y"): 147,
    ("H", "I"): 94, ("H", "K"): 32, ("H", "L"): 99, ("H", "M"): 87, ("H", "N"): 68,
    ("H", "P"): 77, ("H", "Q"): 24, ("H", "R"): 29, ("H", "S"): 89, ("H", "T"): 47,
    ("H", "V"): 84, ("H", "W"): 115, ("H", "Y"): 83,
    ("I", "K"): 102, ("I", "L"): 5, ("I", "M"): 10, ("I", "N"): 149, ("I", "P"): 95,
    ("I", "Q"): 109, ("I", "R"): 97, ("I", "S"): 142, ("I", "T"): 89, ("I", "V"): 29,
    ("I", "W"): 61, ("I", "Y"): 33,
    ("K", "L"): 107, ("K", "M"): 95, ("K", "N"): 94, ("K", "P"): 103, ("K", "Q"): 53,
    ("K", "R"): 26, ("K", "S"): 121, ("K", "T"): 78, ("K", "V"): 97, ("K", "W"): 110,
    ("K", "Y"): 85,
    ("L", "M"): 15, ("L", "N"): 153, ("L", "P"): 98, ("L", "Q"): 113, ("L", "R"): 102,
    ("L", "S"): 145, ("L", "T"): 92, ("L", "V"): 32, ("L", "W"): 61, ("L", "Y"): 36,
    ("M", "N"): 142, ("M", "P"): 87, ("M", "Q"): 101, ("M", "R"): 91, ("M", "S"): 135,
    ("M", "T"): 81, ("M", "V"): 21, ("M", "W"): 67, ("M", "Y"): 36,
    ("N", "P"): 91, ("N", "Q"): 46, ("N", "R"): 86, ("N", "S"): 46, ("N", "T"): 65,
    ("N", "V"): 133, ("N", "W"): 174, ("N", "Y"): 143,
    ("P", "Q"): 76, ("P", "R"): 103, ("P", "S"): 74, ("P", "T"): 38, ("P", "V"): 68,
    ("P", "W"): 147, ("P", "Y"): 110,
    ("Q", "R"): 43, ("Q", "S"): 68, ("Q", "T"): 42, ("Q", "V"): 96, ("Q", "W"): 130,
    ("Q", "Y"): 99,
    ("R", "S"): 110, ("R", "T"): 71, ("R", "V"): 96, ("R", "W"): 101, ("R", "Y"): 77,
    ("S", "T"): 58, ("S", "V"): 124, ("S", "W"): 177, ("S", "Y"): 144,
    ("T", "V"): 69, ("T", "W"): 128, ("T", "Y"): 92,
    ("V", "W"): 88, ("V", "Y"): 55,
    ("W", "Y"): 37,
}


def pair_key(aa1, aa2):
    return tuple(sorted((aa1, aa2)))


def grantham_distance(aa1, aa2):
    if aa1 == aa2:
        return 0
    return GRANTHAM[pair_key(aa1, aa2)]


def aliphatic_index(seq):
    length = len(seq)
    if length == 0:
        return 0
    a = seq.count("A") / length
    v = seq.count("V") / length
    i = seq.count("I") / length
    l = seq.count("L") / length
    return (a * 100) + (2.9 * v * 100) + (3.9 * (i + l) * 100)


def net_charge(seq):
    pos = seq.count("K") + seq.count("R") + seq.count("H")
    neg = seq.count("D") + seq.count("E")
    return pos - neg


def motif_count(seq, motif="PHGGGWGQ"):
    return seq.count(motif)


def n_glyco_motif_count(seq):
    count = 0
    for i in range(len(seq) - 2):
        x1, x2, x3 = seq[i], seq[i + 1], seq[i + 2]
        if x1 == "N" and x2 != "P" and x3 in {"S", "T"}:
            count += 1
    return count


def hydrophobic_stretch(seq, window=5):
    hydrophobic = set("AILMFWYV")
    total = 0
    for i in range(len(seq) - window + 1):
        seg = seq[i : i + window]
        if all(aa in hydrophobic for aa in seg):
            total += 1
    return total


def compute_global_features(seq):
    analysed = ProteinAnalysis(seq)
    length = len(seq)

    return {
        "length": length,
        "molecular_weight": analysed.molecular_weight(),
        "aromaticity": analysed.aromaticity(),
        "instability_index": analysed.instability_index(),
        "gravy": analysed.gravy(),
        "pI": analysed.isoelectric_point(),
        "aliphatic_index": aliphatic_index(seq),
        "net_charge": net_charge(seq),
        "aromatic_density": (seq.count("F") + seq.count("W") + seq.count("Y")) / length,
        "glycine_density": seq.count("G") / length,
        "proline_count": seq.count("P"),
        "octapeptide_repeat_count": motif_count(seq),
        "hydrophobic_stretch_count": hydrophobic_stretch(seq, window=5),
        "n_glyco_motif_count": n_glyco_motif_count(seq),
    }


def translate_cds_records(cds_fasta, reference_seq, min_identity=0.6):
    proteins = []
    for record in SeqIO.parse(cds_fasta, "fasta"):
        seq = str(record.seq.translate(to_stop=True))
        clean = "".join(aa for aa in seq if aa in STANDARD_AA)
        if len(clean) < 200:
            continue

        aln = pairwise2.align.globalxx(reference_seq, clean, score_only=True)
        identity_to_ref = aln / len(reference_seq)
        if identity_to_ref < min_identity:
            continue

        proteins.append(clean)
    return proteins


def build_conservation_profile(reference_seq, homolog_proteins):
    counts_by_pos = defaultdict(Counter)

    for protein in homolog_proteins:
        alignment = pairwise2.align.globalms(reference_seq, protein, 2, -1, -5, -0.5, one_alignment_only=True)
        if not alignment:
            continue
        aligned_ref, aligned_seq = alignment[0].seqA, alignment[0].seqB

        ref_pos = 0
        for r_aa, s_aa in zip(aligned_ref, aligned_seq):
            if r_aa != "-":
                ref_pos += 1
                if s_aa != "-" and s_aa in STANDARD_AA:
                    counts_by_pos[ref_pos][s_aa] += 1

    profile = {}
    for pos in range(1, len(reference_seq) + 1):
        counter = counts_by_pos[pos]
        total = sum(counter.values())
        freqs = {}
        if total > 0:
            for aa, c in counter.items():
                freqs[aa] = c / total
            entropy = -sum(p * math.log2(p) for p in freqs.values())
        else:
            entropy = 0.0

        profile[pos] = {
            "counts": counter,
            "total": total,
            "freqs": freqs,
            "entropy": entropy,
        }

    return profile


def region_flags(position, seq_len):
    return {
        "is_signal_peptide_region": 1 if 1 <= position <= 22 else 0,
        "is_octapeptide_repeat_region": 1 if 51 <= position <= 91 else 0,
        "is_globular_domain_region": 1 if 125 <= position <= 228 else 0,
        "is_gpi_anchor_region": 1 if max(1, seq_len - 15) <= position <= seq_len else 0,
    }


def context_window(seq, position_1based, window=5):
    idx = position_1based - 1
    left = max(0, idx - window)
    right = min(len(seq), idx + window + 1)
    return seq[left:right]


def main():
    parser = argparse.ArgumentParser(description="Generate tier-2 PRNP mutation features.")
    parser.add_argument("--mutations_csv", default="prnp_single_point_mutations.csv")
    parser.add_argument("--reference_fasta", default="prnp_reference_protein.fasta")
    parser.add_argument("--source_cds_fasta", default="prnp_clean_nr.fasta")
    parser.add_argument("--output_csv", default="prnp_mutation_features_tier2.csv")
    args = parser.parse_args()

    if not Path(args.mutations_csv).exists():
        raise FileNotFoundError(f"Mutations CSV not found: {args.mutations_csv}")
    if not Path(args.reference_fasta).exists():
        raise FileNotFoundError(f"Reference FASTA not found: {args.reference_fasta}")
    if not Path(args.source_cds_fasta).exists():
        raise FileNotFoundError(f"Source CDS FASTA not found: {args.source_cds_fasta}")

    mutation_df = pd.read_csv(args.mutations_csv)
    reference_seq = str(next(SeqIO.parse(args.reference_fasta, "fasta")).seq)
    ref_len = len(reference_seq)

    homolog_proteins = translate_cds_records(args.source_cds_fasta, reference_seq, min_identity=0.6)
    conservation = build_conservation_profile(reference_seq, homolog_proteins)
    blosum62 = substitution_matrices.load("BLOSUM62")

    wt_global = compute_global_features(reference_seq)
    wt_octa = wt_global["octapeptide_repeat_count"]
    wt_nglyco = wt_global["n_glyco_motif_count"]

    rows = []
    for row in mutation_df.itertuples(index=False):
        pos = int(row.position)
        wt = row.original_aa
        mut = row.mutated_aa
        wt_seq = row.original_sequence
        mut_seq = row.mutated_sequence

        if wt_seq != reference_seq:
            raise ValueError("Found mutation row with original_sequence different from reference sequence.")

        if pos < 1 or pos > ref_len:
            raise ValueError(f"Mutation position out of reference range: {pos}")

        mut_global = compute_global_features(mut_seq)
        cons = conservation[pos]

        try:
            blosum_score = float(blosum62[(wt, mut)])
        except Exception:
            blosum_score = float("nan")

        wt_freq = cons["freqs"].get(wt, 0.0)
        mut_freq = cons["freqs"].get(mut, 0.0)

        out = {
            "mutation": row.mutation,
            "position": pos,
            "position_norm": pos / ref_len,
            "wt_aa": wt,
            "mut_aa": mut,
            "blosum62": blosum_score,
            "grantham": grantham_distance(wt, mut),
            "delta_hydropathy": HYDROPATHY[mut] - HYDROPATHY[wt],
            "delta_volume": VOLUME[mut] - VOLUME[wt],
            "delta_charge": CHARGE[mut] - CHARGE[wt],
            "wt_is_polar": 1 if wt in POLAR else 0,
            "mut_is_polar": 1 if mut in POLAR else 0,
            "delta_polarity_flag": (1 if mut in POLAR else 0) - (1 if wt in POLAR else 0),
            "local_window_wt": context_window(wt_seq, pos, window=5),
            "conservation_wt_freq": wt_freq,
            "conservation_mut_freq": mut_freq,
            "conservation_delta_freq": mut_freq - wt_freq,
            "position_entropy": cons["entropy"],
            "position_observation_count": cons["total"],
            "delta_instability_index": mut_global["instability_index"] - wt_global["instability_index"],
            "delta_gravy": mut_global["gravy"] - wt_global["gravy"],
            "delta_pI": mut_global["pI"] - wt_global["pI"],
            "delta_net_charge": mut_global["net_charge"] - wt_global["net_charge"],
            "delta_aromaticity": mut_global["aromaticity"] - wt_global["aromaticity"],
            "delta_aliphatic_index": mut_global["aliphatic_index"] - wt_global["aliphatic_index"],
            "delta_hydrophobic_stretch_count": mut_global["hydrophobic_stretch_count"] - wt_global["hydrophobic_stretch_count"],
            "delta_octapeptide_repeat_count": mut_global["octapeptide_repeat_count"] - wt_octa,
            "delta_n_glyco_motif_count": mut_global["n_glyco_motif_count"] - wt_nglyco,
        }

        out.update(region_flags(pos, ref_len))
        rows.append(out)

    out_df = pd.DataFrame(rows)
    out_df.to_csv(args.output_csv, index=False)

    print("Mutations loaded:", len(mutation_df))
    print("Reference length:", ref_len)
    print("Homolog proteins used for conservation:", len(homolog_proteins))
    print("Feature rows generated:", len(out_df))
    print("Feature columns:", out_df.shape[1])
    print("Saved:", args.output_csv)


if __name__ == "__main__":
    main()
