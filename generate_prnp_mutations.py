import argparse
from pathlib import Path

import pandas as pd
from Bio import SeqIO


STANDARD_AA = "ACDEFGHIKLMNPQRSTVWY"


def select_reference_record(records):
    ranked = []

    for record in records:
        header = record.description
        score = 0

        if "[gene=PRNP]" in header:
            score += 3
        if "NM_001080121.3" in header or "NP_001073590.1" in header:
            score += 100
        if "GeneID:5621" in header:
            score += 10
        if "Homo sapiens" in header:
            score += 5

        ranked.append((score, record))

    ranked.sort(key=lambda x: x[0], reverse=True)
    best_score, best_record = ranked[0]

    if best_score <= 0:
        raise ValueError("Could not find a PRNP reference sequence in the FASTA file.")

    return best_record


def cds_to_clean_protein(record):
    protein = str(record.seq.translate(to_stop=True))
    clean = "".join(aa for aa in protein if aa in STANDARD_AA)

    if not clean:
        raise ValueError("Selected reference translated to an empty/invalid protein sequence.")

    return clean


def generate_single_point_mutations(reference_seq):
    rows = []

    for index_1based, original_aa in enumerate(reference_seq, start=1):
        for mutated_aa in STANDARD_AA:
            if mutated_aa == original_aa:
                continue

            mutated_seq_list = list(reference_seq)
            mutated_seq_list[index_1based - 1] = mutated_aa
            mutated_seq = "".join(mutated_seq_list)

            rows.append(
                {
                    "position": index_1based,
                    "original_aa": original_aa,
                    "mutated_aa": mutated_aa,
                    "mutation": f"{original_aa}{index_1based}{mutated_aa}",
                    "original_sequence": reference_seq,
                    "mutated_sequence": mutated_seq,
                }
            )

    return pd.DataFrame(rows)


def main():
    parser = argparse.ArgumentParser(
        description="Generate all possible single-point amino-acid substitutions for a canonical PRNP sequence."
    )
    parser.add_argument(
        "--input_fasta",
        default="prnp_clean_nr.fasta",
        help="Input FASTA file with cleaned PRNP CDS sequences.",
    )
    parser.add_argument(
        "--output_csv",
        default="prnp_single_point_mutations.csv",
        help="Output CSV path for generated mutations.",
    )
    parser.add_argument(
        "--reference_fasta",
        default="prnp_reference_protein.fasta",
        help="Output FASTA path for the chosen translated reference protein.",
    )
    args = parser.parse_args()

    input_path = Path(args.input_fasta)
    if not input_path.exists():
        raise FileNotFoundError(f"Input FASTA not found: {input_path}")

    records = list(SeqIO.parse(str(input_path), "fasta"))
    if not records:
        raise ValueError("Input FASTA is empty.")

    reference_record = select_reference_record(records)
    reference_protein = cds_to_clean_protein(reference_record)

    mutation_df = generate_single_point_mutations(reference_protein)
    mutation_df.to_csv(args.output_csv, index=False)

    with open(args.reference_fasta, "w") as out_fasta:
        out_fasta.write(f">{reference_record.id} | translated_reference\n")
        out_fasta.write(reference_protein + "\n")

    expected = len(reference_protein) * 19
    print("Selected reference:", reference_record.description)
    print("Reference protein length:", len(reference_protein))
    print("Generated mutations:", len(mutation_df))
    print("Expected mutations (length x 19):", expected)
    print("Saved mutation dataset to:", args.output_csv)
    print("Saved reference protein to:", args.reference_fasta)


if __name__ == "__main__":
    main()
