from Bio import SeqIO
from Bio import pairwise2
from Bio.pairwise2 import format_alignment

records = list(SeqIO.parse("prnp_nr_99.fasta", "fasta"))

reference = records[0]

def percent_identity(aln1, aln2):
    matches = sum(a == b for a, b in zip(aln1, aln2))
    return matches / len(aln1) * 100

with open("global_alignment_results.txt", "w") as out:

    for record in records[1:]:
        alignments = pairwise2.align.globalxx(reference.seq, record.seq)
        best = alignments[0]

        identity = percent_identity(best.seqA, best.seqB)

        out.write(f"\n=== Alignment with {record.id} ===\n")
        out.write(format_alignment(*best))
        out.write(f"Percent Identity: {identity:.2f}%\n")
        out.write(f"Alignment Score: {best.score}\n")
