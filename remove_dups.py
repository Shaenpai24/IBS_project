from Bio import SeqIO
unique = []
seen = set()

for record in SeqIO.parse("prnp_len_filtered.fasta", "fasta"):
    seq = str(record.seq)
    if seq not in seen:
        seen.add(seq)
        unique.append(record)

SeqIO.write(unique, "prnp_unique.fasta", "fasta")
print("Unique sequences:", len(unique))
