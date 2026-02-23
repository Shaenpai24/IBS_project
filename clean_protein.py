from Bio import SeqIO

input_file = "sequence.fasta"
filtered = []

for record in SeqIO.parse(input_file, "fasta"):
    if 200 <= len(record.seq) <= 300:
        filtered.append(record)

SeqIO.write(filtered, "prnp_len_filtered.fasta", "fasta")
print("After length filtering:", len(filtered))
