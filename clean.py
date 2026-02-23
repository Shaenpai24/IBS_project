from Bio import SeqIO

input_file = "sequence.txt"
output_file = "prnp_unique_cds.fasta"

seen = set()
unique_records = []

for record in SeqIO.parse(input_file, "fasta"):
    seq = str(record.seq)
    if seq not in seen:
        seen.add(seq)
        unique_records.append(record)

SeqIO.write(unique_records, output_file, "fasta")
print("Unique sequences:", len(unique_records))
