from Bio import SeqIO

input_file = "prnp_nr_99.fasta"
output_file = "prnp_mammals.fasta"

count = 0

with open(output_file, "w") as out:
    for record in SeqIO.parse(input_file, "fasta"):
        header = record.description
        # crude but effective check
        if any(word in header for word in [
            "Homo", "Bos", "Mus", "Rattus", "Sus", "Pan",
            "Canis", "Felis", "Ovis", "Capra", "Macaca",
            "Callithrix", "Panthera", "Pteropus"
        ]):
            SeqIO.write(record, out, "fasta")
            count += 1

print("Total mammalian sequences:", count)
