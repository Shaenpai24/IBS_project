from Bio import Entrez

Entrez.email = "your_email@example.com"

search_term = "PRNP[Gene] AND Mammalia[Organism] AND protein"

# Search once
handle = Entrez.esearch(db="protein", term=search_term, retmax=500)
record = Entrez.read(handle)
ids = record["IdList"]

print(f"Found {len(ids)} sequences")

# Fetch ALL at once (comma separated)
handle = Entrez.efetch(
    db="protein",
    id=",".join(ids),
    rettype="fasta",
    retmode="text"
)

data = handle.read()

with open("prnp_mammalia_raw.fasta", "w") as f:
    f.write(data)

print("Download complete.")
