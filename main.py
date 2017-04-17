# DATABASE INTO STRING

from Bio import SeqIO
record = SeqIO.read("sampleseq.fasta", "fasta")
gene = record.format("fasta")

# ALIGNMENT
# ONLY FINDS THE FIRST MATCH
# STILL TAKES GENE ID (FIRST LINE) INTO ACCOUNT

def findAlignment(query, data):
    if query in data:
        print(data.index(query))
    else:
        print("No match found!")

# SEQUENCE INPUT

query_here = input("Enter your query sequence: ")

findAlignment(query_here, gene)