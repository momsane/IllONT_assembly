import os
from Bio import SeqIO
import pandas as pd

# Input file paths
gff_file = snakemake.input.gff
faa_file = snakemake.input.faa
output_fld = snakemake.params.fld
output_ordered_faa = snakemake.output[0]

print(f"Output file path: {output_ordered_faa}")

print("Parsing files")

# Step 1: Parse the GFF file and extract gene IDs with their start positions
gene_positions = []

with open(gff_file, "r") as gff:
    for line in gff:
        if line.startswith("#") or not line.strip():
            continue  # Skip comments or empty lines
        parts = line.strip().split("\t")
        if len(parts) < 9 or parts[2] != "CDS":  # Focus on coding sequences
            continue
        
        start = int(parts[3])  # Start position of the gene
        attributes = parts[8]
        
        # Extract gene ID (assuming 'ID=' field is present)
        gene_id = None
        for attr in attributes.split(";"):
            if attr.startswith("ID="):
                gene_id = attr.split("=")[1]
                break
        
        if gene_id:
            gene_positions.append((gene_id, start))

print("Sorting genes")

# Sort genes by start position
gene_positions_sorted = sorted(gene_positions, key=lambda x: x[1])

# Step 2: Create a dictionary of sequences from the genes.faa file
faa_sequences = SeqIO.to_dict(SeqIO.parse(faa_file, "fasta"))

print("Reordering genes")

# Step 3: Reorder the sequences based on the sorted gene IDs
ordered_sequences = []
for gene_id, _ in gene_positions_sorted:
    if gene_id in faa_sequences:
        ordered_sequences.append(faa_sequences[gene_id])
    else:
        print(f"Warning: Gene ID {gene_id} not found in {faa_file}")

# Step 4: Write the ordered sequences to a new FASTA file

print("Writing outputs")

if not os.path.exists(output_fld):
      os.makedirs(output_fld)
      print(f"Created output folder: {output_ordered_faa}")
    
SeqIO.write(ordered_sequences, output_ordered_faa, "fasta")
print(f"Ordered FASTA file saved to: {output_ordered_faa}")
