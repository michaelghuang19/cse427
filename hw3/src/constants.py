# describes any global constants or magic numbers we might need

import numpy as np

# n = 10
n = 2
k = 5
genome_file = "GCF_000091665.1_ASM9166v1_genomic"
# regex for filtering what we want to look for 
# when we compare against the database "golden standard"
# r"gene" should consider everything
eval_filter = r"gene_biotype=.*RNA"

init_emissions = np.array([
    # [A, C, T, G]
    [0.25, 0.25, 0.25, 0.25], # state 1
    [0.20, 0.30, 0.30, 0.20]  # state 2
])
init_transitions = np.array([
    # [state 1, state 2]
    [0.9999, 0.0001],  # begin
    [0.9999, 0.0001],  # state 1
    [0.01, 0.99]       # state 2
])
nucleotides = ["A", "C", "G", "T"]

data_folder = "data/"
results_folder = "results/"
fasta_exten = ".fasta"
fna_exten = ".fna"
gff_exten = ".gff"
png_exten = ".png"
text_exten = ".txt"
