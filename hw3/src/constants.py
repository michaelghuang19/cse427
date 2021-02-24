# describes any global constants or magic numbers we might need

import numpy as np

n = 10
k = 5
genome_file = "GCF_000091665.1_ASM9166v1_genomic"

# regex for filtering what we want to look for 
# when we compare against the database "golden standard"

# This is fast, and only keeps the relevant RNAs
eval_filter = r"gene_biotype=.*RNA"

# In our standardized test cases, r"NC" should consider everything,
# since everything relevant has "NC" somewhere in it
# eval_filter = r"gene"

init_emissions = np.array([
    # [A, C, G, T]
    [0.25, 0.25, 0.25, 0.25], # state 1
    [0.20, 0.30, 0.30, 0.20]  # state 2
])
init_transitions = np.array([
    # [state 1, state 2]
    [0.9999, 0.0001],  # state 1
    [0.01, 0.99]       # state 2
])
begin_transitions = np.array(
    [0.9999, 0.0001]  # begin
)
nucleotides = ["A", "C", "G", "T"]

data_folder = "data/"
results_folder = "results/"
fasta_exten = ".fasta"
fna_exten = ".fna"
gff_exten = ".gff"
png_exten = ".png"
text_exten = ".txt"

disclaimer = "***Note that these eval results are not fully representative, since they only consider " + \
    "those elements that are fully contained within the interval itself; " + \
    "see report.txt for select detailed analyses***\n"
