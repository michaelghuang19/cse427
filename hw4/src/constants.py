# Michael Huang (mhuang19)
# 1862567

# describes any global constants or magic numbers we might need

import numpy as np

# kth-order markov model
k = 5
short_threshold = 50
long_threshold = 1400
# long_threshold = 0
genome_file = "GCF_000091665.1_ASM9166v1_genomic"
# genome_file = "test"
genbank_file = "plusgenes-subset"
pseudo_count = 1

# regex for filtering what we want from the "golden standard" database
# This is fast, and only keeps the relevant RNAs
eval_filter = r"gene_biotype=.*RNA"

nucleotides = ["A", "C", "G", "T"]
nuc_map = {"A":"T", "C":"G", "G":"C", "T":"A"}
stop_codons = ["TAA", "TAG", "TGA"]

data_folder = "data/"
results_folder = "results/"
fna_exten = ".fna"
gff_exten = ".gff"
text_exten = ".txt"
png_exten = ".png"
