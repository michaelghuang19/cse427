k = 10
seed_proportion = 0.85
trials = 3

bg_vector = [0.25, 0.25, 0.25, 0.25]
pseudocount_vector = [0.25, 0.25, 0.25, 0.25]
nucleotides = ["A", "C", "G", "T"]

# train_fasta = "hw2-debug-train"
# eval_fasta = "hw2-debug-eval"

train_fasta = "hw2-train"
eval_fasta = "hw2-eval"

file_dict = {
  "hw2-debug-train": "hw2-debug-eval",
  "hw2-train": "hw2-eval"
}

fasta_folder = "fasta/"
results_folder = "results/"
fasta_exten = ".fasta"
png_exten = ".png"
text_exten = ".txt"
