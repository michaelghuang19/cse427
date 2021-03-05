# various helper functions

import collections
import numpy as np
import pandas as pd

import constants as c
import data_structures as ds

# convert fasta data to a list of formatted data structures
def process_fasta(filename, exten):
  print("processing {} for ".format(exten) + filename)

  fasta_array = []

  file = open(c.data_folder + filename + exten, "r")
  text = file.read()
  seq_list = text.split(">")

  # We need to make a few modifications to properly process.
  # First, we need to start from index 1 if we have multiple sequences in file
  # Second, given a single sequence, the description part is ONLY index 1 if
  # we split by \n, and the rest is going to be part of the sequence
  for i in range(1, len(seq_list)):
    seq = seq_list[i]

    info = seq.split("\n")
    description = info[0]
    sequence = ""
    for j in range(1, len(info)):
      sequence += regulate_sequence(info[j])

    fasta_struct = ds.fasta_info(description, sequence)
    fasta_array.append(fasta_struct)

  return fasta_array

# change fasta by replacing non-whitespace to T's, and making all caps
def regulate_sequence(sequence):
  result = ""
  for char in sequence:
    if (char.upper() in c.nucleotides):
      result += char.upper()
    elif (not char.isspace()):
      result += "T"

  return result

# print list to output as 1-base rather than 0-base
def print_1list(list0, output):
  result = np.array(list0)
  result = list(result + 1)
  output.write(str(result))

# convert gff data to a list of known genes
def process_gff(filename, exten):
  print("processing {} for ".format(exten) + filename)

  result = []

  file = open(c.data_folder + filename + exten, "r")
  text = file.read()
  gene_list = text.split("\n")

  assert len(gene_list) > 0

  gene_list = sorted(gene_list)
  while len(gene_list[0]) == 0 or gene_list[0][0] == "#":
    gene_list.pop(0)

  for gene in gene_list:
    col_vals = gene.split("\t")

    info = ds.gene_info(
        col_vals[0], col_vals[2], col_vals[3], col_vals[4], col_vals[6], col_vals[8])

    result.append(info)

  return result

# count kmers in a sequence
def count_kmers(k, seq_list):
  result = []

  for seq in seq_list:
    result += [seq[i - k: i] for i in range(k, len(seq) + 1)]

  result = collections.Counter(result)

  return result

# get AAGxyT counts, given a kmer + 1 count list
def get_AAGxyT_counts(kmer_plusone_counts):
  result = []
  result = df = pd.DataFrame(columns=c.nucleotides, index=c.nucleotides)

  for i in c.nucleotides:
    for j in c.nucleotides:
      comb = "AAG" + i + j + "T"
      result[j][i] = kmer_plusone_counts[comb]

  return result

# construct an orf summary data structure for printing
def construct_orf_summary(ginfo_list, orf_map, scores):
  print("fuck")


