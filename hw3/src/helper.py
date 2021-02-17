# various helper functions

import numpy as np
import re

import constants as c
import data_structures as ds
import main as m

# convert fasta data to a list of formatted data structures
def process_fasta(filename, exten):
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

  print("Successfully processed fasta for " + filename)

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

# extracts sequences from fasta info structs
def get_seq_list(fasta_list):
  result = []

  for fasta in fasta_list:
    result.append(fasta.sequence)

  return result

# convert gff data to a list of known genes
def process_gff(filename, exten):
  result = []

  file = open(c.data_folder + filename + exten, "r")
  text = file.read()
  gene_list = text.split("\n")

  assert len(gene_list) > 0

  gene_list = sorted(gene_list)
  while len(gene_list[0]) == 0 or gene_list[0][0] == "#":
    gene_list.pop(0)

  for gene in gene_list:
    is_nc_rna = re.search(r"gene_biotype=.*RNA", gene)

    if is_nc_rna:
      col_vals = gene.split("\t")

      info = ds.gene_info(
          col_vals[0], col_vals[1], col_vals[2], col_vals[3], col_vals[4], col_vals[6], col_vals[8])

      result.append(info)
  
  return result

# build a count matrix corresponding to a given sequence
def makeCountMatrix(seq_list):
  seq_length = len(seq_list)

  assert (seq_length > 0)

  result = np.zeros((len(c.nucleotides), seq_length), dtype=float)

  for i in range(seq_length):
    nuc = seq_list[i]
    nuc_index = c.nucleotides.index(nuc)

    result[nuc_index][j] = result[nuc_index] + 1

  return result
