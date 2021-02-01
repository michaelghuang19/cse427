# Michael Huang (mhuang19)
# 1862567

import itertools
import pandas as pd
import numpy as np
import urllib as ul
from tabulate import tabulate

import constants as c
import data_structures as ds
import helper as h

# build a count matrix corresponding to a given list of length k sequences.
def makeCountMatrix(seq_list): 
  assert (len(seq_list) > 0)

  seq_length = len(seq_list[0].sequence)

  result = np.zeros((len(c.nucleotides), seq_length), dtype=float)

  for seq in seq_list:
    for i in range(seq_length):
      nuc = seq.sequence[i]
      nuc_index = c.nucleotides.index(nuc)

      result[nuc_index][i] = result[nuc_index][i] + 1

  return result

# given a count matrix, and a pseudocount vector, build a new count matrix by adding them.
def addPseudo(count_matrix, pseudocount_vector): 
  assert (len(count_matrix) == len(pseudocount_vector))

  result = count_matrix

  for i in range(len(pseudocount_vector)):
    count_matrix[i] = count_matrix[i] + pseudocount_vector[i]

  return result


# make a frequency matrix from a count matrix.
def makeFrequencyMatrix(): 
  print("makeFrequencyMatrix")

# calculate the entropy of a frequency matrix relative to background.
def entropy(): 
  print("entropy")

# make a weight matrix from a frequency matrix and background vector.
# (It may be convenient to save the entropy with the WMM, too.)
def makeWMM(): 
  print("makeWMM")

# given a WMM of width k and one or more sequences of varying lengths ≥ k,
# scan/score each position of each sequence(excluding those < k from the rightmost end).
def scanWMM(): 
  print("scanWMM")

# given an(estimated) WMM and a set of sequences, run the E-step of MEME's EM algorithm;
# i.e., what is E[zij], where zij is the zero-one variable indicating whether the
# motif instance in sequence i begins in position j. (scanWMM does most of the required work.)
def Estep(): 
  print("Estep")

# given the Estep result, re-estimate the WMM. 
# (This is similar to makeCountMatrix / addPseudo / makeFrequencyMatrix / makeWMM,
# with the additional wrinkle that the k-mer inputs to makeCountMatrix each have weights,
# and it may be convenient to generalize makeCountMatrix so that its k-mer inputs are successive, 
# overlapping subsequences of some longer strings, 
# rather than having them explicitly presented as distinct, non-overlapping inputs.)
def Mstep(): 
  print("Mstep")

def main():
  print("hello world")

  # 1. Process fasta into individual lists.
  # This results in a list of lists, where each list contains those
  # fasta data structures specified within each fasta file.
  train_files = []
  eval_files = []
  for accession in c.file_dict:
    train_files.append(h.process_fasta(accession))
    eval_files.append(h.process_fasta(c.file_dict[accession]))
  
  # makeCountMatrix test
  count_matrix_list = []
  for fasta_list in eval_files:
    result = makeCountMatrix(fasta_list)
    print(result)
    count_matrix_list.append(result)

  pseudo_matrix_list = []
  for count_matrix in count_matrix_list:
    result = addPseudo(count_matrix, c.pseudocount_vector)
    print(result)
    pseudo_matrix_list.append(result)
  
  print("done")

if __name__ == "__main__":
  main()
