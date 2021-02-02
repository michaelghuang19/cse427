# Michael Huang (mhuang19)
# 1862567

import itertools
import math as m
import numpy as np
import pandas as pd
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

  result = np.copy(count_matrix)

  for i in range(len(pseudocount_vector)):
    result[i] = count_matrix[i] + pseudocount_vector[i]

  return result


# make a frequency matrix from a count matrix.
def makeFrequencyMatrix(count_matrix): 
  assert (len(count_matrix) == len(c.nucleotides))
  
  result = np.copy(count_matrix)
  total = sum(count_matrix)
  
  for i in range(len(result[0])):
    for j in range(len(c.nucleotides)):
      result[j][i] = count_matrix[j][i] / total[i]

  return result

# calculate the entropy of a frequency matrix relative to background.
def entropy(frequency_matrix, bg_vector): 
  assert (len(frequency_matrix) > 0)
  
  total = 0

  for i in range(len(frequency_matrix[0])):
    for j in range(len(bg_vector)):
      val = frequency_matrix[j][i]
      
      if val == 0:
        continue

      entropy = val / bg_vector[j]
      entropy = m.log(entropy, 2)
      entropy = entropy * val

      total = total + entropy
  
  return total

# make a weight matrix from a frequency matrix and background vector.
# (It may be convenient to save the entropy with the WMM, too.)
def makeWMM(frequency_matrix, bg_vector): 
  assert (len(frequency_matrix) > 0)

  result = np.copy(frequency_matrix)

  for i in range(len(result[0])):
    for j in range(len(bg_vector)):
      if frequency_matrix[j][i] == 0:
        result[j][i] = -m.inf
        continue;

      val = frequency_matrix[j][i] / bg_vector[j]
      result[j][i] = m.log(val, 2)

  # result = ds.wmm_info(result, entropy(frequency_matrix, bg_vector))

  return result

# given a WMM of width k and one or more sequences of varying lengths ≥ k,
# scan/score each position of each sequence(excluding those < k from the rightmost end).
def scanWMM(wmm, seq_list): 
  motif_length = len(wmm)
  num_seqs = len(seq_list)

  assert (motif_length > 0 and num_seqs > 0)

  result = []

  for seq in seq_list:
    seq_length = len(seq)

    assert (seq_length - motif_length >= 0)

    num_scores = seq_length - motif_length + 1
    scores = [0] * num_scores

    for i in range(num_scores):
      total = 0
      for j in range(i, i + len(wmm)):
        total = total + wmm[c.nucleotides.index(seq[j])][j]
      scores[i] = total
    result.append(scores)
  
  return result

# given an(estimated) WMM and a set of sequences, run the E-step of MEME's EM algorithm;
# i.e., what is E[zij], where zij is the zero-one variable indicating whether the
# motif instance in sequence i begins in position j. (scanWMM does most of the required work.)
def Estep(wmm, seq_set): 
  motif_length = len(wmm)
  num_seqs = len(seq_set)
  
  assert (motif_length > 0 and num_seqs > 0)

  result = []

  for seq in seq_set:

    expectation = [0] * len(seq)

    scores = scanWMM(wmm, list(seq))

    result.append(expectation)
  
  return result



# given the Estep result, re-estimate the WMM. 
# (This is similar to makeCountMatrix / addPseudo / makeFrequencyMatrix / makeWMM,
# with the additional wrinkle that the k-mer inputs to makeCountMatrix each have weights,
# and it may be convenient to generalize makeCountMatrix so that its k-mer inputs are successive, 
# overlapping subsequences of some longer strings, 
# rather than having them explicitly presented as distinct, non-overlapping inputs.)
def Mstep(): 
  print("Mstep")

def test(train_fasta, eval_fasta):
  
  # makeCountMatrix test
  # Note that we use eval_fasta rather than train_fasta since
  # train_files sequences actually don't have the same length
  count_matrix_list = []
  for fasta_list in eval_fasta:
    result = makeCountMatrix(fasta_list)
    # print(result)
    count_matrix_list.append(result)

  # addPseudo test
  pseudo_matrix_list = []
  for count_matrix in count_matrix_list:
    result = addPseudo(count_matrix, c.pseudocount_vector)
    # print(result)
    pseudo_matrix_list.append(result)

  # makeFrequencyMatrix test
  freq_matrix_list = []
  for pseudo_matrix in count_matrix_list:
    result = makeFrequencyMatrix(pseudo_matrix)
    # print(result)
    freq_matrix_list.append(result)

  # entropy test
  entropy_list = []
  for freq_matrix in freq_matrix_list:
    result = entropy(freq_matrix, c.bg_vector)
    # print(result)
    entropy_list.append(result)
  
  # makeWMM test
  wmm_list = []
  for freq_matrix in freq_matrix_list:
    result = makeWMM(freq_matrix, c.bg_vector)
    # print(result)
    wmm_list.append(result)

  # scanWMM test
  for i in range(len(wmm_list)):
    result = scanWMM(wmm_list[i], h.get_sequence_list(eval_fasta[i]))
    print(result)

  # Estep test
  # for i in range(len(wmm_list)):
    # result = Estep(wmm_list[i], h.get_sequence_list(eval_fasta[i]))
    # print(result)

  # Mstep test


def main():
  print("hello world")

  # 1. Process fasta into individual lists.
  # This results in a list of lists, where each list contains those
  # fasta data structures specified within each fasta file.
  train_fasta = []
  eval_fasta = []
  for accession in c.file_dict:
    train_fasta.append(h.process_fasta(accession))
    eval_fasta.append(h.process_fasta(c.file_dict[accession]))
  
  # TODO: Note that the training files' sequences aren't the same length
  test(train_fasta, eval_fasta)

  print("done")

if __name__ == "__main__":
  main()
