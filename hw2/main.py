# Michael Huang (mhuang19)
# 1862567

import itertools
import math
import numpy as np
import pandas as pd
import urllib as ul
from tabulate import tabulate

import constants as c
import data_structures as ds
import helper as h

"""
build a count matrix corresponding to a given list of length k sequences.
output: count matrix
"""
def makeCountMatrix(seq_list): 
  num_seqs = len(seq_list)

  assert (num_seqs > 0)

  seq_length = len(seq_list[0])

  result = np.zeros((len(c.nucleotides), seq_length), dtype=float)

  for i in range(num_seqs):
    for j in range(seq_length):
      nuc = seq_list[i][j]
      nuc_index = c.nucleotides.index(nuc)

      result[nuc_index][j] = result[nuc_index][j] + 1

  return result

"""
given a count matrix, and a pseudocount vector, build a new count matrix by adding them.
output: modified count matrix
"""
def addPseudo(count_matrix, pseudocount_vector): 
  assert (len(count_matrix) == len(pseudocount_vector))

  result = np.copy(count_matrix)

  for i in range(len(pseudocount_vector)):
    result[i] = count_matrix[i] + pseudocount_vector[i]

  return result

"""
make a frequency matrix from a count matrix.
output: frequency matrix
"""
def makeFrequencyMatrix(count_matrix): 
  assert (len(count_matrix) == len(c.nucleotides))
  
  result = np.copy(count_matrix)
  total = sum(count_matrix)
  
  for i, row in enumerate(count_matrix):
    for j, count in enumerate(row):
      result[i][j] = count / total[j]

  return result

"""
calculate the entropy of a frequency matrix relative to background.
output: entropy value
"""
def entropy(frequency_matrix, bg_vector): 

  assert (len(frequency_matrix) > 0)
  
  total = 0

  for i, row in enumerate(frequency_matrix):
    for val in row:
      
      if val == 0:
        continue

      entropy = val / bg_vector[i]
      entropy = math.log2(entropy)
      entropy = entropy * val

      total += entropy
  
  return total

"""
make a weight matrix from a frequency matrix and background vector.
output: tuple; tuple[0] = wmm, tuple[1] = entropy
"""
def makeWMM(frequency_matrix, bg_vector): 
  assert (len(frequency_matrix) > 0)

  result = np.copy(frequency_matrix)

  for i, row in enumerate(frequency_matrix):
    for j, frequency in enumerate(row):

      if frequency == 0:
        result[i][j] = -math.inf
        continue;

      val = frequency / bg_vector[i]
      result[i][j] = math.log2(val)
  
  return result, entropy(frequency_matrix, bg_vector)

"""
given a WMM of width k and one or more sequences of varying lengths ≥ k,
scan/score each position of each sequence (excluding those < k from the rightmost end).
output: scores of each valid position 
"""
def scanWMM(wmm, seq_list):
  motif_length = len(wmm[0])
  num_seqs = len(seq_list)

  assert (motif_length > 0 and num_seqs > 0)

  result = []

  for seq in seq_list:
    seq_length = len(seq)

    assert (seq_length - motif_length >= 0)

    num_scores = seq_length - motif_length + 1
    scores = [0] * num_scores

    # iterate through start places
    for i in range(num_scores):
      total = 0
      window = seq[i : i + motif_length]

      # iterate through each sequence from start
      for j in range(len(window)):
        total += wmm[c.nucleotides.index(seq[i + j])][j]
      scores[i] = total

    result.append(scores)
  
  return result

"""
given an(estimated) WMM and a set of sequences, run the E-step of MEME's EM algorithm.
output: matrix of probabilities according to the wmm
"""
def Estep(wmm, seq_set):
  assert (len(wmm) > 0 and len(seq_set) > 0)

  result = []

  scores = scanWMM(wmm, seq_set)
  
  for i, row in enumerate(scores):
    prob_list = []
    total_prob = 0

    # get total prob
    for score in row:
      total_prob = total_prob + 2**score

    # normalize and place probs accordingly
    for j, score in enumerate(row):
      prob_list.append(2**score / total_prob)

    result.append(prob_list)
  
  return result

"""
given the Estep result, re-estimate the WMM. 
output: tuple; tuple[0] = frequency matrix, tuple[1] = wmm, tuple[2] = entropy
"""
def Mstep(seq_list, wmm, estep_result, pseudocount_vector, bg_vector):
  num_seqs = len(seq_list)
  wmm_length = len(wmm[0])

  result = np.zeros((len(c.nucleotides), wmm_length), dtype=float)
  # print(estep_result)

  for i, seq in enumerate(seq_list):
    for j in range(len(seq) - wmm_length + 1):
      window = [seq[j : j + wmm_length]]
      prob = estep_result[i][j]

      count_matrix = makeCountMatrix(window)
      window_count = prob * count_matrix
      result = result + window_count
  
  result = addPseudo(result, pseudocount_vector)
  result = makeFrequencyMatrix(result)
  # print(result)
  wmm, entropy = makeWMM(result, bg_vector)

  return result, wmm, entropy

# perform training step
def train_step(wmm_list, seq_list):
  print("training step")
  freq_result = []
  wmm_result = []
  entropy_result = []

  for wmm in wmm_list:
    freq_trials = []
    wmm_trials = []
    entropy_trials = []

    for i in range(c.trials):
      estep_result = Estep(np.copy(wmm), np.copy(seq_list))
      freq, wmm, entropy = Mstep(np.copy(seq_list), np.copy(wmm), np.copy(estep_result), c.pseudocount_vector, c.bg_vector)

      print(entropy)

      freq_trials.append(freq)
      wmm_trials.append(wmm)
      entropy_trials.append(entropy)

    freq_result.append(freq_trials)
    wmm_result.append(wmm_trials)
    entropy_result.append(entropy)

    # median: use int

# perform evaluation step
def eval_step():
  print("evaluation step")

def main():
  print("hello world")

  # 1. Process fasta into list of fasta data structures
  # This results in a list of fasta data structure lists.
  train_data = h.get_seq_list(h.process_fasta(c.train_fasta))
  eval_data = h.get_seq_list(h.process_fasta(c.eval_fasta))
  
  # 1a. Test as necessary
  # TODO: Note that the training files' sequences aren't the same length
  # test(train_fasta, eval_fasta)

  # 2. Run on training data
  wmm_info_train = []
  freq_matrix_list = []

  init_wmm, init_entropy = h.initialize(train_data[0], c.k)
  train_step(init_wmm, train_data)

  # 3. Run on evaluation data
  

  print("done")

if __name__ == "__main__":
  main()
