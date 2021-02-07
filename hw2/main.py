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
  
  for i in range(len(result[0])):
    for j in range(len(c.nucleotides)):
      result[j][i] = count_matrix[j][i] / total[i]

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
output: wmm info struct
"""
def makeWMM(frequency_matrix, bg_vector): 
  assert (len(frequency_matrix) > 0)

  result = np.copy(frequency_matrix)

  for i in range(len(result[0])):
    for j in range(len(bg_vector)):
      if frequency_matrix[j][i] == 0:
        result[j][i] = -math.inf
        continue;

      val = frequency_matrix[j][i] / bg_vector[j]
      result[j][i] = math.log(val, 2)

  result = ds.wmm_info(result, entropy(frequency_matrix, bg_vector))

  return result

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
        total += wmm[c.nucleotides.index(seq[i])][j]
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
  prob_list = []
  for score in scores:
    prob_list.append(2**np.asarray(score))

  for prob in prob_list:
    norm_prob = np.copy(prob)
    
    total_prob = prob.sum()
    norm_prob /= total_prob

    result.append(norm_prob.tolist())

  return result

"""
given the Estep result, re-estimate the WMM. 
output: tuple; tuple[0] = frequency matrix, tuple[1] wmm info struct
"""
def Mstep(seq_list, wmm, estep_result, pseudocount_vector, bg_vector):
  num_seqs = len(seq_list)
  wmm_length = len(wmm[0])

  result = np.zeros((len(c.nucleotides), wmm_length), dtype=float)

  for i, seq in enumerate(seq_list):
    for j in range(len(seq) - wmm_length + 1):
      window = [seq[j : j + wmm_length]]
      prob = estep_result[i][j]

      window_count = prob * makeCountMatrix(window)
      result = result + window_count

  result = addPseudo(result, pseudocount_vector)
  result = makeFrequencyMatrix(result)
  wmm = makeWMM(result, bg_vector)

  return result, wmm

# test suite for working with subroutines
def test(train_fasta, eval_fasta):
  train_list = []
  eval_list = []
  
  for fasta in train_fasta:
    train_list.append(h.get_sequence_list(fasta))
  for fasta in eval_fasta:
    eval_list.append(h.get_sequence_list(fasta))
  
  # makeCountMatrix test
  # Note that we use eval_list rather than train_list since
  # train_files sequences actually don't have the same length
  count_matrix_list = []
  for seq in eval_list:
    result = makeCountMatrix(seq)
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
  wmm_info_list = []
  for freq_matrix in freq_matrix_list:
    result = makeWMM(freq_matrix, c.bg_vector)
    # print(result)
    wmm_info_list.append(result)

  wmm_list = h.get_wmm_list(wmm_info_list)

  # scanWMM test
  score_list = []
  for i in range(len(wmm_list)):
    result = scanWMM(wmm_list[i], eval_list[i])
    # print(result)
    score_list.append(result)

  # Estep test
  estep_list = []
  for i in range(len(wmm_list)):
    result = Estep(wmm_list[i], eval_list[i])
    # print(result)
    estep_list.append(result)

  # Mstep test
  new_wmm_info_list = []
  for i in range(len(eval_fasta)):
    result = Mstep(eval_list[i], wmm_list[i], estep_list[i], c.pseudocount_vector, c.bg_vector)
    # print(result[0])
    # print(result[1].wmm)
    new_wmm_info_list.append(result)

  # initialization test
  result = h.initialize(h.regulate_sequence("ABCDEFGHIJKLMNOPQRSTUVW"), c.k)[3].wmm
  # print(tabulate(result))

# perform training step
def train_step(wmm_list, seq_list):
  print("training step")
  print(len(wmm_list))
  freq_result = []
  wmm_result = []

  # print([wmm.wmm for wmm  in wmm_list])
  print(seq_list[0])
  for struct in wmm_list:
    freq_trials = []
    wmm_trials = []

    wmm_struct = struct

    for i in range(c.trials):
      # print(wmm_struct.wmm)

      estep_result = Estep(wmm_struct.wmm, seq_list)
      mstep_result = Mstep(seq_list, wmm_struct.wmm, estep_result, c.pseudocount_vector, c.bg_vector)
      
      freq = mstep_result[0]
      wmm_struct = mstep_result[1]

      freq_trials.append(freq)
      wmm_trials.append(wmm_struct)

    print([item.entropy for item in wmm_trials])
    wmm_result.append(wmm_trials)

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

  init_list = h.initialize(train_data[0], c.k)
  # print([item.wmm for item in init_list])
  train_step(init_list, train_data)

  # 3. Run on evaluation data
  

  print("done")

if __name__ == "__main__":
  main()
