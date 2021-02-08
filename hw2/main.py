# Michael Huang (mhuang19)
# 1862567

import math
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn import metrics
from tabulate import tabulate

import constants as c
import data_structures as ds
import helper as h

"""
build a count matrix corresponding to a given list of length k sequences.
returns: count matrix
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
returns: modified count matrix
"""
def addPseudo(count_matrix, pseudocount_vector): 
  assert (len(count_matrix) == len(pseudocount_vector))

  result = np.copy(count_matrix)

  for i in range(len(pseudocount_vector)):
    result[i] = count_matrix[i] + pseudocount_vector[i]

  return result

"""
make a frequency matrix from a count matrix.
returns: frequency matrix
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
returns: entropy value
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
returns: tuple; tuple[0] = wmm, tuple[1] = entropy
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
returns: scores of each valid position 
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
returns: matrix of probabilities according to the wmm
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
returns: tuple; tuple[0] = frequency matrix, tuple[1] = wmm, tuple[2] = entropy
"""
def Mstep(seq_list, wmm, estep_result, pseudocount_vector, bg_vector):
  num_seqs = len(seq_list)
  wmm_length = len(wmm[0])

  result = np.zeros((len(c.nucleotides), wmm_length), dtype=float)

  for i, seq in enumerate(seq_list):
    for j in range(len(seq) - wmm_length + 1):
      window = [seq[j : j + wmm_length]]
      prob = estep_result[i][j]

      count_matrix = makeCountMatrix(window)
      window_count = prob * count_matrix
      result = result + window_count
  
  result = addPseudo(result, pseudocount_vector)
  result = makeFrequencyMatrix(result)
  wmm, entropy = makeWMM(result, bg_vector)

  return result, wmm, entropy

# perform training step
def train_step(wmm_list, seq_list, init_entropy):
  print("training step")
  freq_result = []
  wmm_result = []
  entropy_result = []
  
  print("generating ABC from initial 3 E/M pairs")

  for wmm in wmm_list:
    for i in range(c.trials):
      estep_result = Estep(wmm, seq_list)
      freq, wmm, entropy = Mstep(seq_list, wmm, estep_result, c.pseudocount_vector, c.bg_vector)
    
    freq_result.append(freq)
    wmm_result.append(wmm)
    entropy_result.append(entropy)

  ABCD_wmm = []
  ABCD_freq = []

  A_index = np.argmax(entropy_result)
  B_index = np.argmin(entropy_result)
  C_index = np.argsort(entropy_result)[len(entropy_result)//2]

  ABCD_wmm.append(wmm_result[A_index])
  ABCD_wmm.append(wmm_result[B_index])
  ABCD_wmm.append(wmm_result[C_index])

  ABCD_freq.append(freq_result[A_index])
  ABCD_freq.append(freq_result[B_index])
  ABCD_freq.append(freq_result[C_index])

  print("generating D from running 7 additional E/M pairs on ABC")

  freq_result = []
  wmm_result = []
  entropy_result = []

  last_freq_result = []
  last_wmm_result = []
  last_entropy_result = []

  for wmm in wmm_list:
    freq_trials = []
    wmm_trials = []
    entropy_trials = []

    # run 10 iterations
    for i in range(c.k):
      estep_result = Estep(wmm, seq_list)
      freq, wmm, entropy = Mstep(seq_list, wmm, estep_result, c.pseudocount_vector, c.bg_vector)

      freq_trials.append(freq)
      wmm_trials.append(wmm)
      entropy_trials.append(entropy)

    freq_result.append(freq_trials)
    wmm_result.append(wmm_trials)
    entropy_result.append(entropy_trials)

    last_freq_result.append(freq)
    last_wmm_result.append(wmm)
    last_entropy_result.append(entropy)

  D_index = np.argmax(last_entropy_result)
  ABCD_wmm.append(last_wmm_result[D_index])
  ABCD_freq.append(last_freq_result[D_index])
  
  print("writing training data to output")
  
  # write data to output
  for i in range(len(entropy_result)):
    entropy_result[i].insert(0, init_entropy[0])

  output = open(c.results_folder + "train_step" + c.text_exten, "wt")
  output.write("entropies\n")
  output.write(tabulate(entropy_result))
  output.write("\n")
  for i, freq in enumerate(ABCD_freq):
    output.write(h.get_letter(i) + "frequency \n")
    output.write(tabulate(freq))
  output.close()

  return ABCD_wmm, ABCD_freq

# perform evaluation step
# returns: tuple; tuple[0] = ABCD WMMs, tuple[1] = ABCD frequency matrices 
def eval_step(ABCD_wmm, ABCD_freq, seq_list):
  print("evaluation step")

  score_list = []
  start_list = []

  print("evaluating models and plotting histograms")
  for i in range(len(ABCD_wmm)):
    wmm = ABCD_wmm[i]
    freq = ABCD_freq[i]

    scores = scanWMM(wmm, seq_list)
    values = [0] * len(scores[0])

    for j, row in enumerate(scores):
      values[np.argmax(row)] += 1

    start = np.argmax(values)

    score_list.append(scores)
    start_list.append(start)

    # plotting stuff
    plt.bar(range(len(scores[0])), values)
    plt.title("Most common location of best motif hit: " + str(start))
    plt.savefig(c.results_folder + "wmm_plot_{}.png".format(h.get_letter(i)), dpi = 200)
    plt.close()
  
  print("calculating auc and creating roc")
  # do roc and auc stuff
  label_list = []
  auc_list = []
  break_point = [0, 1, 0]
  for i, score_struct in enumerate(score_list):
    ss_pairs = [] 
    counts = []

    for index in score_struct:
      row_list = []
      for element in index:
        row_list.append(element)
      counts.append(row_list)
    
    # we'd rather be able to directly index by start position, so swap axes
    counts = np.swapaxes(counts, 0, 1)

    flat_list = h.flatten_2d_list(counts)
    flat_list = sorted(flat_list)

    for score in flat_list:
      correct_matrix = counts >= score

      motif_section = correct_matrix[start_list[i] - 1]
      nonmotif_section = np.delete(correct_matrix, start_list[i] - 1, 0)

      motif_section = motif_section.flatten()
      nonmotif_section = nonmotif_section.flatten()

      true_positives = np.sum(motif_section)
      false_positives = np.sum(nonmotif_section)
      false_negatives = np.sum(~motif_section)
      true_negatives = np.sum(~nonmotif_section)
      
      tpr = true_positives / (true_positives + false_negatives)
      fpr = false_positives / (false_positives + true_negatives)

      ss_pairs.append([tpr, fpr])

      # see if we can update our special point
      if tpr == 1 and i == 2 and fpr < break_point[1]:
        break_point = [tpr, fpr, score]

    ss_pairs = sorted(ss_pairs, key = lambda pair: [pair[0], pair[1]])
    x_values = [point[0] for point in ss_pairs]
    y_values = [point[1] for point in ss_pairs]
    auc = metrics.auc(x_values, y_values)

    auc_list.append(auc)
    label_list.append(h.get_letter(i) + ": (" + str(auc) + ")")

    plt.plot(x_values, y_values)

  print("creating roc")
  # plot true positive rate vs false positive rate
  plt.plot([0, 1], [0, 1], linestyle = ":")
  plt.scatter(break_point[0], break_point[1])
  plt.legend(label_list)
  plt.ylabel("True Positive Rate")
  plt.xlabel("False Positive Rate")
  plt.savefig(c.results_folder + "roc.png", dpi = 200)
  plt.close()

  # write test output
  output = open(c.results_folder + "eval_step" + c.text_exten, "wt")
  output.write("auc list for [A, B, C, D]\n")
  output.write(str(auc_list))
  output.write("\nconcrete point: [tpr, fpr, score]\n")
  output.write(str(break_point))
  output.close()

def main():
  # 1. Process fasta into list of fasta data structures
  # This results in a list of fasta data structure lists.
  train_data = h.get_seq_list(h.process_fasta(c.train_fasta))
  eval_data = h.get_seq_list(h.process_fasta(c.eval_fasta))

  # 2. Train using training data
  wmm_info_train = []
  freq_matrix_list = []

  init_wmm, init_entropy = h.initialize(train_data[0], c.k)
  ABCD_wmm, ABCD_freq = train_step(init_wmm, train_data, init_entropy)

  # 3. Run on evaluation data
  eval_step(ABCD_wmm, ABCD_freq, eval_data)

  print("done")

if __name__ == "__main__":
  main()
