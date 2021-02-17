# Michael Huang (mhuang19)
# 1862567

import itertools as it
import math
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn import metrics
from tabulate import tabulate

import constants as c
import data_structures as ds
import helper as h

def viterbi(sequence):
  print("performing viterbi algorithm")
  seq_len = len(sequence)

  assert seq_len > 0

  emissions = c.init_emissions
  transitions = c.init_transitions

  prob_list = np.zeros((2, seq_len))
  state_list = np.zeros((2, seq_len))

  # initialize using given emission/transition
  # prob_list = viterbi trellis idea
  # state_list = actual states

  prob_list[0][0] = transitions[0][0] + \
    emissions[0][c.nucleotides.index(sequence[0])]
  prob_list[1][0] = transitions[0][1] + \
    emissions[1][c.nucleotides.index(sequence[0])]

  for i in range(c.n):
    print("---Viterbi Parameter Estimation: Trial {}---".format(i + 1))

    for j in range(1, seq_len):
      # print(tabulate(emissions))
      # print(tabulate(transitions))

      for k in range(len(state_list)):

        prev_scores = prob_list[:, j - 1] + transitions[k]

        state_list[k][j] = np.argmax(prev_scores)
        prob_list[k][j] = emissions[k][c.nucleotides.index(
            sequence[j])] + np.amax(prev_scores)
    
    # one test run
    break;

  path = traceback(prob_list, state_list)
  print(path)
  hit_list = get_hits(path)

def traceback(prob_list, state_list):
  print("doing traceback")
  result = []
  
  last_index = np.argmax(prob_list[:, -1])
  result.append(last_index)

  # decrement to traceback
  for i in range(len(state_list[0]) - 1, 0, -1):
    last_index = int(state_list[last_index][i])
    result.append(last_index)
  
  return np.flip(result)

def get_hits(path):
  print("finding hits")
  result = []

  # find 1's in path
  # path = np.where(path == 1)
  print(path)

  # group by intervals where we had consecutive 1's

  # for _, group in itertools.groupby(enumerate(hit_locations), key=lambda x: x[1]-x[0]):
  #   group = list(group)
  #   result.append((group[0][1], group[-1][1]))
  #     self.intervals = hit_intervals

  return result

def test():
  print("test")

  ginfo_list = h.process_gff(c.genome_file, c.gff_exten)  # "golden standard"

def main():
  print("hello world")

  fasta_list = h.process_fasta(c.genome_file, c.fna_exten)
  seq = h.get_seq_list(fasta_list)[0]

  # remember to store log probabilities vs probabilities
  viterbi(seq)

  print("done")

if __name__ == "__main__":
  main()
