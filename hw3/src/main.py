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

  for i in range(c.n):
    print("---Viterbi Parameter Estimation: Trial {}---".format(i + 1))
    print(tabulate(emissions))
    print(tabulate(transitions))

    prob_list = np.zeros((2, seq_len))
    prev_state_list = np.zeros((2, seq_len))

    # initialize using emissions + transitions
    # prob_list = viterbi trellis idea
    # prev_state_list = previous state

    prob_list[0][0] = transitions[0][0] + \
      emissions[0][c.nucleotides.index(sequence[0])]
    prob_list[1][0] = transitions[0][1] + \
      emissions[1][c.nucleotides.index(sequence[0])]

    for j in range(1, seq_len):
      for k in range(1, len(transitions)):
        prev_scores = prob_list[:, j - 1] + transitions[k]

        prev_state_list[k - 1][j] = np.argmax(prev_scores)
        prob_list[k - 1][j] = emissions[k - 1][c.nucleotides.index(
            sequence[j])] + np.amax(prev_scores)
    
    path = traceback(prob_list, prev_state_list)
    hit_list = get_hits(path)

    emissions = update_emissions(path, sequence)
    transitions = update_transitions(path, hit_list)
    
    # one test run
    break;

def traceback(prob_list, prev_state_list):
  print("doing traceback")
  result = []
  
  last_index = np.argmax(prob_list[:, -1])
  result.append(last_index)

  # decrement to traceback
  for i in range(len(prev_state_list[0]) - 1, 0, -1):
    last_index = int(prev_state_list[last_index][i])
    result.append(last_index)
  
  return np.flip(result)

def get_hits(path):
  print("finding hits")
  result = []

  # find 1's in path
  hits = np.where(path == 1)[0]

  # group by intervals where we had consecutive 1's
  # thanks geeksforgeeks
  for key, group in it.groupby(enumerate(hits), key=lambda t: t[1] - t[0]):
    group = list(group)
    result.append([group[0][1], group[-1][1]])

  return result

def update_emissions(path, sequence):
  result = np.zeros((2, 4))

  # set state 1 and state 2 emission probability

  return result

def update_transitions(path, hit_list):
  result = np.zeros((3, 2))
  
  # how do we set the beginning transition probability?
  result[0] = c.transitions[0]

  # set state 1 and state 2 transition probability

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
