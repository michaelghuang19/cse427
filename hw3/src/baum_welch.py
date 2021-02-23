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
import main as m

"""
perform the baum-welch algorithm + posterior decoding, iterating n times
"""
def baum_welch(sequence, output):
  print("performing baum_welch algorithm")
  seq_len = len(sequence)
  seq_list = list(sequence)

  assert seq_len > 0

  emissions = np.log(c.init_emissions)
  transitions = np.log(c.init_transitions)

  for i in range(c.n):
    output.write(
        "---Baum-Welch Parameter Estimation: Trial {}---\n".format(i + 1))
    # a. the HMM emission/transition parameters used for this pass
    # (e.g., from the tables above for pass 1)
    output.write("emissions\n")
    output.write(tabulate(np.exp(emissions)))
    output.write("\n")
    output.write("transitions\n")
    output.write(tabulate(np.exp(transitions)))
    output.write("\n")

    # initialize using emissions + transitions
    # forward_list = forward probabilities
    # backward_list = backward probabilities
    # prob_list = baum-welch trellis idea

    forward_list, forward_val = get_forward_list(
        emissions, transitions, sequence)
    backward_list = get_backward_list(emissions, transitions, sequence)

    print(forward_list)
    print(forward_val)
    print(backward_list)

    # final_prob_list= forward_list[:, -1]
    # final_prob = 0
    # for prob in final_prob_list:
    #   final_prob += h.log_of_sum_of_logs(final_prob, prob)
    prob_list = (forward_list + backward_list) - forward_val
    # print(prob_list)
    # print(np.exp(prob_list))

    path = np.zeros(seq_len)

    for j in range(0, seq_len):
      if prob_list[1][j] > prob_list[0][j]:
        path[j] = 1

    # print(path)

    # b. the log probability of the genomic input given current params
    output.write("final log-prob: " + str(np.amax(prob_list[:, -1])) + "\n")

    hit_list = m.get_hits(path)
    # c. the total number of "hits" found, where a hit is (contiguous)
    # subsequence assigned to state 2 in the Viterbi path
    output.write("hits: " + str(len(hit_list)) + "\n")

    k = c.k
    if i == c.n - 1:
      k = len(hit_list)
    for interval_counter in range(min(k, len(hit_list))):
      hit = hit_list[interval_counter]
      # d. the lengths and locations(starting and ending positions)
      # of the first k(defined below) "hits." Print all hits, if there are fewer than k of them.
      # (By convention, genomic positions are 1-based, not 0-based, indices.)
      output.write(str(hit) + " of length " + str(hit[1] - hit[0]) + "\n")

    emissions = m.update_emissions(path, seq_list)
    transitions = m.update_transitions(path, hit_list)

  return hit_list

def get_forward_list(emissions, transitions, sequence):
  seq_len = len(sequence)

  assert seq_len > 0

  result = np.zeros((2, seq_len))

  result[0][0] = c.begin_transitions[0] + \
      emissions[0][c.nucleotides.index(sequence[0])]
  result[1][0] = c.begin_transitions[1] + \
      emissions[1][c.nucleotides.index(sequence[0])]

  for i in range(1, len(sequence)):
    for j in range(len(transitions)):
      temp_term0 = result[0][i - 1] + transitions[0][j] + \
          emissions[j][c.nucleotides.index(sequence[i])]
      temp_term1 = result[1][i - 1] + transitions[1][j] + \
          emissions[j][c.nucleotides.index(sequence[i])]

      result[j][i] = h.log_of_sum_of_logs(temp_term0, temp_term1)

  final_vals = result[:, -1]
  forward_val = h.log_of_sum_of_logs(final_vals[0], final_vals[1])

  return result, forward_val


def get_backward_list(emissions, transitions, sequence):
  seq_len = len(sequence)

  assert seq_len > 0

  result = np.zeros((2, seq_len))

  result[0][-1] = np.log(1)
  result[1][-1] = np.log(1)

  for i in range(1, seq_len + 1):
    for j in range(len(transitions)):
      temp_term0 = result[0][-i + 1] + transitions[j][0] + \
          emissions[0][c.nucleotides.index(sequence[-i + 1])]
      temp_term1 = result[1][-i + 1] + transitions[j][1] + \
          emissions[1][c.nucleotides.index(sequence[-i + 1])]

      result[j][-i] = h.log_of_sum_of_logs(temp_term0, temp_term1)

  return result
