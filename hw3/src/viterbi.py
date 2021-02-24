import itertools as it
import math
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn import metrics
from tabulate import tabulate
import time

import constants as c
import data_structures as ds
import helper as h
import main as m

"""
perform the viterbi algorithm, iterating n times
"""
def viterbi(sequence, output):
  print("performing viterbi algorithm")
  start_time = time.time()
  seq_len = len(sequence)
  seq_list = list(sequence)

  assert seq_len > 0

  emissions = np.log(c.init_emissions)
  transitions = np.log(c.init_transitions)

  # TODO: consider refactoring this loop outwards, since it's pretty clunky
  for i in range(c.n):
    output.write(
        "\n---Viterbi Parameter Estimation: Iteration {}---\n".format(i + 1))
    # a. the HMM emission/transition parameters used for this pass
    # (e.g., from the tables above for pass 1)
    output.write("emissions\n")
    output.write(tabulate(np.exp(emissions)))
    output.write("\n")
    output.write("transitions\n")
    output.write(tabulate(np.exp(transitions)))
    output.write("\n")

    prob_list = np.zeros((2, seq_len))
    prev_state_list = np.zeros((2, seq_len))

    # initialize using emissions + transitions
    # prob_list = viterbi trellis idea
    # prev_state_list = previous state

    prob_list[0][0] = c.begin_transitions[0] + \
        emissions[0][c.nucleotides.index(sequence[0])]
    prob_list[1][0] = c.begin_transitions[1] + \
        emissions[1][c.nucleotides.index(sequence[0])]

    for j in range(1, seq_len):
      for k in range(len(transitions)):
        prev_scores = prob_list[:, j - 1] + transitions[k]

        prev_state_list[k][j] = np.argmax(prev_scores)
        prob_list[k][j] = emissions[k][c.nucleotides.index(
            sequence[j])] + np.amax(prev_scores)

    final_prob, path = traceback(prob_list, prev_state_list)
    # b. the log probability(natural log, base-e) of the overall Viterbi path
    output.write("final log-prob: " + str(final_prob) + "\n")

    hit_list = m.get_hits(path, i + 1)
    # c. the total number of "hits" found, where a hit is (contiguous)
    # subsequence assigned to state 2 in the Viterbi path
    output.write("hits: " + str(len(hit_list)) + "\n")

    k = c.k
    if i == c.n - 1:
      k = len(hit_list)
    for interval_counter in range(min(k, len(hit_list))):
      hit = list(np.array(hit_list[interval_counter]) + 1)
      # d. the lengths and locations(starting and ending positions)
      # of the first k(defined below) "hits." Print all hits, if there are fewer than k of them.
      output.write(str(hit) + " of length " + str(hit[1] - hit[0]) + "\n")

    emissions = m.update_emissions(path, seq_list)
    transitions = m.update_transitions(path, hit_list)

  end_time = time.time()
  output.write("\nviterbi completed in " +
               str(end_time - start_time) + " seconds\n")
  return hit_list


"""
trace the most likely probability path
"""
def traceback(prob_list, prev_state_list):
  print("doing traceback")
  result = []

  last_index = int(np.argmax(prob_list[:, -1]))
  result.append(last_index)

  # decrement to traceback
  for i in range(len(prev_state_list[0]) - 1, 0, -1):
    last_index = int(prev_state_list[last_index][i])
    result.append(last_index)

  return prob_list[last_index, -1], np.flip(result)
