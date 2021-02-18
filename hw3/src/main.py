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

"""
perform the viterbi algorithm, iterating n times
"""
def viterbi(sequence):
  print("performing viterbi algorithm")
  seq_len = len(sequence)

  assert seq_len > 0

  emissions = c.init_emissions
  transitions = c.init_transitions

  # TODO: consider refactoring this loop outwards, since it's pretty clunky
  for i in range(c.n):
    print("---Viterbi Parameter Estimation: Trial {}---".format(i + 1))
    # a. the HMM emission/transition parameters used for this pass (e.g., from the tables above for pass 1),
    print("emissions")
    print(tabulate(np.exp(emissions)))
    print("transitions")
    print(tabulate(np.exp(transitions)))

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
    
    final_prob, path = traceback(prob_list, prev_state_list)
    # b. the log probability(natural log, base-e) of the overall Viterbi path,
    print(final_prob)
  
    hit_list = get_hits(path)
    # c. the total number of "hits" found, where a hit is (contiguous) subsequence assigned to state 2 in the Viterbi path, and
    print(len(hit_list))
    
    k = c.k
    if i == c.n - 1:
      k = len(hit_list)
    for hit in hit_list:
      # d. the lengths and locations(starting and ending positions) of the first k(defined below) "hits." Print all hits, if there are fewer than k of them. (By convention, genomic positions are 1-based, not 0-based, indices.)
      print(str(hit) + " of length " + str(hit[1] - hit[0]))

    emissions = update_emissions(path, sequence)
    transitions = update_transitions(path, hit_list)

    # run only once for testing
    break;

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

"""
get hits from the path (where we have state 2)
"""
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

"""
given our path, update our emissions probabilities
"""
def update_emissions(path, sequence):
  state0_seq = ""
  state1_seq = ""

  for i, element in enumerate(path):
    if element == 0:
      state0_seq += sequence[i]
    if element == 1:
      state1_seq += sequence[i]
  
  state0_count_matrix = h.make_count_matrix(state0_seq)
  state1_count_matrix = h.make_count_matrix(state1_seq)

  state0_probs = h.make_freq_matrix(state0_count_matrix)
  state1_probs = h.make_freq_matrix(state1_count_matrix)

  result =  np.array([state0_probs, state1_probs])

  return np.log(result)

"""
given our path, update our transitions probabilities
"""
def update_transitions(path, hit_list):
  result = np.zeros((3, 2))
  num_intervals = len(hit_list)
  
  result[0] = c.init_transitions[0]

  # set state 1 and state 2 transition probability
  # Example: 0011001100, consider boundaries

  state1_total = 0
  state2_total = 0
  
  prev = path[0]
  for i in range(1, len(path)):
    cur = path[i]

    if prev == 0 and cur == 0:
      result[1][0] += 1
      state1_total += 1
    elif prev == 0 and cur == 1:
      result[1][1] += 1
      state1_total += 1
    elif prev == 1 and cur == 0:
      result[2][0] += 1
      state2_total += 1
    elif prev == 1 and cur == 1:
      result[2][1] += 1
      state2_total += 1

    prev = cur

  result[1] = result[1] / state1_total
  result[2] = result[2] / state2_total

  # trans_11_prob = (state1_total - (2 * num_intervals)) / state1_total
  # trans_12_prob = (2 * num_intervals) / state1_total
  # trans_21_prob = (2 * num_intervals) / state2_total
  # trans_22_prob = (state2_total - (2 * num_intervals)) / state2_total

  return np.log(result)

"""
evaluate our list of intervals against the list of 'golden standard' rnas
"""
def evaluate(intervals, ginfo_list):
  print("evaluating against golden standard")
  result = [[] * len(intervals)]

  for ginfo in ginfo_list:
    for i, interval in enumerate(intervals):
      if ginfo.start >= interval[0] and ginfo.end <= interval[1]:
        result[i].append(ginfo)

  for i, item_list in enumerate(result):
    print(intervals[i])
    for item in item_list:
      print(item)

  return result

def main():
  print("hello world")

  fasta_list = h.process_fasta(c.genome_file, c.fna_exten)
  seq = h.get_seq_list(fasta_list)[0]

  intervals = viterbi(seq)

  ginfo_list = h.process_gff(c.genome_file, c.gff_exten)
  ginfo_list = sorted(ginfo_list, key=lambda ginfo: ginfo.start)
  evaluate(intervals, ginfo_list)

  print("done")

if __name__ == "__main__":
  main()
