# Michael Huang (mhuang19)
# 1862567

import itertools as it
import math
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn import metrics
from tabulate import tabulate

import baum_welch as bw
import constants as c
import data_structures as ds
import helper as h
import viterbi as v

"""
get hits from the path (where we have state 2)
"""
def get_hits(path, iteration):
  print("finding hits for iteration {}".format(iteration))
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
def update_emissions(path, seq_list):
  print("updating emissions")

  state0_indices = np.where(path == 0)[0]
  state1_indices = np.where(path == 1)[0]

  seq_list = pd.Series(seq_list)

  state0_seq = seq_list[state0_indices]
  state1_seq = seq_list[state1_indices]

  state0_probs = h.make_freq_matrix(state0_seq)
  state1_probs = h.make_freq_matrix(state1_seq)

  result = np.array([state0_probs, state1_probs])

  return np.log(result)

"""
given our path, update our transitions probabilities
"""
def update_transitions(path, hit_list):
  print("updating transitions")

  result = np.zeros((2, 2))
  num_intervals = len(hit_list)

  # set state 1 and state 2 transition probability
  # Example: 0011001100, consider boundaries

  state1_total = 0
  state2_total = 0
  
  prev = path[0]
  for i in range(len(path)):
    cur = path[i]

    if prev == 0 and cur == 0:
      result[0][0] += 1
      state1_total += 1
    elif prev == 0 and cur == 1:
      result[0][1] += 1
      state1_total += 1
    elif prev == 1 and cur == 0:
      result[1][0] += 1
      state2_total += 1
    elif prev == 1 and cur == 1:
      result[1][1] += 1
      state2_total += 1

    prev = cur

  result[0] = result[0] / state1_total
  result[1] = result[1] / state2_total

  # trans_11_prob = (state1_total - (2 * num_intervals)) / state1_total
  # trans_12_prob = (2 * num_intervals) / state1_total
  # trans_21_prob = (2 * num_intervals) / state2_total
  # trans_22_prob = (state2_total - (2 * num_intervals)) / state2_total
  
  return np.log(result)

"""
evaluate our list of intervals against the list of 'golden standard' rnas
"""
def evaluate(intervals, ginfo_list, output):
  print("evaluating against golden standard")
  result = []

  output.write(c.disclaimer)
  output.write("\n");
  for interval in intervals:
    result_list = []
    for ginfo in ginfo_list:
      if int(ginfo.start) >= int(interval[0]) and ginfo.end <= int(interval[1]):
        result_list.append(ginfo)
    result.append(result_list)

  for i, item_list in enumerate(result):
    output.write("Matches for " + str(list(np.array(intervals[i]) + 1)))
    output.write("\n")
    for item in item_list:
      output.write(str(item))

  return result

def main():
  print("hello world")

  fasta_list = h.process_fasta(c.genome_file, c.fna_exten)
  seq = h.get_seq_list(fasta_list)[0]

  ginfo_list = h.process_gff(c.genome_file, c.gff_exten)
  ginfo_list = sorted(ginfo_list, key=lambda ginfo: ginfo.start)

  viterbi_output = open(c.results_folder + "viterbi" + c.text_exten, "wt")
  viterbi_intervals = v.viterbi(seq, viterbi_output)
  viterbi_output.close()

  viterbi_eval_output = open(c.results_folder + "viterbi_eval" + c.text_exten, "wt")
  evaluate(viterbi_intervals, ginfo_list, viterbi_eval_output)
  viterbi_eval_output.close()
  
  baum_welch_output = open(c.results_folder + "baum_welch" + c.text_exten, "wt")
  baum_welch_intervals = bw.baum_welch(seq, baum_welch_output)
  baum_welch_output.close()

  baum_welch_eval_output = open(c.results_folder + "baum_welch_eval" + c.text_exten, "wt")
  evaluate(baum_welch_intervals, ginfo_list, baum_welch_eval_output)
  baum_welch_eval_output.close()

  print("done")

if __name__ == "__main__":
  main()
