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

    # changes start here
    prob_list = np.zeros((2, seq_len))
    path = np.zeros(seq_len)

    # initialize using emissions + transitions
    # prob_list = baum-welch trellis idea
    # prev_state_list = previous state

    # check 0.5
    prob_list[0][0] = transitions[0][0] + \
      emissions[0][c.nucleotides.index(sequence[0])]
    prob_list[1][0] = transitions[0][1] + \
      emissions[1][c.nucleotides.index(sequence[0])]

    for j in range(1, seq_len):
      for k in range(1, len(transitions)):
        total = 0
        for m in range(1, len(transitions)):
      
      # update path
      # path[k - 1][j] = np.argmax(prev_scores)

      # figure out how to some emissions properly
      # prob_list[k - 1][j] = emissions[k - 1][c.nucleotides.index(sequence[j])] + np.amax(prev_scores)

    # run only once for testing
    break

  return None
