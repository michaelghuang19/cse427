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

def viterbi_training(sequence):
  print("viterbi training")
  seq_len = len(sequence)

  assert seq_len > 0

  emissions = c.init_emissions
  transitions = c.init_transitions

  prob_list = np.zeros((2, seq_len))
  state_list = np.zeros((2, seq_len))

  # initialize using given emission/transition
  # prob_list = viterbi trellis idea
  # state_list = actual states

  prob_list[0][0] = c.transitions[0][0] + \
    c.emissions[0][c.nucleotides.index(sequence[0])]
  prob_list[1][0] = c.transitions[0][1] + \
    c.emissions[1][c.nucleotides.index(sequence[0])]

  for i in range(c.n):
    print("---Viterbi Parameter Estimation: Trial {}---".format(i + 1))

    for j, nuc in enumerate(sequence):
      nuc = 0


    print(c.init_emissions)
    print(c.init_transitions)
    
    # one test run
    break;

def hmm_viterbi():
  print("performing hmm viterbi algorithm")

def traceback():
  print("doing traceback")

def test():
  print("test")

  ginfo_list = h.process_gff(c.genome_file, c.gff_exten)  # "golden standard"

def main():
  print("hello world")

  fasta_list = h.process_fasta(c.genome_file, c.fna_exten)
  seq = h.get_seq_list(fasta_list)[0]

  # remember to store log probabilities vs probabilities
  viterbi_training(seq)

  # test()

  print("done")

if __name__ == "__main__":
  main()
