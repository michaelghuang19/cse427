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

def scan_long_seqs(seq):
  orf_struct_list = []
  trusted_orf_list = []
  hyp_orf_list = []

  for i in range(3):
    long_orf_output = open(
        c.results_folder + "long_orf{}".format(i) + c.text_exten, "wt")

    long_orf = ds.orf_struct(seq)
    long_orf.find_orf_locs(i)

    long_orf_locs = []
    short_orf_locs = []

    long_orf_list = []
    short_orf_list = []

    for j in range(len(long_orf.orf_locs)):
      locj = long_orf.orf_locs[j]
      orfj = long_orf.orf_list[j]

      if (len(orfj) > c.long_threshold):
        h.print_1list(locj, long_orf_output)
        long_orf_output.write("\t" + str(orfj) + "\n")

        long_orf_locs.append(locj)
        long_orf_list.append(orfj)
      else:
        short_orf_locs.append(locj)
        short_orf_list.append(orfj)

    orf_struct_list.append(long_orf)
    trusted_orf_list.append(ds.orf_struct(seq, long_orf_locs, long_orf_list))
    hyp_orf_list.append(ds.orf_struct(seq, short_orf_locs, short_orf_list))

    long_orf_output.close()
  
  return orf_struct_list, trusted_orf_list, hyp_orf_list

def perform_counts(orf_list):
  print("performing counts")

  kmer_counts = h.count_kmers(c.k, orf_list)
  kmer_plusone_counts = h.count_kmers(c.k + 1, orf_list)
  
  kmer_start_counts = h.count_starts(c.k, orf_list)

  return [kmer_counts, kmer_plusone_counts, kmer_start_counts]

def calculate_probs(trusted_list, bg_list):
  print("calculating probabilities")

  # kmer_total = sum(kmer_counts.values())
  # kmer_plusone_total = sum(kmer_plusone_counts.values())
  # kmer_start_total = sum(kmer_start_counts.values())



def main():
  print("hello world")

  fasta_list = h.process_fasta(c.genome_file, c.fna_exten)
  seq = fasta_list[0].sequence

  # With the test file, we expect the following output:
  # ATA-ACG, CCC
  # CGT-GAC
  # AAC-GTG-ACC

  orf_struct_list, trusted_orf_list, hyp_orf_list = scan_long_seqs(seq)

  trusted_seq_list = []
  bg_seq_list = []

  for orf in trusted_orf_list:
    orf.find_bg_seqs()

    trusted_seq_list += orf.orf_list
    bg_seq_list += orf.bg_seqs

  trusted_list = perform_counts(trusted_seq_list)
  bg_list = perform_counts(bg_seq_list)

  calculate_probs(trusted_list, bg_list)

  markov_orf_output = open(c.results_folder + "markov_orf" + c.text_exten, "wt")
    # remember pseudocounts of 1
  markov_orf_output.close()

  # ginfo_list = h.process_gff(c.genome_file, c.gff_exten)
  # ginfo_list = sorted(ginfo_list, key=lambda ginfo: ginfo.start)

  print("done")

if __name__ == "__main__":
  main()
