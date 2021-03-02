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

def perform_counts(orf_struct_list):
  print("performing counts")

  for orf in orf_struct_list:
    orf.find_bg_seqs()

    kmer_counts = h.count_kmers(2, orf.orf_list)
    kmer_total = sum(kmer_counts.values())

    print(kmer_counts)
    print(kmer_total)

    # TODO: count normal and long orfs, and compare respectively

def main():
  print("hello world")

  fasta_list = h.process_fasta(c.genome_file, c.fna_exten)
  seq = fasta_list[0].sequence

  # With the test, we expect the following output:
  # ATA-ACG, CCC
  # CGT-GAC
  # AAC-GTG-ACC

  orf_struct_list, trusted_orf_list, hyp_orf_list = scan_long_seqs(seq)

  perform_counts(trusted_orf_list)

  markov_orf_output = open(c.results_folder + "markov_orf" + c.text_exten, "wt")
    # remember pseudocounts of 1
  markov_orf_output.close()

  # ginfo_list = h.process_gff(c.genome_file, c.gff_exten)
  # ginfo_list = sorted(ginfo_list, key=lambda ginfo: ginfo.start)

  print("done")

if __name__ == "__main__":
  main()
