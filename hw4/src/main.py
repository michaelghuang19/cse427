# Michael Huang (mhuang19)
# 1862567

import collections
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

      if len(orfj) > c.long_threshold:
        h.print_1list(locj, long_orf_output)
        long_orf_output.write("\t" + str(orfj) + "\n")

        long_orf_locs.append(locj)
        long_orf_list.append(orfj)
      elif len(orfj) < c.short_threshold:
        short_orf_locs.append(locj)
        short_orf_list.append(orfj)

    orf_struct_list.append(long_orf)
    trusted_orf_list.append(ds.orf_struct(seq, long_orf_locs, long_orf_list))
    hyp_orf_list.append(ds.orf_struct(seq, short_orf_locs, short_orf_list))

    long_orf_output.close()
  
  return orf_struct_list, trusted_orf_list, hyp_orf_list

def perform_counts(orf_list):

  kmer_counts = h.count_kmers(c.k, orf_list)
  kmer_plusone_counts = h.count_kmers(c.k + 1, orf_list)

  kmer_start_counts = collections.Counter([seq[0:c.k] for seq in orf_list])

  return [kmer_counts, kmer_plusone_counts, kmer_start_counts]

def calculate_score(locs, seq, trusted_list, bg_list):

  p_score = calculate_prob(seq, trusted_list[0], trusted_list[1], trusted_list[2])
  q_score = calculate_prob(seq, trusted_list[0], trusted_list[1], trusted_list[2])
  
  return p_score - q_score

def calculate_prob(seq, kmer_counts, kmer_plusone_counts, kmer_start_counts):
  kmer_key_num = len(kmer_counts.keys())
  kmer_start_total = sum(kmer_start_counts.values())

  # initialize using start probabilities
  result = np.log(c.pseudo_count / (c.pseudo_count * kmer_key_num))
  if seq[0:c.k] in kmer_start_counts.keys():
    result = np.log(kmer_start_counts[seq[0:c.k]] + c.pseudo_count) / \
        (kmer_start_total + c.pseudo_count * kmer_key_num)

  for i in range(c.k + 2, len(seq) + 1):
    term = seq[i - c.k - 1 : i]

    if term[:-1] not in kmer_counts.keys() and term not in kmer_plusone_counts.keys():
      result += np.log(c.pseudo_count / (c.pseudo_count * kmer_key_num))
    elif term[:-1] in kmer_counts.keys() and term not in kmer_plusone_counts.keys():
      result += np.log(c.pseudo_count /
                       (kmer_counts[term[:-1]] + (c.pseudo_count * kmer_key_num)))
    else:
      result += np.log((kmer_plusone_counts[term] + c.pseudo_count) / (
          kmer_counts[term[:-1]] + (c.pseudo_count * kmer_key_num)))
  
  return result

def evaluate(output):
  output.write()

def main():
  print("hello world")

  fasta_list = h.process_fasta(c.genome_file, c.fna_exten)
  seq = fasta_list[0].sequence

  # With the test file, we expect the following output:
  # ATA-ACG, CCC
  # CGT-GAC
  # AAC-GTG-ACC

  # scan in sequences
  print("scanning in sequences")

  orf_struct_list, trusted_orf_list, hyp_orf_list = scan_long_seqs(seq)

  # perform necessary k-mer counts
  print("performing necessary k-mer counts")

  master_orf_locs = []
  master_orf_seq_list = []

  trusted_orf_seq_list = []
  bg_orf_seq_list = []

  short_orf_locs = []
  short_orf_seq_list = []

  for i in range(len(orf_struct_list)):
    master_orf_locs += orf_struct_list[i].orf_locs
    master_orf_seq_list += orf_struct_list[i].orf_list

    trusted_orf_list[i].find_bg_seqs()
    trusted_orf_seq_list += trusted_orf_list[i].orf_list
    bg_orf_seq_list += trusted_orf_list[i].bg_seqs

    short_orf_locs += hyp_orf_list[i].orf_locs
    short_orf_seq_list += hyp_orf_list[i].orf_list

  trusted_list = perform_counts(trusted_orf_seq_list)
  bg_list = perform_counts(bg_orf_seq_list)

  # perform probability calculations
  print("calculating probabilities")
    
  for loc, seq in zip(master_orf_locs, master_orf_seq_list):
    score = calculate_score(loc, seq, trusted_list, bg_list)

  markov_orf_output = open(c.results_folder + "markov_orf" + c.text_exten, "wt")
    # remember pseudocounts of 1
  markov_orf_output.close()

  # ginfo_list = h.process_gff(c.genome_file, c.gff_exten)
  # ginfo_list = sorted(ginfo_list, key=lambda ginfo: ginfo.start)

  # evaluation step
  print("performing evaluation")
  overall_output = open(c.results_folder + "overall_results" + c.text_exten, "wt")
  # evaluate(overall_output)
  overall_output.write("total orfs: " + str(len(master_orf_seq_list)) + "\n")
  overall_output.write("long orfs: " + str(len(trusted_orf_seq_list)) + "\n")
  overall_output.write("short orfs: " + str(len(short_orf_seq_list)) + "\n")
  # for each reading frame: total number, first + length / last + length
  # total number of CDS strands
  # shortest orfs, including start/end, length, score, matches
  # longest orfs, including start/end, length, score, matches
  # p-count
  # q-count
  overall_output.close()

  print("done")

if __name__ == "__main__":
  main()
