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
  q_score = calculate_prob(seq, bg_list[0], bg_list[1], bg_list[2])
  
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

def main():
  fasta_list = h.process_fasta(c.genome_file, c.fna_exten)
  seq = fasta_list[0].sequence

  ginfo_list, ginfo_ends = h.process_gff(c.genbank_file, c.gff_exten)
  ginfo_list = sorted(ginfo_list, key=lambda ginfo: ginfo.start)

  # With the test file, we expect the following output:
  # ATA-ACG, CCC
  # CGT-GAC
  # AAC-GTG-ACC

  # scan in sequences
  print("scanning in sequences")

  orf_struct_list, trusted_orf_list, short_orf_list = scan_long_seqs(seq)

  # perform necessary k-mer counts
  print("performing necessary k-mer counts")  

  master_orf_locs = []
  master_orf_seqs = []

  trusted_orf_locs = []
  trusted_orf_seqs = []
  bg_orf_seqs = []

  short_orf_locs = []
  short_orf_seqs = []

  for i in range(len(orf_struct_list)):

    master_orf_locs += orf_struct_list[i].orf_locs
    master_orf_seqs += orf_struct_list[i].orf_list

    trusted_orf_list[i].find_bg_seqs()
    trusted_orf_locs += trusted_orf_list[i].orf_locs
    trusted_orf_seqs += trusted_orf_list[i].orf_list
    bg_orf_seqs += trusted_orf_list[i].bg_seqs

    short_orf_locs += short_orf_list[i].orf_locs
    short_orf_seqs += short_orf_list[i].orf_list
  
  trusted_orf_seq_map = {key[0] : value for key, value in zip(
      trusted_orf_locs, trusted_orf_seqs)}
  short_orf_seq_map = {key[0] : value for key, value in zip(
      short_orf_locs, short_orf_seqs)}

  trusted_orf_loc_map = {loc[0]: loc for loc in trusted_orf_locs}
  short_orf_loc_map = {loc[0]: loc for loc in short_orf_locs}

  trusted_list = perform_counts(trusted_orf_seqs)
  bg_list = perform_counts(bg_orf_seqs)

  # perform probability calculations
  print("calculating probabilities")
  
  first_long_score_map = {}
  first_short_score_map = {}
  
  for long_start, short_start in zip(sorted(trusted_orf_seq_map.keys())[0:c.k],
      sorted(short_orf_seq_map.keys())[0:c.k]):
    
    long_loc = trusted_orf_loc_map[long_start]
    long_seq = trusted_orf_seq_map[long_start]
    first_long_score_map[long_start] = calculate_score(
        long_loc, long_seq, trusted_list, bg_list)
    
    short_loc = short_orf_loc_map[short_start]
    short_seq = short_orf_seq_map[short_start]
    first_short_score_map[short_start] = calculate_score(
        short_loc, short_seq, trusted_list, bg_list)

  # evaluation step
  print("performing evaluation")
  overall_output = open(c.results_folder + "overall_results" + c.text_exten, "wt")

  # for each reading frame: total number, first + length / last + length
  for i in range(len(orf_struct_list)):
    overall_output.write("\nreading frame at offset {}\n".format(i))
    orf_locs = orf_struct_list[i].orf_locs
    overall_output.write("total orfs: {}\n".format(len(orf_locs)))

    overall_output.write("first orf: {}, length: {}\n".format(
        h.one_format_interval(orf_locs[0]),
        orf_locs[0][1] - orf_locs[0][0] + 1))
    overall_output.write("last orf: {}, length: {}\n".format(
        h.one_format_interval(orf_locs[len(orf_locs) - 1]),
        orf_locs[len(orf_locs) - 1][1] - orf_locs[len(orf_locs) - 1][0] + 1))

  overall_output.write("\ntotal orfs: {}\n".format(len(master_orf_locs)))
  overall_output.write("long orfs: {}\n".format(len(trusted_orf_locs)))
  overall_output.write("short orfs: {}\n".format(len(short_orf_locs)))
  
  overall_output.write("\nsimple plus strand CDSs: {}\n".format(len(ginfo_list)))

  p_counts = h.get_AAGxyT_counts(trusted_list[1])
  q_counts = h.get_AAGxyT_counts(bg_list[1])

  overall_output.write("\np_counts: \n" + str(p_counts) + "\n")
  overall_output.write("\nq_counts: \n" + str(q_counts) + "\n")

  first_long_matches = [loc[1] + 4 in ginfo_ends for loc in
      [trusted_orf_loc_map[key] for key in first_long_score_map.keys()]]
  first_short_matches = [loc[1] + 4 in ginfo_ends for loc in
      [short_orf_loc_map[key] for key in first_short_score_map.keys()]]

  for i, long_start in enumerate(first_long_score_map.keys()):
    overall_output.write("\nstart/end: {}\tlength: {}".format(
        h.one_format_interval(trusted_orf_loc_map[long_start]),
        trusted_orf_loc_map[long_start][1] - trusted_orf_loc_map[long_start][0] + 1))
    overall_output.write("\tscore: {}".format(
        first_long_score_map[long_start]))
    overall_output.write("\tmatch: {}".format(first_long_matches[i]))

  overall_output.write("\n")

  for i, short_start in enumerate(first_short_score_map.keys()):
    overall_output.write("\nstart/end: {}\tlength: {}".format(
        h.one_format_interval(short_orf_loc_map[short_start]),
        short_orf_loc_map[short_start][1] - short_orf_loc_map[short_start][0] + 1))
    overall_output.write("\tscore: {}".format(
        first_short_score_map[short_start]))
    overall_output.write("\tmatch: {}".format(first_short_matches[i]))

  overall_output.close()

  # generating report info step
  print("generating report info")
  report_output = open(c.results_folder + "report" + c.text_exten, "wt")
  report_output.close()

  master_score_map = {}
  master_orf_loc_map = {loc[0]: loc for loc in master_orf_locs}
  master_orf_seq_map = {key[0]: value for key, value in zip(short_orf_locs, short_orf_seqs)}

  for start in sorted(master_orf_seq_map.keys()):
    loc = master_orf_loc_map[start]
    seq = master_orf_seq_map[start]
    master_score_map[start] = calculate_score(loc, seq, trusted_list, bg_list)

  print("done")

if __name__ == "__main__":
  main()
