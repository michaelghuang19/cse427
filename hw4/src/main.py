# Michael Huang (mhuang19)
# 1862567

import collections
import itertools as it
import math
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statistics
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

def calculate_score(seq, trusted_list, bg_list):

  p_score = calculate_prob(seq, trusted_list[0], trusted_list[1], trusted_list[2])
  q_score = calculate_prob(seq, bg_list[0], bg_list[1], bg_list[2])

  return p_score - q_score

def calculate_prob(seq, kmer_counts, kmer_plusone_counts, kmer_start_counts):

  kmer_count_total = sum(kmer_counts.values())
  start = seq[0:c.k]

  # initialize using start probabilities
  if start in kmer_counts.keys():
    result = np.log(kmer_counts[start] + c.pseudo_count) / \
        (kmer_count_total + (c.pseudo_count * len(c.nucleotides)))
  else:
    result = np.log(c.pseudo_count / 
                    (kmer_count_total + (c.pseudo_count * len(c.nucleotides))))

  for i in range(c.k + 2, len(seq) + 1):
    term = seq[i - c.k - 1: i]
    pre_term = term[0:c.k]

    if pre_term not in kmer_counts.keys() and term not in kmer_plusone_counts.keys():
      result += np.log(c.pseudo_count / (c.pseudo_count * len(c.nucleotides)))
    elif pre_term in kmer_counts.keys() and term not in kmer_plusone_counts.keys():
      result += np.log(c.pseudo_count /
                       (kmer_counts[pre_term] + (c.pseudo_count * len(c.nucleotides))))
    else:
      result += np.log((kmer_plusone_counts[term] + c.pseudo_count) / 
          (kmer_counts[pre_term] + (c.pseudo_count * len(c.nucleotides))))

  return result

def get_AAGxyT_counts(count_list):
  kmer_counts, kmer_plusone_counts, kmer_start_counts = count_list

  result = pd.DataFrame(columns=c.nucleotides, index=c.nucleotides)

  for i in c.nucleotides:
    for j in c.nucleotides:
      pre = "AAG" + i + j
      comb = pre + "T"
      result[j][i] = kmer_plusone_counts[comb] / kmer_counts[pre]

  return result

def plot_roc(key_map, color):
  print("plotting " + color)

  df = pd.DataFrame(key_map, columns=["key", "match"])

  fpr, tpr, thresholds = metrics.roc_curve(df["match"], df["key"])

  max_threshold = [0, 0]
  for i in range(len(thresholds)):
    if (tpr[i] >= 0.8 and thresholds[i] > max_threshold[1]):
      max_threshold[0] = i
      max_threshold[1] = thresholds[i]

  threshold_index = max_threshold[0]

  auc = metrics.roc_auc_score(df["match"], df["key"])

  plt.plot(fpr, tpr, color=color, label="auc = {}".format(auc))
  plt.plot(fpr[threshold_index], tpr[threshold_index], "x", color=color,
           label="threshold at {}".format(thresholds[threshold_index]))
  
  tt_pos, tf_pos = h.get_positives(key_map, thresholds[threshold_index])

  return auc, tt_pos, tf_pos, thresholds[threshold_index]

def plot_flashbulb(master_flashbulb_list, long_flashbulb_list, short_flashbulb_list, threshold, pos_color, neg_color):
  long_df = pd.DataFrame(long_flashbulb_list, columns=[
                         "length", "score", "match"])
  short_df = pd.DataFrame(short_flashbulb_list, columns=[
                          "length", "score", "match"])

  long_colors = [
      pos_color if match else neg_color for match in long_df["match"]]
  short_colors = [
      pos_color if match else neg_color for match in short_df["match"]]

  plt.scatter(long_df["length"], long_df["score"], c=long_colors, s=15)
  plt.scatter(short_df["length"], short_df["score"], c=short_colors, s=15)

  m_long_len, m_long_score = statistics.median(
      long_df["length"]), statistics.median(long_df["score"])
  m_short_len, m_short_score = statistics.median(
      short_df["length"]), statistics.median(short_df["score"])

  plt.scatter(m_short_len, m_short_score,
              marker="x", color="red", s=50, label="A: ({},{})".format(m_short_len, m_short_score))
  plt.scatter(m_long_len, m_long_score,
              marker="x", color="red", s=50, label="B: ({},{})".format(m_long_len, m_long_score))
  
  m_slope = (m_long_score - m_short_score) / (m_long_len - m_short_len)
  m_intercept = m_long_score - (m_slope * m_long_len)

  p_x = m_short_len + (threshold * (m_long_len - m_short_len))
  p_y = (m_slope * p_x) + m_intercept
  p_intercept = -(p_y - (1 / m_slope * p_x))

  plt.scatter(p_x, p_y,
              marker="x", color="red", s=50, label="C: ({},{})".format(p_x, p_y))

  x = np.linspace(-500, 9000, 900)
  plt.plot(x, (m_slope * x) + m_intercept)
  plt.plot(x, (-1 / m_slope * x) + p_intercept)

  plt.ylabel("Markov Score")
  plt.xlabel("ORF Length")
  plt.xlim(-500, 9000)
  plt.ylim(-500, 1000)
  plt.legend()
  plt.tight_layout()
  plt.savefig(c.results_folder + "flashbulb" + c.png_exten)
  plt.close()

  # plot the full dataset now
  master_df = pd.DataFrame(master_flashbulb_list, columns=[
                           "length", "score", "match"])
  
  master_colors = [
      pos_color if match else neg_color for match in master_df["match"]]
  
  plt.scatter(master_df["length"], master_df["score"], c=master_colors, s=15)

  plt.axvline(x=c.short_threshold, color='b', linestyle=":")
  plt.axvline(x=c.long_threshold, color='b', linestyle=":")
  plt.plot(x, (m_slope * x) + m_intercept)
  plt.plot(x, (-1 / m_slope * x) + p_intercept)

  plt.ylabel("Markov Score")
  plt.xlabel("ORF Length")
  plt.xlim(-500, 9000)
  plt.ylim(-500, 1000)
  plt.tight_layout()
  plt.savefig(c.results_folder + "flashbulb_master" + c.png_exten)
  plt.close()

  return m_slope, p_intercept

def get_combined_score(master_flashbulb_list, trend_slope, perp_intercept):
  result = []
  df = pd.DataFrame(master_flashbulb_list, columns=[
      "length", "score", "match"])
  
  for row in df.iterrows():
    x = row[1]["length"]
    y = row[1]["score"]
    trend_intercept = y - (trend_slope * x)

    # quick maths to get the closest point on perpendicular line
    perp_slope = (-1 / trend_slope)
    perp_x = (trend_intercept - perp_intercept / (perp_slope - trend_slope))
    perp_y = (perp_slope * perp_x) + perp_intercept
    
    above_line = (perp_slope * x + perp_intercept) <= y

    # use distance here
    if above_line:
      result.append([math.hypot(x - perp_x, y - perp_y), row[1]["match"]])
    else:
      result.append([-math.hypot(x - perp_x, y - perp_y), row[1]["match"]])

  return result

def main():
  fasta_list = h.process_fasta(c.genome_file, c.fna_exten)
  seq = fasta_list[0].sequence

  ginfo_list, ginfo_ends = h.process_gff(c.genbank_file, c.gff_exten)
  ginfo_list = sorted(ginfo_list, key=lambda ginfo: ginfo.start)

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
  master_orf_seq_map = {key[0]: value for key, value in zip(
      master_orf_locs, master_orf_seqs)}

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
    
    long_seq = trusted_orf_seq_map[long_start]
    first_long_score_map[long_start] = calculate_score(
        long_seq, trusted_list, bg_list)
    
    short_seq = short_orf_seq_map[short_start]
    first_short_score_map[short_start] = calculate_score(
        short_seq, trusted_list, bg_list)

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

  p_probs = get_AAGxyT_counts(trusted_list)
  q_probs = get_AAGxyT_counts(bg_list)

  overall_output.write("\np_probs: \n" + str(p_probs) + "\n")
  overall_output.write("\nq_probs: \n" + str(q_probs) + "\n")

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
  overall_output.write("\n")

  # generating report info step
  print("generating rocs")

  roc_length_map = []
  roc_score_map = []

  long_flashbulb_list = []
  short_flashbulb_list = []

  for start in trusted_orf_seq_map.keys():
    loc = trusted_orf_loc_map[start]
    length = loc[1] - loc[0] + 1
    seq = trusted_orf_seq_map[start]
    score = calculate_score(seq, trusted_list, bg_list)
    end = loc[1] + 4
    match = end in ginfo_ends

    roc_length_map.append([length, match])
    roc_score_map.append([score, match])
    long_flashbulb_list.append([length, score, match])

  for start in short_orf_seq_map.keys():
    loc = short_orf_loc_map[start]
    length = loc[1] - loc[0] + 1
    if length < 4:
      continue
    seq = short_orf_seq_map[start]
    score = calculate_score(seq, trusted_list, bg_list)
    end = loc[1] + 4
    match = end in ginfo_ends

    roc_length_map.append([length, match])
    roc_score_map.append([score, match])
    short_flashbulb_list.append([length, score, match])

  print("plotting flashbulb")

  master_length_map = []
  master_score_map = []
  master_flashbulb_list = []

  for loc in master_orf_locs:
    length = loc[1] - loc[0] + 1
    if length < 4:
      continue
    seq = master_orf_seq_map[loc[0]]
    score = calculate_score(seq, trusted_list, bg_list)
    end = loc[1] + 4
    match = end in ginfo_ends

    master_length_map.append([length, match])
    master_score_map.append([score, match])
    master_flashbulb_list.append([length, score, match])

  trend_slope, perp_intercept = plot_flashbulb(
      master_flashbulb_list, long_flashbulb_list, short_flashbulb_list, 0.20, "orange", "blue")

  master_combined_map = get_combined_score(
      master_flashbulb_list, trend_slope, perp_intercept)

  # do rocs for training (short + long orfs) data

  length_auc, length_t, length_f, length_thresh = plot_roc(roc_length_map, "red")
  score_auc, score_t, score_f, score_thresh = plot_roc(roc_score_map, "green")

  overall_output.write("\ntraining length auc: {}".format(length_auc))
  overall_output.write(
      "\nnumber of length true, false positives: {}, {}\nthreshold:{}\n".format(length_t, length_f, length_thresh))
  overall_output.write("\ntraining score auc: {}".format(score_auc))
  overall_output.write(
      "\nnumber of score true, false positives: {}, {}\nthreshold:{}\n".format(score_t, score_f, score_thresh))
  
  plt.title("ROCs for Training (Short + Long) ORFs")
  plt.ylabel("True Positive Rate")
  plt.xlabel("False Positive Rate")
  plt.legend()
  plt.tight_layout()
  plt.savefig(c.results_folder + "training_rocs" + c.png_exten)

  plt.xlim(-0.05, 0.10)
  plt.ylim(0.95, 1.05)
  plt.tight_layout()
  plt.savefig(c.results_folder + "training_rocs_zoom" + c.png_exten)
  plt.close()

  # now do it for full orf data

  master_length_auc, master_length_t, master_length_f, master_length_thresh = plot_roc(
      master_length_map, "red")
  master_score_auc, master_score_t, master_score_f, master_score_thresh = plot_roc(
      master_score_map, "green")
  master_combined_auc, master_combined_t, master_combined_f, master_combined_thresh = plot_roc(
      master_combined_map, "blue")

  overall_output.write("\nfull length auc: {}".format(master_length_auc))
  overall_output.write("\nnumber of score true, false positives: {}, {}\nthreshold:{}\n".format(
      master_length_t, master_length_f, master_length_thresh))
  overall_output.write("\nfull score auc: {}".format(master_score_auc))
  overall_output.write("\nnumber of length true, false positives: {}, {}\nthreshold:{}\n".format(
      master_score_t, master_score_f, master_score_thresh))
  overall_output.write("\nfull combined auc: {}".format(master_combined_auc))
  overall_output.write("\nnumber of combined true, false positives: {}, {}\nthreshold:{}\n".format(
      master_combined_t, master_combined_f, master_combined_thresh))

  plt.title("ROCs for Full ORFs")
  plt.ylabel("True Positive Rate")
  plt.xlabel("False Positive Rate")
  plt.legend()
  plt.tight_layout()
  plt.savefig(c.results_folder + "full_roc" + c.png_exten)

  plt.xlim(-0.05, 0.10)
  plt.ylim(0.80, 1.05)
  plt.tight_layout()
  plt.savefig(c.results_folder + "full_rocs_zoom" + c.png_exten)
  plt.close()
  
  overall_output.close()

  print("done")

if __name__ == "__main__":
  main()
