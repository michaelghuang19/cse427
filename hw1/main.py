import pandas as pd
import numpy as np
import csv
import urllib as ul
# import Bio

import fasta
import constants

# remember to re-set these constants
element_map = constants.test_map
score_matrix = constants.test_matrix
num_epochs = 100
gap_score = -1;

# extra credit j: automated fasta retrieval + write into fasta folder
def get_fasta():
  for item in constants.accession_set:
    url = constants.uniprot_url + item + constants.fasta_exten
    f = ul.request.urlopen(url)
    fasta = f.read()

    output = open(constants.fasta_folder + item + constants.fasta_exten, "wb")
    output.write(fasta)
    output.close()

    print("Successfully retrieved the fasta file for " + item)

# convert fasta data to a list of formatted data structures
def process_fasta(filename):
  fasta_array = []

  file = open(constants.fasta_folder + filename + constants.fasta_exten, "r")
  text = file.read()
  seq_list = text.split(">")

  # We need to make a few modifications to properly process.
  # First, we need to start from index 1 if we have multiple sequences in file
  # Second, given a single sequence, the description part is ONLY index 1 if
  # we split by \n, and the rest is going to be part of the sequence
  for i in range(1, len(seq_list)):
    seq = seq_list[i]

    info = seq.split("\n")
    species, accession, description = info[0].split("|")
    sequence = ""
    for j in range(1, len(info)):
      sequence += info[j]

    fasta_struct = fasta.fasta_info(species, accession, description, sequence)
    fasta_array.append(fasta_struct)

  print("Successfully processed fasta for " + filename)

  return fasta_array

# Perform the comparison between lists of fasta structures
def compare_seqs(flist1, flist2, num_permutations):
  for f1 in flist1:
    for f2 in flist2:
      # TODO: open file for printing

      print("Comparing " + f1.accession + " and " + f2.accession)
      
      # Create matrix of 0's
      len1 = len(f1.sequence)
      len2 = len(f2.sequence)
      matrix = np.zeros((len1 + 1, len2 + 1), dtype = int)

      # Fill out the matrix
      align_info = local_align(list(f1.sequence), list(f2.sequence), matrix)
      print(matrix)

      # print optimal score
      print(align_info[0])
      
      # print optimal alignment
      # start from last location, and backtrack
      print(align_info[1], align_info[2], sep=",")

# helper function for performing the local alignment between 2 sequences
def local_align(seq1, seq2, matrix):
  best = 0;
  bestx = 0;
  besty = 0;
  len1 = len(seq1) + 1
  len2 = len(seq2) + 1

  # TODO: dp here
  for i in range(0, len1):
    for j in range(0, len2):
      if (i == 0 or j == 0):
        continue

      cur_score = score_matrix[element_map[seq1[i-1]]][element_map[seq2[j-1]]]

      matrix[i][j] = max(matrix[i-1][j-1] + cur_score,
                        matrix[i-1][j] + gap_score,
                        matrix[i][j-1] + gap_score,
                        0)

      if (matrix[i][j] > best):
        best = matrix[i][j]
        bestx = i
        besty = j
  
  return (best, bestx, besty)

# Helper function for generating permutations, as detailed in lecture
def generate_permutations(sequence):
  seq_letters = list(sequence)

  for i in range(len(seq_letters) - 1, 0, -1):
    j = np.random.randint(i)
    first = seq_letters[i]
    seq_letters[i] = seq_letters[j]
    seq_letters[j] = first

  return ''.join(seq_letters)

def main():
  # 0. Ensure our print is working, and set our constants at top of file
  print("hello world!")

  # 1. Get our desired fasta files. Uncomment as necessary
  # get_fasta()

  # 2. Process fasta into individual lists. This results in a list of lists,
  # where each list contains corresponding fasta data structures for each
  # individual accession file, should a file contain multiple sequences.
  # fasta_list = []
  # for accession in constants.accession_set:
  #   fasta_list.append(process_fasta(accession))
  
  # 2a. Uncomment to see what sequences we have processed, grouped by file
  # for file in fasta_list:
  #   print("\n")
  #   for seq in file:
  #     print(seq.species)
  #     print(seq.accession)
  #     print(seq.description)
  #     print(seq.sequence)
  
  # 3. Perform analysis using the Smith-Waterman sequence alignment algorithm
  # k = len(fasta_list)
  # for i in range(0, k):
  #   for j in range(i + 1, k):
  #     compare_seqs(fasta_list[i], fasta_list[j], num_epochs)

  # example: xxxcde and abcxdex
  compare_seqs([fasta.fasta_info("", "abc", "", "abcxdex")],
  [fasta.fasta_info("", "xxx", "", "xxxcde")], 0)

  # . Generate random permutations for sequences, and then compare them to the
  # original sequence to determine the p-value.
  # for i in range(0, num_epochs):
  #   random_seq = generate_permutations("abcd")
  #   compare_seqs(sequence, random_seq)
  #   if compare_seqs > og_val:
  #     p_count = p_count + 1
  # empirical_p = ((float) p_count) / ((float) num_epochs)

if __name__ == "__main__":
  main()
