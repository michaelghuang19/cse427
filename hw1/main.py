import pandas as pd
import numpy as np
import urllib as ul

import fasta
import constants
import tests

# 0. remember to set these constants
element_map = constants.blosum_map
score_matrix = constants.blosum_matrix
# num_epochs = 1000
num_epochs = 50
gap_score = -4

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
      sequence += info[j].upper()

    fasta_struct = fasta.fasta_info(species, accession, description, sequence)
    fasta_array.append(fasta_struct)

  print("Successfully processed fasta for " + filename)

  return fasta_array

# Perform the comparison between lists of fasta structures
def compare_seqs(flist1, flist2, num_permutations):
  for f1 in flist1:
    for f2 in flist2:
      print("Comparing " + f1.accession + " and " + f2.accession)

      # Open the output file and print the respect identifiers
      output = open(constants.results_folder
                    + f1.accession
                    + "_"
                    + f2.accession 
                    + constants.text_exten, "wt")
      output.write(str(f1.accession) + " to " + str(f2.accession) + "\n")
      
      # Create matrix of 0's
      len1 = len(f1.sequence)
      len2 = len(f2.sequence)
      matrix = np.zeros((len1 + 1, len2 + 1), dtype = int)

      # Fill out the matrix
      best_align = local_align(list(f1.sequence), list(f2.sequence), matrix)

      # print optimal score
      output.write("\nscore: " + str(best_align.best_score) + "\n\n")
      
      # print optimal alignment information
      alignment = backtrack(matrix.tolist(), best_align, list(f1.sequence), list(f2.sequence))

      # Line up our final strings for printing
      print_alignment(str(f1.accession), str(f2.accession), alignment, output)

      # Print the matrix if it's short enough
      if (len(alignment[0]) < 15 and len(alignment[1]) < 15):
        print_matrix(matrix, output)
        output.write("\n")

      # Generate random permutations for sequences, and then compare them to the
      # original sequence to determine the p-value.
      p_count = 0
      for i in range(0, num_permutations):
        random_matrix = np.zeros((len1 + 1, len2 + 1), dtype=int)
        random_seq = generate_permutations(f2.sequence)

        random_score = local_align(f1.sequence, random_seq, random_matrix).best_score
        
        if (random_score > best_align.best_score):
          p_count = p_count + 1
      
      empirical_p = (p_count + 1) / (num_permutations + 1)
      output.write("p-value: " + "{:e}".format(empirical_p))

      output.close()

# Helper function for performing the local alignment between 2 sequences
def local_align(seq1, seq2, matrix):
  best = 0;
  bestx = 0;
  besty = 0;
  len1 = len(seq1) + 1
  len2 = len(seq2) + 1

  for i in range(0, len1):
    for j in range(0, len2):
      if (i == 0 or j == 0):
        continue

      cur_score = score_matrix[element_map.index(seq1[i-1])][element_map.index(seq2[j-1])]

      matrix[i][j] = max(matrix[i-1][j-1] + cur_score,
                        matrix[i-1][j] + gap_score,
                        matrix[i][j-1] + gap_score,
                        0)

      if (matrix[i][j] > best):
        best = matrix[i][j]
        bestx = i
        besty = j
  
  return fasta.align_info(best, bestx, besty)

# Helper function for backtracking through a matrix from best location
# Extra credit: Processes and later prints the similarities in between strings
def backtrack(matrix, best_align, seq1, seq2):
  alignment1 = []
  alignment2 = []
  similarities = []
  curx = best_align.best_coord1
  cury = best_align.best_coord2
  xchange = True
  ychange = True

  while curx > 0 and cury > 0 and matrix[curx][cury] > 0:
    print(str(curx) + ", " + str(cury) + ": " + str(matrix[curx][cury]))

    if xchange:
      alignment1.append(seq1[curx-1])
    else:
      alignment1.append("-")

    if ychange:
      alignment2.append(seq2[cury-1])
    else:
      alignment2.append("-")

    last1 = alignment1[len(alignment1) - 1]
    last2 = alignment2[len(alignment2) - 1]
    if (last1 == "-" or last2 == "-"):
      similarities.append(" ")
    elif (last1 == last2):
      similarities.append(last1)
    elif (score_matrix[element_map.index(last1)][element_map.index(last2)] > 0):
      similarities.append("+")
    else:
      similarities.append(" ")
    
    oldScore = score_matrix[element_map.index(seq1[curx-1])][element_map.index(seq2[cury-1])]

    maxVal = max(matrix[curx-1][cury-1] + oldScore,
            matrix[curx-1][cury] + gap_score,
            matrix[curx][cury-1] + gap_score)
    print(maxVal)

    if (matrix[curx-1][cury-1] + oldScore == maxVal):
      curx = curx - 1
      cury = cury - 1
      xchange = True
      ychange = True
      continue
    elif (matrix[curx][cury - 1] + gap_score == maxVal):
      cury = cury - 1
      xchange = False
      ychange = True
      continue
    elif (matrix[curx-1][cury] + gap_score == maxVal):
      curx = curx - 1
      xchange = True
      ychange = False
      continue

  return (''.join(alignment1)[::-1], ''.join(alignment2)[::-1], curx, cury, ''.join(similarities)[::-1])

# Helper function for printing out alignment
def print_alignment(f1_id, f2_id, alignment, output):
  cur = 0
  end = 0
  length = len(alignment[0])

  seq1 = alignment[0]
  seq2 = alignment[1]

  index1 = alignment[2]
  index2 = alignment[3]

  similarities = alignment[4]

  seq1_section = ""
  seq2_section = ""

  while cur < length:
    if (length - cur < 60):
      end = length
    else:
      end = cur + 60

    seq1_section = seq1[cur:end]
    seq2_section = seq2[cur:end]

    seq1_info = f1_id + ":  " + str(index1 + 1) + "  "
    seq2_info = f2_id + ":  " + str(index2 + 1) + "  "
    while (len(seq1_info) < len(seq2_info)):
      seq1_info = seq1_info + " "
    while (len(seq2_info) < len(seq1_info)):
      seq2_info = seq2_info + " "

    output.write(seq1_info + seq1_section + "\n")
    output.write((" " * len(seq1_info)) + similarities[cur:end] + "\n")
    output.write(seq2_info + seq2_section + "\n")
    output.write("\n")

    index1 = index1 + 60 - seq1_section.count("-")
    index2 = index2 + 60 - seq2_section.count("-")
    cur = length + 60

# Helper function for printing out a matrix, since numpy is being annoying
def print_matrix(matrix, output):
  if (matrix.shape[0] > 0 and matrix.shape[1] > 0):
    for line in matrix:
      output.write("[")
      output.write(str(line[0]))
      for i in range(1, len(line)):
        output.write(" " + str(line[i]))
      output.write("]")
      output.write("\n")

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
  # 1. Get our desired fasta files. Comment if you don't need to.
  # get_fasta()

  # 2. Process fasta into individual lists. This results in a list of lists,
  # where each list contains corresponding fasta data structures for each
  # individual accession file, should a file contain multiple sequences.
  # You shouldn't really be commenting this in normal cases.
  # fasta_list = []
  # for accession in constants.accession_set:
  #   fasta_list.append(process_fasta(accession))
  
  # 2a. Uncomment to see what sequences we have processed in console,
  # grouped by file in fasta_list. This won't work unless 2. is uncommented.
  #   print("\n")
  #   for seq in file:
  #     print(seq.species)
  #     print(seq.accession)
  #     print(seq.description)
  #     print(seq.sequence)
  
  # 3. Perform analysis using the Smith-Waterman sequence alignment algorithm
  # This is the actual analysis work, so it shouldn't be commented.
  # k = len(fasta_list)
  # for i in range(0, k):
  #   for j in range(i + 1, k):
  #     compare_seqs(fasta_list[i], fasta_list[j], num_epochs)

  # 3a. Uncomment for a example test seqs, that will be output as separate files.
  tests.site_tests()

  # 4. Perform test suite as defined in the assignment spec
  tests.test1()

if __name__ == "__main__":
  main()
