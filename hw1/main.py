# Michael Huang (mhuang19)
# 1862567

import pandas as pd
import numpy as np
import urllib as ul
from tabulate import tabulate
import itertools

import fasta
import constants
import tests

# 0. Remember to set these constants as necessary
element_map = constants.blosum_map
score_matrix = constants.blosum_matrix
final_matrix = []
num_epochs = 0
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
def compare_seqs(flist1, flist2, num_permutations, output):
  best_score = 0

  for f1 in flist1:
    for f2 in flist2:
      print("Comparing " + f1.accession + " and " + f2.accession)

      # Open the output file and print the respect identifiers.
      # Uncomment this and the lower .close() out if you want to print to individual files

      # output.close()
      # output = open(constants.results_folder
      #               + f1.accession
      #               + "-"
      #               + f2.accession 
      #               + constants.text_exten, "wt")
      output.write(str(f1.accession) + " to " + str(f2.accession) + "\n")
      
      # Create matrix of 0's
      len1 = len(f1.sequence)
      len2 = len(f2.sequence)
      matrix = np.zeros((len1 + 1, len2 + 1), dtype = int)

      # Fill out the matrix
      best_align = local_align(list(f1.sequence), list(f2.sequence), matrix)

      # print optimal score
      best_score = best_align.best_score
      output.write("\nscore: " + str(best_score) + "\n\n")
      
      # print optimal alignment information
      alignment = backtrack(matrix.tolist(), best_align, list(f1.sequence), list(f2.sequence))

      # Line up our final strings for printing
      print_alignment(str(f1.accession), str(f2.accession), alignment, output)

      # Print the matrix if it's short enough
      if (len(f1.sequence) < 15 and len(f2.sequence) < 15):
        print_matrix(matrix, output, list("  " + f2.sequence), list(" " + f1.sequence))
        output.write("\n")

      # Generate random permutations for sequences, and then compare them to the
      # original sequence to determine the p-value.
      if (num_permutations > 0):
        print_pvalue(f1, f2, best_align, num_permutations, output)

      # Uncomment this if you uncommented the block above to print to individual files.
      # output.close()
  
  return best_score

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
  coords = []
  alignment1 = []
  alignment2 = []
  similarities = []
  curx = best_align.best_coord1
  cury = best_align.best_coord2

  while curx > 0 and cury > 0 and matrix[curx][cury] > 0:
    coords.insert(0, (curx, cury))
    
    oldScore = score_matrix[element_map.index(seq1[curx-1])][element_map.index(seq2[cury-1])]

    maxVal = max(matrix[curx-1][cury-1] + oldScore,
            matrix[curx-1][cury] + gap_score,
            matrix[curx][cury-1] + gap_score)

    if (matrix[curx-1][cury-1] + oldScore == maxVal):
      curx = curx - 1
      cury = cury - 1
    elif (matrix[curx][cury - 1] + gap_score == maxVal):
      cury = cury - 1
    elif (matrix[curx-1][cury] + gap_score == maxVal):
      curx = curx - 1

  # for loop here to complete the sequence
  for i in range(len(coords)):

    if (i != 0 and coords[i-1][0] == coords[i][0]):
      alignment1.append("-")
    else:
      alignment1.append(seq1[coords[i][0]-1])

    if (i != 0 and coords[i-1][1] == coords[i][1]):
      alignment2.append("-")
    else:
      alignment2.append(seq2[coords[i][1]-1])

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

  return (''.join(alignment1), ''.join(alignment2), curx, cury, ''.join(similarities))

# Helper function for printing out alignments (relatively) neatly
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
    cur = cur + 60

# Helper function for printing out a matrix, since numpy is being annoying
def print_matrix(matrix, output, columns, rows):
  if (matrix.shape[0] > 0 and matrix.shape[1] > 0):
    matrix = matrix.tolist()

    for i in range(len(matrix)):
      matrix[i].insert(0, rows[i])

    matrix.insert(0, list("-" * len(columns)))
    matrix.insert(0, columns)
    

    output.write(tabulate(matrix))

    output.write("\n")

# Helper function for printing out the p-value
def print_pvalue(f1, f2, best_align, num_permutations, output):
  p_count = 0
  for i in range(0, num_permutations):
    random_matrix = np.zeros((len(f1.sequence) + 1, len(f2.sequence) + 1), dtype=int)
    random_seq = generate_permutations(f2.sequence)

    random_score = local_align(f1.sequence, random_seq, random_matrix).best_score
    
    if (random_score > best_align.best_score):
      p_count = p_count + 1

  empirical_p = (p_count + 1) / (num_permutations + 1)
  output.write("p-value: " + "{:e}".format(empirical_p) + "\n\n")

# Helper function for generating permutations, using as detailed in lecture
def generate_permutations(sequence):
  seq_letters = list(sequence)

  for i in range(len(seq_letters) - 1, 0, -1):
    j = np.random.randint(i)
    first = seq_letters[i]
    seq_letters[i] = seq_letters[j]
    seq_letters[j] = first

  return ''.join(seq_letters)

# Main function for running the alignment algorithm
def main():
  # 1. Get our desired fasta files. Comment if you don't need to/already did,
  # this can take a bit (check if the fasta folder contains fasta files already).
  get_fasta()

  # 2. Process fasta into individual lists. This results in a list of lists,
  # where each list contains corresponding fasta data structures for each
  # individual accession file, should a file contain multiple sequences.
  # You shouldn't really be commenting this out in normal cases.
  fasta_list = []
  for accession in constants.accession_set:
    fasta_list.append(process_fasta(accession))
  
  final_matrix = np.zeros((len(fasta_list), len(fasta_list)), dtype=int)
  
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
  # For the purposes of the assignment, we open output.txt rather than individual files.
  output = open("output" + constants.text_exten, "wt")

  # Test 1.
  tests.test1(output, 999)

  # Tests 2. and 3.
  k = len(fasta_list)
  for i in range(0, k):
    for j in range(i + 1, k):
      if ((fasta_list[i][0].accession == "P15172" and fasta_list[j][0].accession == "Q10574")
          or (fasta_list[i][0].accession == "P15172" and fasta_list[j][0].accession == "O95363")):
        final_matrix[i][j] = compare_seqs(fasta_list[i], fasta_list[j], 999, output)
      else:
        final_matrix[i][j] = compare_seqs(fasta_list[i], fasta_list[j], num_epochs, output)
      
      output.write("\n")

  # 4. Write the final matrix as well before closing output
  seqnum_list = [*range(1, len(final_matrix) + 1)]
  columns = seqnum_list.insert(0, " ")
  print_matrix(final_matrix, output, columns, seqnum_list)

  output.close()

# Extra credit f: Helper for finding out exactly how many possible permutations
# outscore the given deadly/ddgearlyk pair
def ec_exact():
  
  output = open(constants.results_folder + "ec_exact" + constants.text_exten, "wt")

  threshold = tests.test1(output, 0)

  count = 0
  num_perms = 0

  s1_list = list(itertools.permutations("DEADLY"))
  s2_list = list(itertools.permutations("DDGEARLYK"))

  for s1 in s1_list:
    for s2 in s2_list:
        score = compare_seqs([fasta.fasta_info("", "s1", "", ''.join(s1))],
                          [fasta.fasta_info("", "s2", "", ''.join(s2))],
                          0,
                          output)
        
        num_perms = num_perms + 1
        if (score > threshold):
          count = count + 1

  output.write(str(count))
  output.write(str(count / num_perms))

  output.close()

# Extra credit f: Helper for testing how often random permutations outscore the 
# deadly/ddgearlyk pair
def ec_random():
  output = open(constants.results_folder + "ec_random" + constants.text_exten, "wt")

  threshold = tests.test1(output, 0)

  # this will print p-values to the file, which we can then accumulate in report
  for i in range(20):
    tests.test1(output, 10000)

  output.close()


if __name__ == "__main__":
  # main()
  ec_exact()
  ec_random()
