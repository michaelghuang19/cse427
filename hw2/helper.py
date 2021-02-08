# Various helper functions
import numpy as np

import constants as c
import data_structures as ds
import main as m

# convert fasta data to a list of formatted data structures
def process_fasta(filename):
  fasta_array = []

  file = open(c.fasta_folder + filename + c.fasta_exten, "r")
  text = file.read()
  seq_list = text.split(">")

  # We need to make a few modifications to properly process.
  # First, we need to start from index 1 if we have multiple sequences in file
  # Second, given a single sequence, the description part is ONLY index 1 if
  # we split by \n, and the rest is going to be part of the sequence
  for i in range(1, len(seq_list)):
    seq = seq_list[i]

    info = seq.split("\n")
    description = info[0]
    sequence = ""
    for j in range(1, len(info)):
      sequence += regulate_sequence(info[j])

    fasta_struct = ds.fasta_info(description, sequence)
    fasta_array.append(fasta_struct)

  print("Successfully processed fasta for " + filename)

  return fasta_array

# change fasta by replacing non-whitespace to T's, and making all caps
def regulate_sequence(sequence):
  result = ""
  for char in sequence:
    if (char.upper() in c.nucleotides):
      result += char.upper()
    elif (not char.isspace()):
      result += "T"

  return result

# extracts sequences from fasta info structs
def get_seq_list(fasta_list):
  result = []

  for fasta in fasta_list:
    result.append(fasta.sequence)
  
  return result

# em initialization step
# output: tuple; tuple[0] = list of wmms, tuple[1] = list of entropies
def initialize(sequence, k):
  i = 0
  increment = int(k / 2)
  init_pseudo_count = np.full((4), ((1 / c.seed_proportion) - 1) / 4)

  wmm_result = []
  entropy_result = []

  while i + k <= len(sequence):
    window = [sequence[i : i + k]]

    count_matrix = m.makeCountMatrix(window)
    count_matrix = m.addPseudo(count_matrix, init_pseudo_count)

    freq_matrix = m.makeFrequencyMatrix(count_matrix)
    wmm, entropy = m.makeWMM(freq_matrix, c.bg_vector)
    
    wmm_result.append(wmm)
    entropy_result.append(entropy)

    i += increment

  # check if we missed a possible spot, and slide back accordingly
  diff = (i - len(sequence))
  if diff % increment != 0:
    i = len(sequence) - k
    window = [sequence[i : i + k]]

    count_matrix = m.makeCountMatrix(window)
    count_matrix = m.addPseudo(count_matrix, init_pseudo_count)

    freq_matrix = m.makeFrequencyMatrix(count_matrix)
    wmm, entropy = m.makeWMM(freq_matrix, c.bg_vector)

    wmm_result.append(wmm)
    entropy_result.append(wmm)

  return wmm_result, entropy_result

# flattens a 2d list into a long singular list
def flatten_2d_list(input_list):
  return [item for sublist in input_list for item in sublist]

# get corresponding letter to model number
def get_letter(num):
  return chr(ord('@') + num + 1)
