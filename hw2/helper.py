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

# extracts wmms from wmm info structs
def get_wmm_list(wmm_list):
  result = []

  for wmm in wmm_list:
    result.append(wmm.wmm)

  return result

# em initialization step
def initialize(sequence, k):
  i = 0
  increment = int(k / 2)
  init_pseudo_count = np.full((4), ((1 / c.seed_proportion) - 1) / 4)

  result = []

  while i + k <= len(sequence):
    window = [sequence[i : i + k]]

    count_matrix = m.makeCountMatrix(window)
    count_matrix = m.addPseudo(count_matrix, init_pseudo_count)

    freq_matrix = m.makeFrequencyMatrix(count_matrix)
    wmm_struct = m.makeWMM(freq_matrix, c.bg_vector)
    
    result.append(wmm_struct)
    i += increment

  # check if we missed a possible spot, and slide back accordingly
  diff = (i - len(sequence))
  if diff % increment != 0:
    i = len(sequence) - k
    window = [sequence[i : i + k]]

    count_matrix = m.makeCountMatrix(window)
    count_matrix = m.addPseudo(count_matrix, init_pseudo_count)

    freq_matrix = m.makeFrequencyMatrix(count_matrix)
    wmm_struct = m.makeWMM(freq_matrix, c.bg_vector)

    result.append(wmm_struct)

  return result
