# Various helper functions
import constants as c
import data_structures as ds

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
    if (char.upper() == "A" or char.upper() == "C"
        or char.upper() == "G" or char.upper() == "T"):
      result += char.upper()
    elif (not char.isspace()):
      result += "T"

  return result


def create_empty_count_matrix():
  print("create_empty_count_matrix")

def create_empty_frequency_matrix():
  print("create_empty_frequency_matrix")

def create_empty_weight_matrix():
  print("create_empty_weight_matrix")


