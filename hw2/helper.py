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
    species, accession, description = info[0].split("|")
    sequence = ""
    for j in range(1, len(info)):
      sequence += info[j].upper()

    fasta_struct = ds.fasta_info(species, accession, description, sequence)
    fasta_array.append(fasta_struct)

  print("Successfully processed fasta for " + filename)

  return fasta_array