import pandas as pd
import numpy as np
import csv
import urllib as ul
# import Bio

import fasta
import constants

def process_fasta(fasta):
  print("we process a fasta")

# convert fasta data to data structure format
def create_fasta(filename):
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
    species,accession,description = info[0].split("|")
    sequence = ""
    for j in range(1, len(info)):
      sequence += info[j]

    fasta_struct = fasta.fasta_info(species, accession, description, sequence)
    fasta_array.append(fasta_struct)

  return fasta_array

# extra credit j: automated fasta retrieval + write into fasta folder
def get_fasta():
  for item in constants.accession_set:
    url = constants.uniprot_url + item + constants.fasta_exten
    f = ul.request.urlopen(url)
    fasta = f.read()
    
    output = open(constants.fasta_folder + item + constants.fasta_exten, "wb")
    output.write(fasta)
    output.close()

def main():
  print("hello world!")
  # get_fasta()

  fasta_list = []
  for accession in constants.accession_set:
    fasta_list.append(create_fasta(accession))
  
  for file in fasta_list:
    for seq in file:
      print(seq.species)
      print(seq.accession)
      print(seq.description)
      print(seq.sequence)

if __name__ == "__main__":
  main()
