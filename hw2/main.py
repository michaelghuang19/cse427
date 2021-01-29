# Michael Huang (mhuang19)
# 1862567

import itertools
import pandas as pd
import numpy as np
import urllib as ul
from tabulate import tabulate

import constants as c
import data_structures as ds
import helper as h

# build a count matrix corresponding to a given list of length k sequences.
def makeCountMatrix(): 
  print("makeCountMatrix")

# given a count matrix, and a pseudocount vector, build a new count matrix by adding them.
def addPseudo(): 
  print("addPseudo")

# make a frequency matrix from a count matrix.
def makeFrequencyMatrix(): 
  print("makeFrequencyMatrix")

# calculate the entropy of a frequency matrix relative to background.
def entropy(): 
  print("entropy")

# make a weight matrix from a frequency matrix and background vector.
# (It may be convenient to save the entropy with the WMM, too.)
def makeWMM(): 
  print("makeWMM")

# given a WMM of width k and one or more sequences of varying lengths ≥ k,
# scan/score each position of each sequence(excluding those < k from the rightmost end).
def scanWMM(): 
  print("scanWMM")

# given an(estimated) WMM and a set of sequences, run the E-step of MEME's EM algorithm;
# i.e., what is E[zij], where zij is the zero-one variable indicating whether the
# motif instance in sequence i begins in position j. (scanWMM does most of the required work.)
def Estep(): 
  print("Estep")

# given the Estep result, re-estimate the WMM. 
# (This is similar to makeCountMatrix / addPseudo / makeFrequencyMatrix / makeWMM,
# with the additional wrinkle that the k-mer inputs to makeCountMatrix each have weights,
# and it may be convenient to generalize makeCountMatrix so that its k-mer inputs are successive, 
# overlapping subsequences of some longer strings, 
# rather than having them explicitly presented as distinct, non-overlapping inputs.)
def Mstep(): 
  print("Mstep")

def main():
  print("hello world")

  # 2. Process fasta into individual lists. This results in a list of lists,
  # where each list contains corresponding fasta data structures for each
  # individual accession file, should a file contain multiple sequences.
  # You shouldn't really be commenting this out in normal cases.
  fasta_files = []
  for accession in c.file_dict.keys():
    fasta_files.append(h.process_fasta(accession))

  for fasta_list in fasta_files:
    for fasta in fasta_list:
      # if (fasta.description == "mm9_chr1:80073137-80073149(+)"):
      print(fasta.sequence)

if __name__ == "__main__":
  main()
