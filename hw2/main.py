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

def main():
  print("hello world")

  # 2. Process fasta into individual lists. This results in a list of lists,
  # where each list contains corresponding fasta data structures for each
  # individual accession file, should a file contain multiple sequences.
  # You shouldn't really be commenting this out in normal cases.
  fasta_list = []
  for accession in c.file_set:
    fasta_list.append(h.process_fasta(accession[0]))

if __name__ == "__main__":
  main()
