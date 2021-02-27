# Michael Huang (mhuang19)
# 1862567

import itertools as it
import math
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn import metrics
from tabulate import tabulate

import constants as c
import data_structures as ds
import helper as h



def main():
  print("hello world")

  fasta_list = h.process_fasta(c.genome_file, c.fna_exten)
  seq = h.get_seq_list(fasta_list)[0]

  ginfo_list = h.process_gff(c.genome_file, c.gff_exten)
  ginfo_list = sorted(ginfo_list, key=lambda ginfo: ginfo.start)

  long_orf_output = open(c.results_folder + "long_orf" + c.text_exten, "wt")
  
  long_orf_output.close()

  # known_orf_output = open(c.results_folder + "known_orf" + c.text_exten, "wt")

  # known_orf_output.close()

  print("done")

if __name__ == "__main__":
  main()
