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
  seq = fasta_list[0].sequence

  # With the test, we expect the following output:
  # ATAACG, CCC
  # CGTGAC
  # AACGTGACC
  
  for i in range(3):
    long_orf_output = open(c.results_folder + "long_orf{}".format(i) + c.text_exten, "wt")
  
    long_orf = ds.orf_struct(seq)
    long_orf.find_orf_locs(i)

    for j in range(len(long_orf.orf_locs)):
      h.print_1list(long_orf.orf_locs[j], long_orf_output)
      long_orf_output.write("\t" + str(long_orf.orf_list[j]) + "\n")

    long_orf_output.close()

  # known_orf_output = open(c.results_folder + "known_orf" + c.text_exten, "wt")

  # known_orf_output.close()

  # ginfo_list = h.process_gff(c.genome_file, c.gff_exten)
  # ginfo_list = sorted(ginfo_list, key=lambda ginfo: ginfo.start)

  print("done")

if __name__ == "__main__":
  main()
