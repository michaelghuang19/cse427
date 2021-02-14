# Michael Huang (mhuang19)
# 1862567

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

  fasta_list = h.process_fasta(c.genome_file + c.fna_exten)
  seq = h.get_seq_list(fasta_list)[0]



  print("done")

if __name__ == "__main__":
  main()
