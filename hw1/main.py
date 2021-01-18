import pandas as pd
import numpy as np
import csv
import urllib

# import Bio

import fasta

name_set = {"MYOD1_HUMAN", "TAL1_HUMAN", "MYOD1_MOUSE", "MYOD1_CHICK", 
"MYODA_XENLA", "MYOD1_DANRE", "Q8IU24_BRABE", "MYOD_DROME", "LIN32_CAEEL", 
"SYFM_HUMAN"}

accession_set = {"P15172", "P17542", "P10085", "P16075", "P13904", "Q90477", 
"Q8IU24", "P22816", "Q10574", "O95363"}

def process_fasta(fasta):
  print(fasta)

def get_fasta():
  uniprot_url = "www.uniprot.org/uniprot/"
  for item in accession_set:
    print(uniprot_url + item + ".fasta");

def main():
  print ("hello world!")
  get_fasta()

if __name__ == "__main__":
  main()