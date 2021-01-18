import pandas as pd
import numpy as np
import csv
import urllib as ul
# import Bio

import fasta
import constants

def process_fasta(fasta):
  print(fasta)

def get_fasta():
  for item in constants.accession_set:
    url = constants.uniprot_url + item + "." + constants.fasta_const
    f = ul.request.urlopen(url)
    fasta = f.read()
    
    output = open(constants.fasta_const + "/" + item + "." +
    constants.fasta_const, 'wb')
    output.write(fasta)
    output.close()

def main():
  print ("hello world!")
  get_fasta()

if __name__ == "__main__":
  main()
