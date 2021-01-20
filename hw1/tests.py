import pandas as pd
import numpy as np
import urllib as ul

import main
import fasta
import constants

def test1():
  main.compare_seqs([fasta.fasta_info("", "s1", "", "DEADLY")],
                    [fasta.fasta_info("", "s2", "", "DDGEARLYK")],
                    999)

def test2():
  print("")

def test3():
  print("")

def site_tests():
  main.compare_seqs([fasta.fasta_info("", "seqid001", "", "KEVLAR")],
               [fasta.fasta_info("", "seqid002", "", "KNIEVIL")],
               main.num_epochs)
  main.compare_seqs([fasta.fasta_info("", "seqid003", "", "MELLSLCSWFAAATTYDADFYDDP")],
               [fasta.fasta_info("", "seqid004", "", "MSNWTATSSDSTS")],
               main.num_epochs)
