import pandas as pd
import numpy as np
import urllib as ul

import main
import fasta
import constants

def test1(output):
  main.compare_seqs([fasta.fasta_info("", "s1", "", "DEADLY")],
                    [fasta.fasta_info("", "s2", "", "DDGEARLYK")],
                    999,
                    output)

def site_tests(output):
  main.compare_seqs([fasta.fasta_info("", "seqid001", "", "KEVLAR")],
               [fasta.fasta_info("", "seqid002", "", "KNIEVIL")],
               main.num_epochs,
               output)
  main.compare_seqs([fasta.fasta_info("", "seqid003", "", "MELLSLCSWFAAATTYDADFYDDP")],
               [fasta.fasta_info("", "seqid004", "", "MSNWTATSSDSTS")],
               main.num_epochs,
               output)
