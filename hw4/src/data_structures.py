# describes any custom data structures that we might need

import numpy as np

import constants as c

"""
fasta info struct, with description and sequence fields
"""
class fasta_info:
  def __init__(self, description, sequence):
    self.description = description
    self.sequence = sequence

class orf_struct:
  def __init__(self, *args):
    if len(args) == 0:
      self.sequence = ""
    else:
      self.sequence = args[0]
    
    self.stop_list = []

    # be careful! these are formatted as [first nuc, last nuc]
    # i.e. you CAN'T directly use these to find sequence
    if len(args) == 3:
      self.orf_locs = args[1]
    else:
      self.orc_locs = []

    if len(args) == 3:
      self.orf_list = args[2]
    else:
      self.orf_list = []
    
    self.bg_seqs = []

  def find_orf_locs(self, offset):
    print("finding orfs with offset {}".format(offset))
    locs = []
    orfs = []

    self.find_stops(offset)

    start = offset
    for end in self.stop_list:
      if start < end - 1:
        locs.append([start, end - 1])
        orfs.append(self.sequence[start : end])
      
      start = end + 3

    # see if there is a sequence left to process at the end
    last_len = len(self.sequence) - start
    if (last_len / 3) >= 1:
      last_len -= (last_len % 3)
      end = start + last_len
      locs.append([start, end - 1])
      orfs.append(self.sequence[start : end])

    self.orf_locs = locs
    self.orf_list = orfs
    return locs, orfs

  def find_bg_seqs(self):
    result = []

    for orf in self.orf_list:
      rev_orf = orf[::-1]
      rev_orf_complements = [""] * len(rev_orf)

      for i, nuc in enumerate(rev_orf):
        rev_orf_complements[i] = c.nuc_map[nuc]
      
      result.append(''.join(rev_orf_complements))

    self.bg_seqs = result
    return result

  def find_stops(self, offset):
    result = []

    # TODO: check if we can perform indexing better here
    for i in range(offset, len(self.sequence) - 1, 3):
      codon = self.sequence[i : i + 3]

      if (codon in c.stop_codons):
        result.append(i)

    self.stop_list = result
    return result

class gene_info:
  def __init__(self, seqid, ftype, start, end, strand, attributes):
    self.seqid = seqid
    self.ftype = ftype
    self.start = int(start)
    self.end = int(end)
    self.strand = strand
    self.attributes = attributes

  def __str__(self):
    result = "\t" + self.seqid + " " + self.ftype + " " + self.strand + "\n\t"
    result += str([self.start, self.end]) + "\n\t"
    attributes_list = self.attributes.split(";")
    for attribute in attributes_list:
      result += attribute + " "
    result += "\n"

    return result
