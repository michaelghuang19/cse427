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
  def __init__(self, sequence):
    self.sequence = sequence
    self.stop_list = np.array([])
    # be careful! these are formatted as [first nuc, last nuc]
    # i.e. you CAN'T directly use these to find sequence
    self.orc_locs = []
    self.orf_list = []
  
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

  def find_stops(self, offset):
    result = []

    # TODO: check if we can perform indexing better here
    for i in range(offset, len(self.sequence) - 1, 3):
      codon = self.sequence[i : i + 3]

      if (codon in c.stop_codons):
        result.append(i)
        
    result = np.array(result)

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
