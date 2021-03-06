# describes any custom data structures that we might need

"""
fasta info struct, with description and sequence fields
"""
class fasta_info:
  def __init__(self, description, sequence):
    self.description = description
    self.sequence = sequence

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
