# describes any custom data structures that we might need

"""
fasta info struct, with description and sequence fields
"""
class fasta_info:
  def __init__(self, description, sequence):
    self.description = description
    self.sequence = sequence

class gene_info:
  def __init__(self, seqid, source, ftype, start, end, strand, attributes):
    self.seqid = seqid
    self.source = source
    self.ftype = ftype
    self.start = start
    self.end = end
    self.strand = strand
