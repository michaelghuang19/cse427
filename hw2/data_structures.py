# describes any custom data structures that we might need

"""
fasta info struct, with description and sequence fields
"""
class fasta_info:
  def __init__(self, description, sequence):
    self.description = description
    self.sequence = sequence

class score_info:
  def __init__(self, score, true):
    self.score = score
    self.true = true
