# describes any custom data structures that we might need

class fasta_info:
  def __init__(self, description, sequence):
    self.description = description
    self.sequence = sequence

class wmm_info:
  def __init__(self, entropy, wmm):
    self.entropy = entropy
    self.wmm = wmm
