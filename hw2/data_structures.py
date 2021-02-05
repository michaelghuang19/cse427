# describes any custom data structures that we might need


"""
fasta info struct, with description and sequence fields
"""
class fasta_info:
  def __init__(self, description, sequence):
    self.description = description
    self.sequence = sequence

"""
wmm info struct, with wmm and entropy fields
"""
class wmm_info:
  def __init__(self, wmm, entropy):
    self.wmm = wmm
    self.entropy = entropy
