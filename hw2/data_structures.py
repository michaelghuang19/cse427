# describes the fasta info data structure, and best alignment data structure
class fasta_info:
  def __init__(self, species, accession, description, sequence):
    self.species = species
    self.accession = accession
    self.description = description
    self.sequence = sequence
