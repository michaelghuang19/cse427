# describes the fasta info data structure
class fasta_info:
  def __init__(self, species, accession, description, sequence):
    self.species = species
    self.accession = accession
    self.description = description
    self.sequence = sequence

class align_info:
  def __init__(self, best, bestx, besty):
    self.best_score = best
    self.best_coord1 = bestx
    self.best_coord2 = besty
