# describes the fasta info data structure
class fasta_info:
  def __init__(self, species, accession, description, sequence):
    self.species = species
    self.accession = accession
    self.description = description
    self.sequence = sequence
  def add(self, n1, n2):
    return n1 + n2
  def sub(self, n1, n2):
    return n1 - n2
 
#Bonus Function       
def display():
    print("Hello World")
