class fasta:
  def __init__(self, species, accession, description):
    self.species = species
    self.accession = accession
    self.description = description
  def add(self, n1, n2):
    return n1 + n2
  def sub(self, n1, n2):
    return n1 - n2
 
#Bonus Function       
def display():
    print("Hello World")