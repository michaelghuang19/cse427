# various helper functions

import constants as c
import data_structures as ds

# convert gff data to a list of known genes
def process_gff(filename, exten):
  result = []

  file = open(c.data_folder + filename + exten, "r")
  text = file.read()
  gene_list = text.split("\n")

  assert len(gene_list) > 0

  gene_list = sorted(gene_list)
  while len(gene_list[0]) == 0 or gene_list[0][0] == "#":
    gene_list.pop(0)

  for gene in gene_list:
    col_vals = gene.split("\t")

    info = ds.gene_info(
        col_vals[0], col_vals[2], col_vals[3], col_vals[4], col_vals[6], col_vals[8])

    result.append(info)

  print("Successfully processed {} for ".format(exten) + filename)

  return result
