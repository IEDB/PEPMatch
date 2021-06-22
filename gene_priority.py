#!/usr/bin/env python3

import sys

from Bio import SeqIO

# takes in a proteome and a gene priority proteome from
# UniProt and adds GP as part of the FASTA headers
# to indentify protein records as beloning to the gene priority
def append_gp(proteome, gene_priority_proteome):
  gpp = list(SeqIO.parse(gene_priority_proteome, 'fasta'))
  p = list(SeqIO.parse(proteome, 'fasta'))
  ids = []
  for record in gpp:
    ids.append(record.id)
  
  with open(proteome, 'w') as fout:
    for record in p:
      if 'GP' in record.description:
        continue
      if record.id in ids:
        GP = 'GP=1'
      else:
        GP = 'GP=0'
      record.description += ' ' + GP
      SeqIO.write(record, fout, 'fasta')

if __name__ == '__main__':
  proteome = sys.argv[1]
  gene_priority_proteome = sys.argv[2]
  append_gp(proteome, gene_priority_proteome)
