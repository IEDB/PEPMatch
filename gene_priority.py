#!/usr/bin/env python3

from Bio import SeqIO


proteome = '/home/dan/Desktop/9606.fasta'
gene_priority_proteome = '/home/dan/Desktop/9606_small.fasta'

def append_gp(proteome, gene_priority_proteome):
  gpp = list(SeqIO.parse(gene_priority_proteome, 'fasta'))
  p = list(SeqIO.parse(proteome, 'fasta'))
  ids = []
  for record in gpp:
    ids.append(record.id)
  
  with open(proteome, 'w') as fout:
    for record in p:
      if record.id in ids:
        GP = 'GP=1'
      else:
        GP = 'GP=0'
      record.description += ' ' + GP
      SeqIO.write(record, fout, 'fasta')

if __name__ == '__main__':
  append_gp(proteome, gene_priority_proteome)
