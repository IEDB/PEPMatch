#!/usr/bin/env python3

from Bio import SeqIO

def parse_fasta(file):
  '''Return a parsed Biopython SeqRecord object from a FASTA file.'''
  return SeqIO.parse(file, 'fasta')