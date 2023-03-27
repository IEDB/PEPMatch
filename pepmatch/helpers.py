#!/usr/bin/env python3

from Bio import SeqIO


def split_sequence(sequence, k):
  """
  Splits a peptide into equal sized k-mers on a rolling basis.
  Ex: k = 4, NSLFLTDLY --> ['NSLF', 'SLFL', 'LFLT', 'FLTD', 'LTDL', 'TDLY']
  """
  return [sequence[i:i+k] for i in range(len(sequence) - k + 1)]

def parse_fasta(file):
  """Return a parsed Biopython SeqRecord object from a FASTA file."""
  return list(SeqIO.parse(file, 'fasta'))