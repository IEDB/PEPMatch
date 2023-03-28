#!/usr/bin/env python3

import re

from Bio import SeqIO


def parse_fasta(file):
  """Return a parsed Biopython SeqRecord object from a FASTA file."""
  return list(SeqIO.parse(file, 'fasta'))

def split_sequence(sequence, k):
  """
  Splits a peptide into equal sized k-mers on a rolling basis.
  Ex: k = 4, NSLFLTDLY --> ['NSLF', 'SLFL', 'LFLT', 'FLTD', 'LTDL', 'TDLY']
  """
  return [sequence[i:i+k] for i in range(len(sequence) - k + 1)]

def extract_metadata(record):
  """Extract metadata from a FASTA header."""
  regexes = {
      'protein_id': re.compile(r"\|([^|]*)\|"),      # between | and |
      'protein_name': re.compile(r"\s(.+?)\sOS"),    # between first space and space before OS
      'species': re.compile(r"OS=(.+?)\sOX"),        # between OS= and space before OX (species can have spaces)
      'taxon_id': re.compile(r"OX=(.+?)\s"),         # between OX= and space
      'gene': re.compile(r"GN=(.+?)\s"),             # between GN= and space
      'pe_level': re.compile(r"PE=(.+?)\s"),         # between PE= and space
      'sequence_version': re.compile(r"SV=(.+?)\s"), # between SV= and space
      'gene_priority': re.compile(r"GP=(.+?)\s"),    # between GP= and space
  }
  metadata = []
  for key in regexes: # loop through compiled regexes to extract metadata
    match = regexes[key].search(str(record.description))
    
    if match:
      metadata.append(match.group(1))
    else:
      if key == 'protein_id':
        metadata.append(str(record.id)) # get record.id from FASTA header instead
      elif key == 'sequence_version':
        metadata.append('1')
      elif key in ['pe_level', 'gene_priority']:
        metadata.append('0') # zeros for integer columns
      else:
        metadata.append('')  # empty strings for string columns

  return metadata