#!/usr/bin/env python3

import pandas as pd
from Bio import SeqIO


def parse_fasta(file):
  return SeqIO.parse(file, 'fasta')


class KnuthMorrisPratt(object):
  def __init__(self, query, proteome):
    self.query = query
    self.proteome = proteome

  def build_lps(self, pattern):
    m = len(pattern)
    lps_array = [0] * m
    i, j = 1, 0  # start from the 2nd character in pattern
    while i < m:
      if pattern[i] == pattern[j]:
        lps_array[i] = j + 1
        j += 1
        i += 1
      else:
        if j > 0:
          j = lps_array[j - 1]
        else:
          lps_array[i] = 0
          i += 1
    return lps_array


  def exact_search(self):
    query = list(parse_fasta(self.query))
    proteome = list(parse_fasta(self.proteome))

    all_matches = []
    for peptide in query:
      for protein in proteome:

        pattern = str(peptide.seq)
        text = str(protein.seq)

        if not text and not pattern:
          return 0
        elif not pattern:
          return 0

        # build longest proper suffix array for pattern
        lps_array = self.build_lps(pattern)

        n, m = len(text), len(pattern)
        i, j = 0, 0
        while i < n:
          # current characters match, move to the next characters

          if text[i] == pattern[j]:
            i += 1
            j += 1
          # current characters don't match
          else:
            if j > 0:  # try start with previous longest prefix
              j = lps_array[j - 1]
            # 1st character of pattern doesn't match character in text
            # go to the next character in text
            else:
              i += 1

          if j == m:
            all_matches.append((
              str(peptide.seq), str(peptide.seq), str(protein.id), i - m + 1
            ))
            i += m + 1
            j = lps_array[j - 1]


    return all_matches


class Benchmarker(KnuthMorrisPratt):
  def __init__(
    self, benchmark: str, query: str, proteome: str, lengths: list, max_mismatches: int,
    method_parameters: dict
  ):    
    if max_mismatches > 0:
      raise ValueError(self.__str__() + ' cannot do any mismatching.\n')
    elif max_mismatches == -1:
      raise ValueError(self.__str__() + ' does not have a best match feature.\n')

    self.benchmark = benchmark
    self.query = query
    self.proteome = proteome
    self.lengths = lengths
    self.max_mismatches = max_mismatches
    self.method_parameters = method_parameters
    
    super().__init__(query, proteome)


  def __str__(self):
    return 'Knuth-Morris-Pratt algorithm'


  def preprocess_proteome(self):
    raise TypeError(self.__str__() + ' does not preprocess proteomes.\n')


  def preprocess_query(self):
    raise TypeError(self.__str__() + ' does not preprocess queries.\n')


  def search(self):
    matches = self.exact_search()

    all_matches = []
    for match in matches:
      match = list(match)
      match[2] = match[2].split('|')[1]
      all_matches.append(','.join([str(i) for i in match]))

    all_matches = [match.split(',') for match in all_matches]
    columns = ['Query Sequence', 'Matched Sequence', 'Protein ID', 'Index start']
    
    return pd.DataFrame(all_matches, columns = columns)