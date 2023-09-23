#!/usr/bin/env python3

import os
import pandas as pd
from Bio import SeqIO


directory = os.path.dirname(os.path.abspath(__file__))


def parse_fasta(file):
  return SeqIO.parse(file, 'fasta')


class DIAMOND(object):
  def __init__(self, query, proteome, max_mismatches, method_parameters):
    if max_mismatches == -1:
      max_mismatches = 7

    self.query = query
    self.proteome = proteome
    self.proteome_name = str(proteome).replace('.fasta', '')

    self.max_mismatches = max_mismatches

    bin_directory = method_parameters['bin_directory']
    self.bin_file = os.path.join(bin_directory, 'diamond')


  def __str__(self):
    return 'DIAMOND'
  

  def preprocesss(self):
    os.system(f"{self.bin_file} makedb --in {self.proteome} -d {self.proteome_name}")
  

  def diamond_search(self):
    os.system(
      f"{self.bin_file} blastp -d {self.proteome_name} -q {self.query} -o matches.m8 "
      f"-e 10000 -k 100 --ultra-sensitive --masking 0 -f 6 "
       "full_qseq sseq sseqid mismatch sstart"
    )

    all_matches = []
    with open('matches.m8', 'r') as file:
      lines = file.readlines()
      for line in lines:
        match = []
        result = line.split('\t')
        for i in range(len(result)):
          if i == 3: 
            continue  # skip the max_mismatches column
          if i == 4:
            match.append(int(result[i].replace('\n', '')))
          else:
            match.append(result[i])

        all_matches.append(match)

    return all_matches


class Benchmarker(DIAMOND):
  def __init__(
    self, benchmark: str, query: str, proteome: str, lengths: list, max_mismatches: int,
    method_parameters: dict
  ):
    self.benchmark = benchmark
    self.query = query
    self.proteome = proteome
    self.lengths = lengths
    self.max_mismatches = max_mismatches
    self.method_parameters = method_parameters
    
    super().__init__(query, proteome, max_mismatches, method_parameters)


  def __str__(self):
    return 'DIAMOND'


  def preprocess_proteome(self):
    return self.preprocesss()


  def preprocess_query(self):
    raise TypeError(self.__str__() + ' does not preprocess queries.\n')


  def search(self):
    matches = self.diamond_search()

    all_matches = []
    for match in matches:
      match = list(match)
      try: # get the UniProt ID or do nothing 
        match[2] = match[2].split('|')[1]
      except IndexError:
        pass
      all_matches.append(','.join([str(i) for i in match]))
    
    all_matches = [match.split(',') for match in all_matches]
    columns = ['Query Sequence', 'Matched Sequence', 'Protein ID', 'Index start']
    
    os.remove('matches.m8')
    os.remove(
      os.path.dirname(self.proteome) + '/%s.dmnd' % self.proteome_name.split('/')[-1]
    )

    return pd.DataFrame(all_matches, columns = columns)