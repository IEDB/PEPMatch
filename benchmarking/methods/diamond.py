#!/usr/bin/env python3

from Bio import SeqIO
import os


directory = os.path.dirname(os.path.abspath(__file__))


def parse_fasta(file):
  return SeqIO.parse(file, 'fasta')


class DIAMOND(object):
  def __init__(self, query, proteome, max_mismatches, method_parameters):
    if max_mismatches == -1:
      max_mismatches = 7

    self.query = query
    self.proteome = proteome
    self.proteome_name = proteome.replace('.fasta', '')

    self.max_mismatches = max_mismatches

    bin_directory = method_parameters['bin_directory']

    self.bin_file = os.path.join(bin_directory, 'diamond')


  def __str__(self):
    return 'DIAMOND'
  

  def preprocesss(self):
    os.system(f"{self.bin_file} makedb --in {self.proteome} -d {self.proteome_name}")
  

  def diamond_search(self, query, proteome):
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
          if i == 4:
            match.append(int(result[i].replace('\n', '')) - 1)
          else:
            match.append(result[i])

        all_matches.append(match)

    return all_matches


class Benchmarker(DIAMOND):
  def __init__(self, query, proteome, lengths, max_mismatches, method_parameters):
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
    matches = self.diamond_search(self.query, self.proteome)

    all_matches = []
    for match in matches:
      match = list(match)
      try: # get the UniProt ID or do nothing 
        match[2] = match[2].split('|')[1]
      except IndexError:
        pass
      all_matches.append(','.join([str(i) for i in match]))

    os.remove('matches.m8')
    os.remove(
      os.path.dirname(self.proteome) + '/%s.dmnd' % self.proteome_name.split('/')[-1]
    )

    return all_matches