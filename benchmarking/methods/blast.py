#!/usr/bin/env python3

import os 
import glob
import pandas as pd

from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastpCommandline


directory = os.path.dirname(os.path.abspath(__file__))


def parse_fasta(file):
  return SeqIO.parse(file, 'fasta')


class BLAST(object):
  def __init__(self, query, proteome, max_mismatches, method_parameters):
    if max_mismatches == -1:
      max_mismatches = 7

    self.query = query
    self.proteome = proteome

    self.max_mismatches = max_mismatches

    bin_directory = method_parameters['bin_directory']
    self.makeblastdb_path = os.path.join(bin_directory, 'makeblastdb')
    self.blastp_path = os.path.join(bin_directory, 'blastp')


  def __str__(self):
    return 'BLAST'
  

  def preprocess(self):
    os.system(self.makeblastdb_path + ' -in ' + str(self.proteome) + ' -dbtype prot')
  

  def blast_search(self):
    peptides = parse_fasta(self.query)
    proteins = parse_fasta(self.proteome)

    peptide_dict = {}
    for peptide in peptides:
      peptide_dict[str(peptide.id)] = str(peptide.seq)

    protein_dict = {}
    for protein in proteins:
      try:
        protein_dict[str(protein.id).split('|')[1]] = str(protein.seq)
      except IndexError:
        protein_dict[str(protein.id)] = str(protein.seq)

    if self.max_mismatches == 0:
      blastx_cline = NcbiblastpCommandline(
        cmd=self.blastp_path,
        query = self.query, 
        db = self.proteome, 
        evalue=100, outfmt=10, out='output.csv')

      stdout, stderr = blastx_cline()

    else:
      blastx_cline = NcbiblastpCommandline(
        cmd=self.blastp_path, 
        query = self.query, 
        db = self.proteome, 
        evalue=10000, outfmt=10, out='output.csv')

      stdout, stderr = blastx_cline()

    df = pd.read_csv(
      'output.csv', 
      names = [
        'Peptide Sequence', 'Protein ID', 'Sequence Identity', 
        'Length', 'Mismatches', 'Gap Openings', 'Query start', 
        'Query end', 'Index start', 'Index end', 'e value', 'bit score'
      ]
    )
    df['Peptide Sequence'] = df['Peptide Sequence'].apply(str)
    df = df.replace({'Peptide Sequence': peptide_dict})

    df.to_csv('blast_results.csv')

    all_matches = []
    for i, row in df.iterrows():
      try:
        # for UniProt IDs - betacoronaviruses have different NCBI IDs
        row['Protein ID'] = row['Protein ID'].split('|')[1]
      except IndexError:
        pass
      
      peptide_sequence = row['Peptide Sequence']
      protein_id = row['Protein ID']
      index_start = int(row['Index start']) - 1 # BLAST is 1-indexed
      
      all_matches.append((
        peptide_sequence, 
        protein_dict[protein_id][index_start:(index_start+len(peptide_sequence))],
        protein_id,
        index_start + 1,
      ))

    return all_matches


class Benchmarker(BLAST):
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
    return 'BLAST'


  def preprocess_proteome(self):
    return self.preprocess()


  def preprocess_query(self):
    raise TypeError(self.__str__() + ' does not preprocess queries.\n')


  def search(self):
    matches = self.blast_search()

    all_matches = []
    for match in matches:
      match = list(match)
      try: # try taking the UniProt ID - else do nothing 
        match[2] = match[2].split('|')[1]
      except IndexError:
        pass
      all_matches.append(','.join([str(i) for i in match]))

    all_matches = [match.split(',') for match in all_matches]
    columns = ['Query Sequence', 'Matched Sequence', 'Protein ID', 'Index start']

    for extension in ['pdb', 'phr', 'pin', 'psq', 'ptf', 'pot', 'pto']:
      os.remove(glob.glob(os.path.dirname(self.proteome) + '/*.' + extension)[0])

    os.remove('output.csv')
    os.remove('blast_results.csv')

    return pd.DataFrame(all_matches, columns = columns)