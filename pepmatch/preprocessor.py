import pickle
import sqlite3
import re
import os
import argparse

from .parser import parse_fasta


class Preprocessor(object):
  '''
  Object class that takes in a proteome FASTA file, k for k-mer size (split), and format to
  store preprocessed data.

  With the preprocess method, it will break the proteome into equal size k-mers and map
  them to locations within the individual proteins. The mapped keys and values will then
  be stored in either pickle files or a SQLite database.

  The proteins within the proteome will be assigned numbers along with the index position of
  each k-mer within the protein.

  Optional:
  preprocessed_files_path - default is current directory. This is the location where either
  the .db file will be stored (for exact matching) or the .pickle files will be stored (for
  mismatching).

  gene_priority_proteome - used by UniProt proteomes which contain a prioritized protein
  sequence per gene. This will be used in the search and prioritize these sequences when
  one match is requested.

  versioned_ids - protein IDs can be versioned, so the versioned_ids argument can be passed
  as True to store them as versioned.
  '''
  def __init__(self,
               proteome,
               preprocess_format,
               preprocessed_files_path='.',
               versioned_ids = True,
               gene_priority_proteome=''):

    if not preprocess_format in ('sql', 'pickle'):
      raise AssertionError('Unexpected value of preprocessing format:', preprocess_format)

    if not os.path.isdir(preprocessed_files_path):
      raise ValueError('Directory specified does not exist: ', preprocessed_files_path)

    self.proteome = proteome
    self.proteome_name = proteome.split('/')[-1].split('.')[0]
    self.preprocess_format = preprocess_format
    self.preprocessed_files_path = preprocessed_files_path
    self.gene_priority_proteome = gene_priority_proteome
    self.versioned_ids = versioned_ids

  def split_protein(self, seq, k):
    '''
    Splits a protein into equal sized k-mers on a rolling basis.
    Ex: k = 4, NSLFLTDLY --> ['NSLF', 'SLFL', 'LFLT', 'FLTD', 'LTDL', 'TDLY']
    '''
    kmers = []
    for i in range(len(seq)-k + 1):
      kmers.append(seq[i:i+k])
    return kmers

  def pickle_proteome(self, kmer_dict, names_dict, k):
    '''
    Takes the preprocessed proteome (below) and creates a pickle file for
    both k-mer and names dictionaries created. This is for compression and
    for being able to load the data in when a query is called.
    '''
    with open(os.path.join(self.preprocessed_files_path, 
      f'{self.proteome_name}_{str(k)}mers.pickle'), 'wb') as f:

      pickle.dump(kmer_dict, f)

    with open(os.path.join(self.preprocessed_files_path, 
      f'{self.proteome_name}_names.pickle'), 'wb') as f:

      pickle.dump(names_dict, f)

  def sql_proteome(self, kmer_dict, names_dict, k):
    '''
    Takes the preprocessed proteome (below) and creates SQLite tables for both the
    k-mer and names dictionaries created. These SQLite tables can then be used
    for searching. This is much faster for exact matching.
    '''
    kmers_table = f'{self.proteome_name}_{str(k)}mers'
    names_table = f'{self.proteome_name}_names'

    conn = sqlite3.connect(os.path.join(self.preprocessed_files_path, self.proteome_name + '.db'))
    c = conn.cursor()

    c.execute(f'CREATE TABLE IF NOT EXISTS "{kmers_table}"(kmer TEXT, position INT)')
    c.execute(f'CREATE TABLE IF NOT EXISTS "{names_table}"(protein_number INT, taxon INT, species TEXT, gene TEXT, protein_id TEXT, protein_name TEXT, pe_level INT, gene_priority INT)')

    # make a row for each unique k-mer and position mapping
    for kmer, positions in kmer_dict.items():
      for position in positions:
        c.execute(f'INSERT INTO "{kmers_table}" (kmer, position) VALUES (?, ?)', (str(kmer), position,))

    # make a row for each number to protein ID mapping
    for protein_number, protein_data in names_dict.items():
      c.execute(f'INSERT INTO "{names_table}"(protein_number, taxon, species, gene, protein_id, protein_name, pe_level, gene_priority) VALUES(?, ?, ?, ?, ?, ?, ?, ?)',
        (protein_number, protein_data[0], protein_data[1], protein_data[2], protein_data[3], protein_data[4], protein_data[5], protein_data[6]))

    # create indexes for both k-mer, unique position, and name tables
    c.execute(f'CREATE INDEX IF NOT EXISTS "{kmers_table}_kmer_id" ON "{kmers_table}"(kmer)')
    c.execute(f'CREATE INDEX IF NOT EXISTS "{names_table}_id" ON "{names_table}"(protein_number)')

    conn.commit()
    c.close()
    conn.close()

  def preprocess(self, k):
    '''
    Method which preprocessed the given proteome, by splitting each protein into k-mers
    and assigninga unique index to each unique k-mer within each protein. This is done by
    assigning a number to each protein and for each k-mer, multiplying the protein number
    by 100,000 and adding the index position of the index within the protein. This
    guarantees a unique index for each and every possible k-mer. Also, each protein #
    assigned is also mappedto the protein ID to be read back later after searching.
    '''
    assert k >= 2, 'k-sized split is invalid. Cannot be less than 2.'

    proteome = parse_fasta(self.proteome)
    kmer_dict = {}
    names_dict = {}
    protein_count = 1

    # get all gene priority IDs to check for them
    if self.gene_priority_proteome:
      gene_priority_proteome_ids = []
      gene_priority_proteome = parse_fasta(self.gene_priority_proteome)
      for protein in gene_priority_proteome:
        protein_id = protein.id.split('|')[1]
        gene_priority_proteome_ids.append(protein_id)

    for protein in proteome:
      kmers = self.split_protein(str(protein.seq), k)
      for i in range(len(kmers)):
        if kmers[i] in kmer_dict.keys():
          kmer_dict[kmers[i]].append(protein_count * 100000 + i) # add index to k-mer list
        else:
          kmer_dict[kmers[i]] = [protein_count * 100000 + i]     # create entry for new k-mer

      # grab UniProt ID which is usually in the middle of two vetical bars
      try:
        protein_id = protein.id.split('|')[1]
      except IndexError:
        protein_id = protein.id

      # use regex to get all data from the UniProt FASTA header
      try:
        taxon = int(re.search('OX=(.*?) ', protein.description).group(1))
      except AttributeError:
        taxon = None

      try:
        species = re.search('OS=(.*) OX=', protein.description).group(1)
      except AttributeError:
        species = None

      try:
        gene = re.search('GN=(.*?) ', protein.description).group(1)
      except AttributeError:
        gene = None

      try:
        protein_name = re.search(' (.*) OS', protein.description).group(1)
      except AttributeError:
        protein_name = None

      try:
        pe_level = int(re.search('PE=(.*?) ', protein.description).group(1))
      except AttributeError:
        pe_level = 0

      # label protein record as being in gene priority proteome if it's found there
      if self.gene_priority_proteome:
        gene_priority = 1 if protein_id in gene_priority_proteome_ids else 0
      else:
        gene_priority = 0

      if self.versioned_ids:
        try:
          versioned_id = re.search('SV=(.*)', protein.description).group(1)
          protein_id += '.' + versioned_id
        except AttributeError:
          versioned_id = None

      names_dict[protein_count] = (taxon, species, gene, protein_id, protein_name, pe_level, gene_priority)
      protein_count += 1

    # store data based on format specified
    if self.preprocess_format == 'pickle':
      self.pickle_proteome(kmer_dict, names_dict, k)
    elif self.preprocess_format == 'sql':
      self.sql_proteome(kmer_dict, names_dict, k)

    return 0


# run via command line

def parse_arguments():
  parser = argparse.ArgumentParser()

  parser.add_argument('-p', '--proteome', required=True)
  parser.add_argument('-k', '--kmer_size', type=int, required=True)
  parser.add_argument('-f', '--format', required=True)
  parser.add_argument('-P', '--preprocessed_files_path', default='.')
  parser.add_argument('-v', '--versioned_ids', type=bool, default=True)
  parser.add_argument('-g', '--gene_priority_proteome', default='')

  args = parser.parse_args()

  return args

def run():
  args = parse_arguments()

  Preprocessor(args.proteome, args.format, args.preprocessed_files_path,
               args.versioned_ids, args.gene_priority_proteome).preprocess(args.kmer_size)
