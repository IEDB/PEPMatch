#!/usr/bin/env python3

import _pickle as pickle
import sqlite3

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

  Optional: protein IDs can be versioned, so the versioned_protein_ids argument can be passed
  as True to store them as versioned.
  '''
  def __init__(self, proteome, split, preprocess_format, database='', versioned_protein_ids = False):
    if split < 2:
      raise ValueError('k-sized split is invalid. Cannot be less than 2.')

    if preprocess_format == 'sql' and database == '':
      raise ValueError('SQL format selected but database path not specified.')

    self.proteome = proteome
    self.split = split
    self.preprocess_format = preprocess_format
    self.database = database
    self.versioned_protein_ids = versioned_protein_ids

  def split_protein(self, seq, k):
    '''
    Splits a protein into equal sized k-mers on a rolling basis.
    Ex: k = 4, NSLFLTDLY --> ['NSLF', 'SLFL', 'LFLT', 'FLTD', 'LTDL', 'TDLY']
    '''
    kmers = []
    for i in range(len(seq)-k + 1):
      kmer = seq[i:i+k]
      kmers.append(kmer)
    return kmers

  def pickle_proteome(self, kmer_dict, names_dict):
    '''
    Takes the preprocessed proteome (below) and creates a pickle file for 
    both k-mer and names dictionaries created. This is for compression and
    for being able to load the data in when a query is called.
    '''
    name = self.proteome.split('/')[-1].split('.')[0]
    with open(name + '_kmers' + '_' + str(self.split) + '.pickle', 'wb') as f:
      pickle.dump(kmer_dict, f)
    with open(name + '_names.pickle', 'wb') as f:
      pickle.dump(names_dict, f)

  def sql_proteome(self, kmer_dict, names_dict):
    '''
    Takes the preprocessed proteome (below) and creates SQLite tables for both the 
    k-mer and names dictionaries created. These SQLite tables can then be used 
    for searching. This is much faster for exact matching.
    '''
    name = self.proteome.split('/')[-1].split('.')[0]
    kmers_table = name + '_kmers' + '_' + str(self.split)
    names_table = name + '_names'

    conn = sqlite3.connect(self.database)
    c = conn.cursor()
    c.execute('CREATE TABLE IF NOT EXISTS "{k}"(kmer TEXT, position INT)'.format(k = kmers_table))
    c.execute('CREATE TABLE IF NOT EXISTS "{n}"(protein_number INT, protein_ID TEXT)'.format(n = names_table))

    # make a row for each unique k-mer and position mapping
    for kmer, positions in kmer_dict.items():
      for position in positions:
        c.execute('INSERT INTO "{k}" (kmer, position) VALUES (?, ?)'.format(k = kmers_table), (str(kmer), position,))

    # make a row for each number to protein ID mapping
    for protein_number, protein_ID in names_dict.items():
      c.execute('INSERT INTO "{n}"(protein_number, protein_ID) VALUES(?, ?)'.format(n = names_table), 
        (protein_number, protein_ID))

    # create indexes for both k-mer and name tables
    c.execute('CREATE INDEX IF NOT EXISTS "{id}" ON "{k}"(kmer)'.format(id = kmers_table + '_id', k = kmers_table))
    c.execute('CREATE INDEX IF NOT EXISTS "{id}" ON "{n}"(protein_number)'.format(id = names_table + '_id', n = names_table))
    conn.commit()
    c.close()
    conn.close()

  def preprocess(self):
    '''
    Method which preprocessed the given proteome, by splitting each protein into k-mers 
    and assigninga unique index to each unique k-mer within each protein. This is done by 
    assigning a number to each protein and for each k-mer, multiplying the protein number 
    by 100,000 and adding the index position of the index within the protein. This 
    guarantees a unique index for each and every possible k-mer. Also, each protein # 
    assigned is also mappedto the protein ID to be read back later after searching.
    '''
    proteome = parse_fasta(self.proteome)
    kmer_dict = {}
    names_dict = {}
    protein_count = 1

    for protein in proteome:
      kmers = self.split_protein(str(protein.seq), self.split)
      for i in range(len(kmers)):
        if kmers[i] in kmer_dict.keys():
          kmer_dict[kmers[i]].append(protein_count * 100000 + i) # add index to k-mer list 
        else: 
          kmer_dict[kmers[i]] = [protein_count * 100000 + i]     # create entry for new k-mer

      # create names mapping # to protein ID (include versioned if argument is passed) 
      try:
        if self.versioned_protein_ids:
          protein_ID = str(protein.description).split(' ')[0]
          names_dict[protein_count] = protein_ID.split('|',1)[0] + '|' + protein_ID.split('|')[1] + '.' + str(protein.description).split('SV=')[1] + '|' + protein_ID.split('|')[2]
        else:
          names_dict[protein_count] = str(protein.description).split(' ')[0]
      except IndexError:
        names_dict[protein_count] = str(protein.description).split(' ')[0]

      protein_count += 1

    if self.preprocess_format == 'pickle':
      self.pickle_proteome(kmer_dict, names_dict)
    elif self.preprocess_format == 'sql':
      self.sql_proteome(kmer_dict, names_dict)
    else:
      raise AssertionError('Unexpected value of preprocessing format', self.preprocess_format)
    return kmer_dict, names_dict