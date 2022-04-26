import pickle
import sqlite3
import re
import os

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
               split,
               preprocess_format,
               preprocessed_files_path='.',
               gene_priority_proteome='',
               versioned_ids = True):

    assert split >= 2, 'k-sized split is invalid. Cannot be less than 2.'

    if not preprocess_format in ('sql', 'pickle'):
      raise AssertionError('Unexpected value of preprocessing format:', preprocess_format)

    if not os.path.isdir(preprocessed_files_path):
      raise ValueError('Directory specified does not exist: ', preprocessed_files_path)

    self.proteome = proteome
    self.proteome_name = proteome.split('/')[-1].split('.')[0]
    self.split = split
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

  def pickle_proteome(self, kmer_dict, names_dict):
    '''
    Takes the preprocessed proteome (below) and creates a pickle file for
    both k-mer and names dictionaries created. This is for compression and
    for being able to load the data in when a query is called.
    '''
    with open(os.path.join(self.preprocessed_files_path, self.proteome_name + '_' +
              str(self.split) + 'mers.pickle'), 'wb') as f:

      pickle.dump(kmer_dict, f)

    with open(os.path.join(self.preprocessed_files_path, self.proteome_name +
              '_names.pickle'), 'wb') as f:

      pickle.dump(names_dict, f)

  def sql_proteome(self, kmer_dict, names_dict):
    '''
    Takes the preprocessed proteome (below) and creates SQLite tables for both the
    k-mer and names dictionaries created. These SQLite tables can then be used
    for searching. This is much faster for exact matching.
    '''
    kmers_table = self.proteome_name + '_' + str(self.split) + 'mers'
    names_table = self.proteome_name + '_names'

    conn = sqlite3.connect(os.path.join(self.preprocessed_files_path, self.proteome_name + '.db'))
    c = conn.cursor()

    c.execute('CREATE TABLE IF NOT EXISTS "{k}"(kmer TEXT, position INT)'.format(k = kmers_table))
    c.execute('CREATE TABLE IF NOT EXISTS "{n}"(protein_number INT, taxon INT, species TEXT, gene TEXT, protein_id TEXT, protein_name TEXT, pe_level INT, gene_priority INT)'.format(n = names_table))

    # make a row for each unique k-mer and position mapping
    for kmer, positions in kmer_dict.items():
      for position in positions:
        c.execute('INSERT INTO "{k}" (kmer, position) VALUES (?, ?)'.format(k = kmers_table), (str(kmer), position,))

    # make a row for each number to protein ID mapping
    for protein_number, protein_data in names_dict.items():
      c.execute('INSERT INTO "{n}"(protein_number, taxon, species, gene, protein_id, protein_name, pe_level, gene_priority) VALUES(?, ?, ?, ?, ?, ?, ?, ?)'.format(n = names_table),
        (protein_number, protein_data[0], protein_data[1], protein_data[2], protein_data[3], protein_data[4], protein_data[5], protein_data[6]))

    # create indexes for both k-mer, unique position, and name tables
    c.execute('CREATE INDEX IF NOT EXISTS "{id}" ON "{k}"(kmer)'.format(id = kmers_table + '_kmer_id', k = kmers_table))
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

    if self.gene_priority_proteome:
      gene_priority_proteome_ids = []
      gene_priority_proteome = parse_fasta(self.gene_priority_proteome)
      for protein in gene_priority_proteome:
        protein_id = protein.id.split('|')[1]
        gene_priority_proteome_ids.append(protein_id)

    for protein in proteome:
      kmers = self.split_protein(str(protein.seq), self.split)
      for i in range(len(kmers)):
        if kmers[i] in kmer_dict.keys():
          kmer_dict[kmers[i]].append(protein_count * 100000 + i) # add index to k-mer list
        else:
          kmer_dict[kmers[i]] = [protein_count * 100000 + i]     # create entry for new k-mer

      if self.gene_priority_proteome:
        gene_priority = 1 if protein_id in gene_priority_proteome_ids else 0
      else:
        gene_priority = None

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
        pe_level = None

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
      self.pickle_proteome(kmer_dict, names_dict)
    elif self.preprocess_format == 'sql':
      self.sql_proteome(kmer_dict, names_dict)

    return kmer_dict, names_dict
