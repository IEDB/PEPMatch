import _pickle as pickle
import sqlite3
import os
import argparse

from .helpers import parse_fasta, split_sequence, extract_metadata


class Preprocessor:
  """
  Class that takes in a proteome FASTA file and preprocesses the file to be
  used in the matching portion of PEPMatch.

  Two tables are stored after the preprocessing:
  1. kmers - stores the k-mers and the index location of the k-mer within the
             entire proteome. The index is the protein_number * 1000000 + 
             the index within the protein.
  2. metadata - stores the metadata of the proteome. The metadata includes
                the protein ID, protein name, species, taxon ID, gene,
                protein existence level, sequence version, and gene priority
                (1 if the protein is in the gene priority proteome, 0 otherwise).

  The tables are stored in a SQLite database file or in pickle files.

  ** SQLite for exact matching and pickle files for mismatching. **

  Optional:
  preprocessed_files_path - location of where to put the preprocessed files.

  gene_priority_proteome - another proteome FASTA file that contains the proteins
                           the canonical sequence for every gene of the species
                           of the proteome. If this is specified, the GP=1 will
                           be appended to the FASTA header of the proteins in
                           this proteome. Otherwise, GP=0.
  """
  def __init__(self,
               proteome,
               preprocessed_files_path='.',
               gene_priority_proteome=''):

    if not os.path.isdir(preprocessed_files_path):
      raise ValueError('Directory specified does not exist: ', preprocessed_files_path)

    self.preprocessed_files_path = preprocessed_files_path

    # appending GP=# to the header of the proteome
    self.proteome = self._append_gp_to_header(proteome, gene_priority_proteome)

    # extracting the proteome name from the file path
    self.proteome_name = proteome.split('/')[-1].split('.')[0]

    # extract all the data from the proteome
    self.all_seqs, self.all_metadata = self._get_data_from_proteome()

  def _append_gp_to_header(self, proteome, gene_priority_proteome):
    """
    Appends GP=1 or GP=0 to the FASTA header of a proteome depending on
    if the protein is in the gene priority proteome.
    """
    proteome_records = parse_fasta(proteome) # get regular proteome
    try: # try to get gene priority proteome
      gp_proteome_records = parse_fasta(gene_priority_proteome)
    except FileNotFoundError:
      gp_proteome_records = []

    gp_ids = [] # get list of gene priority protein IDs
    for record in gp_proteome_records:
      gp_ids.append(record.id)
    
    records = []
    for record in proteome_records:
      if 'GP=' in record.description:
        records.append(record)
        continue
      
      if record.id in gp_ids:
        GP = 'GP=1 '
      else:
        GP = 'GP=0 '

      record.description += f' {GP}' # append notation to FASTA header
      records.append(record)

    return records

  def _get_data_from_proteome(self):
    """
    Extract all the data from a FASTA file and returns two lists:

    Use extract_metadata_from_fasta_header in helpers.py.

    1. A list of sequences
    2. A list of metadata in tuples. The metadata includes the protein ID,
       protein name, species, taxon ID, gene, protein existence level,
       sequence version, and gene priority label. 
    """
    all_seqs = []
    all_metadata = []
    protein_number = 1
    for record in self.proteome:
      all_seqs.append(str(record.seq))

      metadata = [protein_number]
      metadata.extend(extract_metadata(record))

      all_metadata.append(tuple(metadata))
      protein_number += 1

    return all_seqs, all_metadata

  def sql_proteome(self, k):
    """
    Writes the kmers_table and metadata_table to a SQLite database.
    """
    # create table names
    kmers_table = f'{self.proteome_name}_{str(k)}mers'
    metadata_table = f'{self.proteome_name}_metadata'
    
    # connect to database
    db_path = os.path.join(self.preprocessed_files_path, f'{self.proteome_name}.db')
    with sqlite3.connect(db_path) as conn:
      cursor = conn.cursor()

      # check if kmers table already exists and exit if it does
      # for some reason writing to the same table multiple times messes up results
      if cursor.execute(f'SELECT name FROM sqlite_master WHERE type="table" AND name="{kmers_table}";').fetchone():
        return

      self._create_tables(cursor, kmers_table, metadata_table)
      self._insert_kmers(cursor, kmers_table, k)
      self._insert_metadata(cursor, metadata_table)
      self._create_indexes(cursor, kmers_table, metadata_table)

      conn.commit()

  def _create_tables(self, cursor, kmers_table, metadata_table):
    cursor.execute(
      f'CREATE TABLE IF NOT EXISTS {kmers_table} ('\
        'kmer TEXT NOT NULL,'\
        'idx  INTEGER NOT NULL)'
      )
    cursor.execute(
      f'CREATE TABLE IF NOT EXISTS {metadata_table} ('\
        'protein_number   INTEGER NOT NULL,'\
        'protein_id       TEXT NOT NULL,'\
        'protein_name     TEXT NOT NULL,'\
        'species          TEXT NOT NULL,'\
        'taxon_id         INTEGER NOT NULL,'\
        'gene             TEXT NOT NULL,'\
        'pe_level         INTEGER NOT NULL,'\
        'sequence_version INTEGER NOT NULL,'\
        'gene_priority    INTEGER NOT NULL)'\
    )

  def _insert_kmers(self, cursor, kmers_table, k):
    kmer_rows = []
    for protein_count, seq in enumerate(self.all_seqs):
        for j, kmer in enumerate(split_sequence(seq, k)):
            kmer_rows.append((kmer, (protein_count + 1) * 1000000 + j))
    cursor.executemany(f'INSERT INTO {kmers_table} VALUES (?, ?)', kmer_rows)

  def _insert_metadata(self, cursor, metadata_table):
    cursor.executemany(f'INSERT INTO {metadata_table} VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)', self.all_metadata)

  def _create_indexes(self, cursor, kmers_table, metadata_table):
    cursor.execute(f'CREATE INDEX IF NOT EXISTS {kmers_table}_kmer_idx ON {kmers_table} (kmer)')
    cursor.execute(f'CREATE INDEX IF NOT EXISTS {metadata_table}_protein_number_idx ON {metadata_table} (protein_number)')

  def pickle_proteome(self, k):
    """
    Pickles the proteome into a dictionary of k-mers and a dictionary of metadata.
    """
    # create kmer_dict and metadata_dict out of self.all_seqs and self.all_metadata
    kmer_dict = {}
    for protein_count, seq in enumerate(self.all_seqs):
      for j, kmer in enumerate(split_sequence(seq, k)):
        if kmer in kmer_dict.keys():
          kmer_dict[kmer].append((protein_count + 1) * 1000000 + j) # add index to k-mer list
        else:
          kmer_dict[kmer] = [(protein_count + 1) * 1000000 + j] # create entry for new k-mer
    
    metadata_dict = {}
    for data in self.all_metadata:
      metadata_dict[data[0]] = data[1:]
    
    # write kmer_dict and metadata_dict to pickle files
    with open(os.path.join(self.preprocessed_files_path, 
      f'{self.proteome_name}_{str(k)}mers.pickle'), 'wb') as f:
      pickle.dump(kmer_dict, f)

    with open(os.path.join(self.preprocessed_files_path, 
      f'{self.proteome_name}_metadata.pickle'), 'wb') as f:
      pickle.dump(metadata_dict, f)

  def redis_proteome(self, k):
    """
    Writes the k-mers and metadata to a Redis database.
    """
    import redis
    r = redis.Redis(host='localhost', port=6379, db=0)

    # store k-mers
    for protein_count, seq in enumerate(self.all_seqs):
      for j, kmer in enumerate(split_sequence(seq, k)):
        idx = (protein_count + 1) * 1000000 + j
        key = f'kmer:{kmer}'
        r.set(key, str(idx))

    # store metadata
    for data in self.all_metadata:
      protein_number = data[0]
      metadata_values = ','.join([str(val) for val in data[1:]])
      key = f'metadata:{protein_number}'
      r.set(key, metadata_values)

  def preprocess(self, preprocess_format, k):
    """
    Preprocesses the proteome and stores it in the specified format.
    """
    if not preprocess_format in ('sql', 'pickle', 'redis'):
      raise AssertionError('Unexpected value of preprocessing format:', preprocess_format)

    assert k >= 2, 'k-sized split is invalid. Cannot be less than 2.'

    # store data based on format specified
    if preprocess_format == 'pickle':
      self.pickle_proteome(k)
    elif preprocess_format == 'sql':
      self.sql_proteome(k)
    elif preprocess_format == 'redis':
      self.redis_proteome(k)



# run via command line
def parse_arguments():
  parser = argparse.ArgumentParser()

  parser.add_argument('-p', '--proteome', required=True)
  parser.add_argument('-k', '--kmer_size', type=int, required=True)
  parser.add_argument('-f', '--preprocess_format', required=True)
  parser.add_argument('-P', '--preprocessed_files_path', default='.')
  parser.add_argument('-g', '--gene_priority_proteome', default='')

  args = parser.parse_args()

  return args

def run():
  args = parse_arguments()

  Preprocessor(
    args.proteome, 
    args.preprocessed_files_path,
    args.gene_priority_proteome).preprocess(args.preprocess_format, args.kmer_size)
