import _pickle as pickle
import sqlite3
import os

from .helpers import parse_fasta, split_sequence, extract_metadata


class Preprocessor:
  """Class that takes in a proteome FASTA file and preprocesses the file to be
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
                           this proteome. Otherwise, GP=0."""
  def __init__(
    self,
    proteome,
    proteome_name='',
    preprocessed_files_path='.',
    gene_priority_proteome=''
  ):

    if not os.path.isdir(preprocessed_files_path):
      raise ValueError('Directory specified does not exist: ', preprocessed_files_path)

    self.preprocessed_files_path = preprocessed_files_path

    # appending GP=# to the header of the proteome
    self.proteome = self._append_gp_to_header(proteome, gene_priority_proteome)

    # extracting the proteome name from the file path if not specified
    if not proteome_name:
      self.proteome_name = str(proteome).split('/')[-1].split('.')[0]
    else:
      self.proteome_name = proteome_name

    # extract all the data from the proteome
    self.all_seqs, self.all_metadata = self._get_data_from_proteome()


  def _append_gp_to_header(
    self, proteome: str, gene_priority_proteome: str
  ) -> list:
    """Appends GP=1 or GP=0 to the FASTA header of a proteome depending on
    if the protein is in the gene priority proteome.

    Args:
      proteome: path to proteome FASTA file.
      gene_priority_proteome: path to gene priority proteome FASTA file."""
    
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


  def _get_data_from_proteome(self) -> tuple:
    """Extract all the data from a FASTA file and returns two lists:

    Use extract_metadata_from_fasta_header in helpers.py.

    1. A list of sequences
    2. A list of metadata in tuples. The metadata includes the protein ID,
       protein name, species, taxon ID, gene, protein existence level,
       sequence version, and gene priority label."""
    
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


  def sql_proteome(self, k: int) -> None:
    """Writes the kmers_table and metadata_table to a SQLite database.
    
    Args:
      k: k-mer length to split the proteome into."""
  
    # create table names
    kmers_table = f'{self.proteome_name}_{str(k)}mers'
    metadata_table = f'{self.proteome_name}_metadata'
    
    # connect to database
    db_path = os.path.join(self.preprocessed_files_path, f'{self.proteome_name}.db')
    with sqlite3.connect(db_path) as conn:
      cursor = conn.cursor()

      # check if kmers table already exists and exit if it does
      # for some reason writing to the same table multiple times messes up results
      if cursor.execute(
        f'SELECT name FROM sqlite_master WHERE type="table" AND name="{kmers_table}";'
      ).fetchone():
        return

      self._create_tables(cursor, kmers_table, metadata_table)
      self._insert_kmers(cursor, kmers_table, k)
      self._insert_metadata(cursor, metadata_table)
      self._create_indexes(cursor, kmers_table, metadata_table)

      conn.commit()


  def _create_tables(
    self, cursor: sqlite3.Cursor, kmers_table: str, metadata_table: str
  ) -> None:
    """Creates the kmers_table and metadata_table in the SQLite database.

    Args:
      cursor: cursor object to execute SQL commands.
      kmers_table: name of the k-mers table.
      metadata_table: name of the metadata table."""

    cursor.execute(
      f'CREATE TABLE IF NOT EXISTS "{kmers_table}" ('\
        'kmer TEXT NOT NULL,'\
        'idx  INTEGER NOT NULL)'
      )
    cursor.execute(
      f'CREATE TABLE IF NOT EXISTS "{metadata_table}" ('\
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


  def _insert_kmers(self, cursor: sqlite3.Cursor, kmers_table: str, k: int) -> None:
    """Inserts the k-mers into the kmers_table.

    Args:
      cursor: cursor object to execute SQL commands.
      kmers_table: name of the k-mers table.
      k: k-mer length to split the proteome into."""
    
    kmer_rows = []
    for protein_count, seq in enumerate(self.all_seqs):
        for j, kmer in enumerate(split_sequence(seq, k)):
            kmer_rows.append((kmer, (protein_count + 1) * 1000000 + j))
    cursor.executemany(f'INSERT INTO "{kmers_table}" VALUES (?, ?)', kmer_rows)


  def _insert_metadata(self, cursor: sqlite3.Cursor, metadata_table: str) -> None:
    """Inserts the metadata into the metadata_table.

    Args:
      cursor: cursor object to execute SQL commands.
      metadata_table: name of the metadata table."""
    
    cursor.executemany(
      f'INSERT INTO "{metadata_table}" VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)', 
      self.all_metadata
    )


  def _create_indexes(
    self, cursor: sqlite3.Cursor, kmers_table: str, metadata_table: str
  ) -> None:
    """Creates indexes for the kmers_table and metadata_table.

    Args:
      cursor: cursor object to execute SQL commands.
      kmers_table: name of the k-mers table.
      metadata_table: name of the metadata table."""
    
    cursor.execute(
      f'CREATE INDEX IF NOT EXISTS "{kmers_table}_kmer_idx" ON "{kmers_table}" (kmer)'
    )
    cursor.execute(
      f'CREATE INDEX IF NOT EXISTS "{metadata_table}_protein_number_idx" '
      f'ON "{metadata_table}" (protein_number)'
    )


  def pickle_proteome(self, k: int) -> None:
    """Pickles the proteome into a dictionary of k-mers and a dictionary of metadata.

    Args:
      k: k-mer length to split the proteome into."""
    
    # create kmer_dict and metadata_dict out of self.all_seqs and self.all_metadata
    kmer_dict = {}
    for protein_count, seq in enumerate(self.all_seqs):
      for j, kmer in enumerate(split_sequence(seq, k)):
        if kmer in kmer_dict.keys(): # add index to k-mer list
          kmer_dict[kmer].append((protein_count + 1) * 1000000 + j) 
        else: # create entry for new k-mer
          kmer_dict[kmer] = [(protein_count + 1) * 1000000 + j] 
    
    metadata_dict = {}
    for data in self.all_metadata:
      metadata_dict[data[0]] = data[1:]
    
    # write kmer_dict and metadata_dict to pickle files
    with open(os.path.join(self.preprocessed_files_path, 
      f'{self.proteome_name}_{str(k)}mers.pkl'), 'wb') as f:
      pickle.dump(kmer_dict, f)

    with open(os.path.join(self.preprocessed_files_path, 
      f'{self.proteome_name}_metadata.pkl'), 'wb') as f:
      pickle.dump(metadata_dict, f)


  def preprocess(self, preprocess_format: str, k: int) -> None:
    """Preprocesses the proteome and stores it in the specified format.

    Args:
      preprocess_format: format to store the proteome in (sql or pickle).
      k: k-mer length to split the proteome into."""
    
    if preprocess_format not in ('sql', 'pickle'):
      raise AssertionError(
        'Unexpected value of preprocessing format:', preprocess_format
      )

    assert k >= 2, 'k-sized split is invalid. Cannot be less than 2.'

    # store data based on format specified
    if preprocess_format == 'pickle':
      self.pickle_proteome(k)
    elif preprocess_format == 'sql':
      self.sql_proteome(k)