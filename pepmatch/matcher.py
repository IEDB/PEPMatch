import os
import _pickle as pickle
import sqlite3
import pandas as pd
import numpy as np

from typing import Union
from collections import Counter, defaultdict

from .helpers import parse_fasta, split_sequence, extract_metadata
from .preprocessor import Preprocessor
from .hamming import hamming


NUM_OUTPUT_COLUMNS = 14
VALID_OUTPUT_FORMATS = ['dataframe', 'csv', 'tsv', 'xlsx', 'json', 'html']


class Matcher:
  """
  Object class that class takes in a query (either FASTA file or Python list)
  of peptides, a proteome file, and a maximum # of mismatches. It can either
  find exact matches, matches with mismatches, or the best match for a peptide,
  depending on these parameters.

  A k value can also be passed if the proteome is already preprocessed, which
  most of the time it will be - so pass the k value to save time.

  If best match is selected and mismatching is to be done, one of two things
  will happen: 
  1. A k is specified already - just run mismatching and get the best match.
  2. A k is not specified - run best match code which preprocesses the proteome
  multiple times to get the best match quickly.

  Optional: output and output_format arguments to write results to file.
  Supported formats are "csv", "tsv", "xlsx", "json", and "html"."""

  def __init__(
    self,
    query,
    proteome_file,
    max_mismatches=-1,
    k=0,
    preprocessed_files_path='.',
    best_match=False,
    output_format='csv',
    output_name='',
    sequence_version=True
  ):
    if k < 0 or k == 1:
      raise ValueError('Invalid k value given. k cannot be negative or 1.')

    if output_format not in VALID_OUTPUT_FORMATS:
      raise ValueError(
        'Invalid output format, please choose `dataframe`, '
        '`csv`, `tsv`, `xlsx`, `json`, or `html`.'
      )

    # initialize query and output name based on input type
    self.query = self._initialize_query(query, proteome_file, output_name)
    self.proteome = proteome_file
    self.proteome_name = str(proteome_file).split('/')[-1].split('.')[0]

    # discontinuous epitopes and linear epitopes handling - store separately
    self.discontinuous_epitopes = self._find_discontinuous_epitopes()
    self.query = self._clean_query()
    
    assert self.query or self.discontinuous_epitopes, 'Query is empty.'
    
    self.lengths = sorted(set(len(peptide) for peptide in self.query))
    self.max_mismatches = max_mismatches
    self.preprocessed_files_path = preprocessed_files_path
    self.best_match = best_match
    self.output_format = output_format
    self.sequence_version = sequence_version

    # select format based on # of mismatches to pass to Preprocessor
    # SQLite for exact matching - pickle for mismatching
    self.preprocess_format = 'sql' if not max_mismatches else 'pickle'
    
    # initialize k and k_specified
    if self.query:
      self.k, self.k_specified = self._initialize_k(k)

    # for mismatching, if no k is specified, batch the peptides by length
    # then later we will use the ideal k to search them
    if max_mismatches > 0:
      self.batched_peptides = self._batch_query()

    # best match where no mismatches is specified means we will preprocess
    # the proteome by the smallest length, try to find exact matches,
    # then preprocess by half the length and search for more and repeat
    # until every peptide in the query has a match
    if max_mismatches == -1:
      self.ks = self._best_match_ks()


  def _initialize_query(
    self, query: Union[str, list], proteome_file: str, output_name: str
  ) -> list:
    """Initialize query and output name based on input type.
    
    Args:
      query: either a FASTA file or a Python list of peptides.
      proteome_file: the proteome FASTA file.
      output_name: the name of the output file."""
    if isinstance(query, list):
      if output_name:
        self.output_name = output_name
      else:
        self.output_name = 'PEPMatch_results'
      return [seq.upper() for seq in query]
    
    else: # parse from FASTA if not Python list
      parsed_query = [str(record.seq) for record in parse_fasta(query)]
      if output_name:
        self.output_name = output_name
      else: # output_name = query_name_to_proteome_name
        query_name = str(query).split('/')[-1].split('.')[0]
        proteome_name = str(proteome_file).split('/')[-1].split('.')[0]
        self.output_name = f'{query_name}_to_{proteome_name}'

      return [seq.upper() for seq in parsed_query]


  def _find_discontinuous_epitopes(self) -> list:
    """Find discontinuous epitopes in query and store separately."""
    discontinuous_epitopes = []
    for peptide in self.query:
      try:
        discontinuous_epitope = [(x[0], int(x[1:])) for x in peptide.split(', ')]
        discontinuous_epitopes.append(discontinuous_epitope)
      except ValueError:
        continue

    return discontinuous_epitopes


  def _clean_query(self) -> list:
    """Remove discontinous epitopes from query."""
    epitopes = self.discontinuous_epitopes
    discontinuous_epitope_strings = [
      ', '.join([x[0] + str(x[1]) for x in epitopes[i]]) for i in range(len(epitopes))
    ]
    query_filtered = set(self.query) - set(discontinuous_epitope_strings)
    return list(query_filtered)
  

  def _initialize_k(self, k: int) -> tuple:
    """Initialize k and k_specified values based on k and max_mismatches input.
    
    Args:
      k: the k-mer length to use for matching, 0 if not specified."""
    
    if k > 1:
      return k, True
    else: # use the length of the shortest peptide for exact match 
      if not self.max_mismatches and not k:
        return self.lengths[0], False
      else:
        return 0, False


  def _batch_query(self) -> dict:
    """Batch peptides together by ideal k so it can be used when searching.
    If k is specified, just return the query as a dictionary with k as the key."""

    if self.k_specified:
      return {self.k: self.query}
    else:
      batched_peptides = defaultdict(list)
      for seq in self.query:
        key = len(seq) // (self.max_mismatches + 1)
        batched_peptides[key].append(seq)
      
      return dict(batched_peptides)


  def _best_match_ks(self) -> list:
    """For the special case where mismatching is to be done, k is not 
    specified and best match is selected, we need to get all the k values
    to preprocess the proteome multiple times. Starting with the length 
    of the smallest peptide and then halving until we get to 2."""

    initial_k = self.lengths[0]
    ks = [initial_k]

    # halve the length and add k values until we get to 2 
    # and make sure 2 is a k in the list
    while initial_k > 2:
      initial_k //= 2

      if initial_k == 1:
        ks.append(2)
      else:
        ks.append(initial_k)

    # mismatching function uses batched peptides, for best match, batch into one
    self.batched_peptides = {0: self.query}

    return ks


  def match(self) -> list:
    """Overarching function that calls the appropriate search matching function
    based on the parameters."""

    if self.query:
      
      if self.max_mismatches == -1:
        query_df = self._dataframe_matches(self.best_match_search())

      elif self.max_mismatches > 0:
        query_df = self._dataframe_matches(self.mismatch_search())
      
      else:
        query_df = self._dataframe_matches(self.exact_match_search())
    
    else:
      query_df = pd.DataFrame()

    # search for discontinuous epitopes if they exist
    if self.discontinuous_epitopes:
      discontinuous_df = self._dataframe_matches(self.discontinuous_search())
    else:
      discontinuous_df = pd.DataFrame()

    # combine the dataframes if both query and discontinuous epitopes exist
    df = pd.concat([query_df, discontinuous_df], ignore_index=True)

    # return a dataframe instead of outputting file if specified
    if self.output_format == 'dataframe':
      return df
    
    else:
      self._output_matches(df)


  def exact_match_search(self) -> list:
    """Using preprocessed data within a SQLite database and the query peptides,
    find the peptide matches within the proteome without any residue 
    substitutions."""

    preprocessed_db = os.path.join(
      self.preprocessed_files_path, self.proteome_name + '.db'
    )
    if not os.path.isfile(preprocessed_db):
      Preprocessor(self.proteome).sql_proteome(self.k)

    kmers_table_name = f'{self.proteome_name}_{str(self.k)}mers'
    metadata_table_name = f'{self.proteome_name}_metadata'

    conn = sqlite3.connect(
      os.path.join(self.preprocessed_files_path, f'{self.proteome_name}.db')
    )
    cursor = conn.cursor()

    all_matches = []
    for peptide in self.query:
      
      if len(peptide) < self.k:
        continue

      # split peptide into kmers - only use kmers necessary that overlap entire peptide
      all_kmers = split_sequence(peptide, self.k)
      target_kmers = self._get_target_kmers(all_kmers)

      # SQL fetch
      sql_placeholders = ', '.join('?' * len(target_kmers))
      sql_query = f"""
                   SELECT kmer, idx FROM "{kmers_table_name}" 
                   WHERE kmer IN ({sql_placeholders})
                   """
      cursor.execute(sql_query, target_kmers)
      kmer_indexes = cursor.fetchall()

      kmer_hit_list = []
      for kmer, index in kmer_indexes:
        kmer_positions = [i for i, k in enumerate(all_kmers) if k == kmer] # this accounts for duplicate k-mers
        for pos in kmer_positions:
          kmer_hit_list.append(index - pos)

      matches = []
      sum_hits = Counter(kmer_hit_list)
      for hit, count in sum_hits.items():
        if count == len(target_kmers): # number of index recordings that 
          matches.append(hit)  # agree with the number of kmers used

      processed_matches = self._process_exact_matches(
        peptide, matches, cursor, metadata_table_name
      )
      all_matches.extend(processed_matches)

    cursor.close()
    conn.close()

    return all_matches


  def _get_target_kmers(self, all_kmers: list) -> list:
    """Return the target kmers that overlap the entire peptide.
    
    Args:
      all_kmers: all possible kmers of the query peptide for a given k."""
    
    if len(all_kmers) == self.k:
        return all_kmers

    target_kmers = all_kmers[::self.k]
    if all_kmers[-1] != target_kmers[-1]:
        target_kmers.append(all_kmers[-1])
    return target_kmers


  def _process_exact_matches(
    self, peptide: str, matches: list, cursor: sqlite3.Cursor, metadata_table_name: str
  ) -> list:
    """Extract all metadata for the exact matches and return as a list of tuples.
    
    Args:
      peptide: the query peptide.
      matches: the list of exact matches for the peptide.
      cursor: the cursor object to execute SQL queries.
      metadata_table_name: the name of the metadata table in the database."""
    
    all_matches = []
    if not matches:
      all_matches.append((peptide,) + (np.nan,) * NUM_OUTPUT_COLUMNS)
    else:
      for match in matches:
        protein_number = (match - (match % 1000000)) // 1000000
        query = f"""SELECT * 
                 FROM "{metadata_table_name}" 
                 WHERE protein_number = "{protein_number}"
                 """
        cursor.execute(query)
        protein_data = cursor.fetchone()
        match_data = (
          peptide,                          # query peptide
          peptide,                          # matched peptide (same as query)
          protein_data[1],                  # protein ID
          protein_data[2],                  # protein name
          protein_data[3],                  # species
          protein_data[4],                  # taxon ID
          protein_data[5],                  # gene symbol
          0,                                # 0 mismatches for exact match
          [],                               # mutated positions (none)
          (match % 1000000) + 1,            # index start
          (match % 1000000) + len(peptide), # index end
          protein_data[6],                  # protein existence level
          protein_data[7],                  # sequence version
          protein_data[8],                  # gene priority flag
          protein_data[9])                  # swissprot flag
        
        all_matches.append(match_data)

    return all_matches


  def mismatch_search(self) -> list:
    """Using preprocessed data within a serialized pickle files, the query
    peptides, and a maximum number of residue substitutions, find all peptide
    matches up to and including the maximum number of residue substitutions.

    This will utilize the pigeon hole principle to find all possible matches.
    We first search for any k-mer exact matches in the proteome and then
    check the left and right k-mer neighbors for mismatches. If the number of
    mismatches is less than or equal to the max, then it's a match."""

    all_matches = []
    for k, peptides in self.batched_peptides.items(): # iterate through each batch
      if not self.k_specified:                        # k: [peptides]
        self.k = k

      try:
        kmer_dict, rev_kmer_dict, metadata_dict = self._read_pickle_files()
      except FileNotFoundError: # do preprocessing if pickle files don't exist
        Preprocessor(self.proteome).pickle_proteome(self.k)
        kmer_dict, rev_kmer_dict, metadata_dict = self._read_pickle_files()

      for peptide in peptides:

        all_kmers = split_sequence(peptide, self.k)

        if len(peptide) % self.k == 0: # faster search if possible
          matches = self._find_even_split_matches(
            all_kmers, kmer_dict, rev_kmer_dict, len(peptide)
          )

        else: # slower search if necessary
          matches = self._find_uneven_split_matches(
            all_kmers, kmer_dict, rev_kmer_dict, len(peptide)
          )

        processed_matches = self._process_mismatch_matches(
          peptide, matches, metadata_dict
        )
        all_matches.extend(processed_matches)

    return all_matches


  def _read_pickle_files(self) -> tuple:
    """Read in the already created pickle files for each dictionary in the
    preprocessing step."""

    with open(os.path.join(self.preprocessed_files_path, 
      f'{self.proteome_name}_{str(self.k)}mers.pkl'), 'rb') as f:
      kmer_dict = pickle.load(f)

    with open(os.path.join(self.preprocessed_files_path, 
      f'{self.proteome_name}_metadata.pkl'), 'rb') as f:
      metadata_dict = pickle.load(f)

    rev_kmer_dict = {i: k for k, v in kmer_dict.items() for i in v}

    return kmer_dict, rev_kmer_dict, metadata_dict


  def _find_even_split_matches(
    self, kmers: list, kmer_dict: dict, rev_kmer_dict: dict, peptide_length: int
  ) -> list:
    """If the peptide length is evenly divisible by k, perform faster search. 
    Check the associated k-mers for the left and right neighbors of any 
    exact matching k-mers in the proteome.
    
    EXAMPLE: "YLLDLHSYL", k=3 -> check only ["YLL", "DLH", "SYL"]. 
    If, say "DLH" is found somewhere in the proteome, check the neighboring
    k-mers to the left and right to see if "YLL" (left neighbor) and "SYL" 
    (right neighbor) have mismatches with the proteome k-mers, return any
    matches that are less than or equal to self.max_mismatches.

    Args:
      kmers: the k-mers of the query peptide.
      kmer_dict: the k-mer -> index dictionary.
      rev_kmer_dict: the index -> k-mer dictionary.
      peptide_length: the length of the query peptide."""
    
    matches = set()
    for idx in range(0, len(kmers), self.k): # gets only the k-mers to check
      kmer_hits = kmer_dict.get(kmers[idx], [])
      for kmer_hit in kmer_hits:

        mismatches = 0
        mismatches = self._check_left_neighbors(
          kmers, idx, kmer_hit, rev_kmer_dict, mismatches
        )
        
        mismatches = self._check_right_neighbors(
          kmers, idx, kmer_hit, rev_kmer_dict, mismatches
        )
        
        if mismatches <= self.max_mismatches:
          matched_peptide = '' # add k-mers to get the matched peptide
          for i in range(0, peptide_length, self.k): 
            matched_peptide += rev_kmer_dict[kmer_hit - idx + i]

          matches.add((matched_peptide, mismatches, kmer_hit - idx))

          if self.best_match and not mismatches: # can't have a better match
            return list(matches)
        
        else:
            continue
    
    return list(matches)


  def _find_uneven_split_matches(
    self, kmers: list, kmer_dict: dict, rev_kmer_dict: dict, peptide_length: int
  ) -> list:
    """If the peptide length is NOT evenly divisible by k, perform slow search. 
    Check the associated residues for the left and right neighbors of any 
    exact matching k-mers in the proteome.
    
    EXAMPLE: "YLLDLHSYL", k=3 -> check all k-mers: 
    ["YLL", "LLD", "LDL", "DLH", "LHS", "HSY", "SYL"].
    
    If, say "DLH" is found somewhere in the proteome, check left "L" from "LDL",
    left "L" from "LLD", and "Y" from "YLL" to for any mismatches. Then do the 
    same for the residues of the right neighbors. Return any matches that are
    less than or equal to self.max_mismatches.

    Args:
      kmers: the k-mers of the query peptide.
      kmer_dict: the k-mer -> index dictionary.
      rev_kmer_dict: the index -> k-mer dictionary.
      peptide_length: the length of the query peptide."""

    matches = set()
    for idx in range(0, len(kmers)): # every k-mer
      kmer_hits = kmer_dict.get(kmers[idx], [])
      for kmer_hit in kmer_hits:

        mismatches = 0
        mismatches = self._check_left_residues(
          kmers, idx, kmer_hit, rev_kmer_dict, mismatches
        )

        mismatches = self._check_right_residues(
          kmers, idx, kmer_hit, rev_kmer_dict, mismatches
        )
        
        if mismatches <= self.max_mismatches:
          matched_peptide = ''
          for i in range(0, peptide_length, self.k):
            if i + self.k >= peptide_length: # handle last k-mer when uneven split
              remaining_res = peptide_length - i
              last_kmer = rev_kmer_dict[kmer_hit - idx + i - (self.k - remaining_res)]
              matched_peptide += last_kmer[-remaining_res:] # attach remaining residues
            else:
              matched_peptide += rev_kmer_dict[kmer_hit - idx + i]

          matches.add((matched_peptide, mismatches, kmer_hit - idx))

          if self.best_match and not mismatches: # can't have a better match
            return matches

    return matches


  def _check_left_neighbors(
    self, kmers: list, idx: int, kmer_hit: int, rev_kmer_dict: dict, mismatches: int
  ) -> int:
    """Get mismatches of left k-mer neighbors in proteome.
    
    Args:
      kmers: the k-mers of the query peptide.
      idx: the index of the k-mer in the query peptide.
      kmer_hit: the index of the k-mer hit in the proteome.
      rev_kmer_dict: the index -> k-mer dictionary.
      mismatches: the number of mismatches so far."""
    
    for i in range(0, idx, self.k):
      if mismatches > self.max_mismatches:
        return 100

      kmer_to_check = rev_kmer_dict.get(kmer_hit + i - idx)
      if kmer_to_check is not None:
        mismatches += hamming(kmer_to_check, kmers[i], self.max_mismatches)      
      else:
        return 100
    
    return mismatches


  def _check_right_neighbors(
    self, kmers: list, idx: int, kmer_hit: int, rev_kmer_dict: dict, mismatches: int
  ) -> int:
    """Get mismatches of right k-mer neighbors in proteome.
    
    Args:
      kmers: the k-mers of the query peptide.
      idx: the index of the k-mer in the query peptide.
      kmer_hit: the index of the k-mer hit in the proteome.
      rev_kmer_dict: the index -> k-mer dictionary.
      mismatches: the number of mismatches so far."""
    
    for i in range(self.k + idx, len(kmers), self.k):
      if mismatches > self.max_mismatches:
        return 100
      
      kmer_to_check = rev_kmer_dict.get(kmer_hit + i - idx)
      if kmer_to_check is not None:
        mismatches += hamming(kmer_to_check, kmers[i], self.max_mismatches)
      else:
        return 100
    
    return mismatches


  def _check_left_residues(
    self, kmers: list, idx: int, kmer_hit: int, rev_kmer_dict: dict, mismatches: int
  ) -> int:
    """Get mismatches of left residues of left k-mer neighbors in proteome.
    
    Args:
      kmers: the k-mers of the query peptide.
      idx: the index of the k-mer in the query peptide.
      kmer_hit: the index of the k-mer hit in the proteome.
      rev_kmer_dict: the index -> k-mer dictionary.
      mismatches: the number of mismatches so far."""
    
    for i in range(0, idx):
      if mismatches > self.max_mismatches:
        return 100

      kmer_to_check = rev_kmer_dict.get(kmer_hit + i - idx)
      if kmer_to_check is not None:
        if kmer_to_check[0] != kmers[i][0]:
          mismatches += 1
      else:
        return 100
    
    return mismatches


  def _check_right_residues(
    self, kmers: list, idx: int, kmer_hit: int, rev_kmer_dict: dict, mismatches: int
  ) -> int:
    """Get mismatches of right residues of right k-mer neighbors in proteome.
    
    Args:
      kmers: the k-mers of the query peptide.
      idx: the index of the k-mer in the query peptide.
      kmer_hit: the index of the k-mer hit in the proteome.
      rev_kmer_dict: the index -> k-mer dictionary.
      mismatches: the number of mismatches so far."""
    
    for i in range(idx + 1, len(kmers)):
      if mismatches > self.max_mismatches: # check before since we checked left first
        return 100

      kmer_to_check = rev_kmer_dict.get(kmer_hit + i - idx)
      if kmer_to_check is not None:
        if kmer_to_check[-1] != kmers[i][-1]:
          mismatches += 1
      else:
        return 100
    
    return mismatches


  def _process_mismatch_matches(
    self, peptide: str, matches: list, metadata_dict: dict
  ) -> list:
    """Extract all metadata for the mismatch matches and return as a list of tuples.
    
    Args:
      peptide: the query peptide.
      matches: the list of mismatch matches for the peptide.
      metadata_dict: the protein number -> metadata dictionary."""
    
    all_matches = []
    if not matches:
      all_matches.append((peptide,) + (np.nan,) * NUM_OUTPUT_COLUMNS)
    else:
      for match in matches:
        metadata_key = (match[2] - (match[2] % 1000000)) // 1000000
        metadata = metadata_dict[metadata_key]

        mutated_positions = [
          i+1 for i in range(len(peptide)) if peptide[i] != match[0][i]
        ]
        index_start = int((match[2] % 1000000) + 1)
        index_end = int((match[2] % 1000000) + len(peptide))

        taxon_id = int(metadata[3]) if metadata[3] else np.nan
        pe_level = int(metadata[5]) if metadata[5] else np.nan

        match_data = (
          peptide,                      # query peptide
          match[0],                     # matched peptide
          metadata[0],                  # protein ID
          metadata[1],                  # protein name
          metadata[2],                  # species
          taxon_id,                     # taxon ID
          metadata[4],                  # gene symbol
          int(match[1]),                # mismatches count
          mutated_positions,            # mutated positions
          index_start,                  # index start
          index_end,                    # index end
          pe_level,                     # protein existence level
          metadata[6],                  # sequence version
          metadata[7],                  # gene priority flag
          metadata[8]                   # swissprot flag
        )
        all_matches.append(match_data)

    return all_matches


  def best_match_search(self) -> list:
    """After calculating the splits we would need (starting with lowest peptide
    length in query), we can then call the mismatch_search function. We use
    the maximum # of mismatches calculated from the each length and split."""

    all_matches = []
    for k in self.ks:
      self.k = k
      self.k_specified = True

      # optimal # of mismatches for k and smallest peptide
      self.max_mismatches = (self.lengths[0] // self.k) - 1 

      if self.max_mismatches == 0:
        matches = self.exact_match_search()
      else:
        matches = self.mismatch_search()

      # separate out peptides that did not match into a new query
      self.new_query = []
      for match in matches:
        if not pd.isna(match[1]):
          all_matches.append(match)
        else:
          self.new_query.append(match[0])

      if not self.new_query: # all peptides have a match
        return all_matches
      else: # reset query with peptides that did not match
        self.query = self.new_query

      # last k - keep searching until all peptides have a match
      if self.k == 2:
        self.max_mismatches += 1
        while self.new_query:
          self.query = self.new_query
          matches = self.mismatch_search()
          self.new_query = []
          for match in matches:
            if not pd.isna(match[1]):
              all_matches.append(match)
            else:
              self.new_query.append(match[0])
          self.max_mismatches += 1
      
      # reset batched peptides
      self.batched_peptides = {0: self.query}

    return all_matches


  def discontinuous_search(self) -> list:
    """Find matches for discontinuous epitopes. Loops through every protein
    in the proteome and checks if the residues at the given positions match
    the query epitope residues up to the maximum number of mismatches."""

    all_matches = []
    for dis_epitope in self.discontinuous_epitopes:
      match = False
      for protein_record in parse_fasta(self.proteome):
        try:
          residue_matches = sum(
            [x[0] == protein_record.seq[x[1] - 1] for x in dis_epitope]
          )
          if residue_matches >= (len(dis_epitope) - self.max_mismatches):
            match = True
            metadata = extract_metadata(protein_record, False)
            match_data = (
              ', '.join(                                    # query epitope
                [x[0] + str(x[1]) for x in dis_epitope]),
              ', '.join(                                    # matched epitope
                [protein_record.seq[x[1] - 1] + str(x[1]) for x in dis_epitope]), 
              metadata[0],                                  # protein ID
              metadata[1],                                  # protein name
              metadata[2],                                  # species
              metadata[3],                                  # taxon ID
              metadata[4],                                  # gene symbol
              len(dis_epitope) - residue_matches,           # mismatches count
              # mutated positions
              [x[1] for x in dis_epitope if x[0] != protein_record.seq[x[1] - 1]],
              dis_epitope[0][1],                            # index start
              dis_epitope[-1][1],                           # index end
              metadata[5],                                  # protein existence level
              metadata[6],                                  # sequence version
              metadata[7],                                  # gene priority flag
              metadata[8])                                  # swissprot flag
            
            all_matches.append(match_data)
        
        except IndexError:
          continue

      if not match:
        all_matches.append(
          (', '.join([x[0] + str(x[1]) for x in dis_epitope]),) + (np.nan,) * NUM_OUTPUT_COLUMNS
        )

    return all_matches

  def _take_matches(self, df: pd.DataFrame, search: str, min_max: str = 'min'):
    return (df.groupby(['Query Sequence'])[search].transform(min_max) == df[search]) | df[search].isna()


  def _dataframe_matches(self, all_matches: list) -> pd.DataFrame:
    """Return Pandas dataframe of the results.
    
    Args:
      all_matches: the list of all matches for all peptides."""
    df = pd.DataFrame(
      all_matches,
      columns=[
        'Query Sequence', 'Matched Sequence', 'Protein ID', 'Protein Name', 'Species',
        'Taxon ID', 'Gene', 'Mismatches', 'Mutated Positions','Index start', 'Index end',
        'Protein Existence Level', 'Sequence Version', 'Gene Priority', 'SwissProt Reviewed'
      ]
    )

    if self.best_match:
      def filter_fragments(group: pd.DataFrame) -> pd.DataFrame:
        """Takes out matches with 'Fragment' in the protein name if there are other
        matches without 'Fragment'.

        Args:
          group: the group of matches for a given query peptide.
        """
        group['Protein Name'] = group['Protein Name'].astype(str)
        no_fragments = group[~group['Protein Name'].str.contains('Fragment')]
        if len(no_fragments) > 0:
          return no_fragments
        else:
          return group

      # take matches with the least number of mismatches
      df = df[self._take_matches(df, "Mismatches")]

      # fill NaN values with 0 - otherwise pandas will drop the rows
      df[['Gene Priority', 'Protein Existence Level']] = df[
         ['Gene Priority', 'Protein Existence Level']
      ].fillna(value=0)

      # take matches that are in the gene priority proteome
      df = df[self._take_matches(df, "Gene Priority", "max")]

      # and with the best protein existence leve
      df = df[self._take_matches(df, "Protein Existence Level")]

      # filter fragments
      df = df.groupby(
        'Query Sequence', group_keys=False
      ).apply(filter_fragments).reset_index(drop=True)

      # sort values by protein ID and drop duplicates, guaranteeing same results 
      df.sort_values(by=['Query Sequence', 'Protein ID', 'Index start'], inplace=True)
      df.drop_duplicates(['Query Sequence'], inplace=True)

    if self.sequence_version:
      # combine protein ID and sequence version
      df['Sequence Version'] = df['Sequence Version'].apply(
        lambda x: f'.{int(x)}' if not pd.isna(x) else ''
      )   
      df['Protein ID'] = df['Protein ID'] + df['Sequence Version']
    
    # drop "Sequence Version" column
    df.drop(columns=['Sequence Version'], inplace=True)

    # force integers on some columns
    int_cols = [
      'Mismatches', 'Index start', 'Index end', 'Protein Existence Level'
    ]
    df[int_cols] = df[int_cols].astype('Int64')
    # TODO: convert taxon ID to int if possible
    return df

  def _get_output_function(self, df: pd.DataFrame):
    """Returns the df conversion function, given the required output format
    
    Args:
      df: the dataframe of the matches."""
    return [df.to_csv, df.to_csv, df.to_excel, df.to_json, df.to_html][
      ['csv', 'tsv', 'xlsx', 'json', 'html'].index(self.output_format)
    ]

  def _output_matches(self, df: pd.DataFrame) -> None:
    """Write Pandas dataframe to format that is specified

    Args:
      df: the dataframe of the matches."""
    
    # appends '.' + filetype if the name does not already contain it
    path = self.output_name.__str__()
    if path.split('.')[-1] != self.output_format:
      path += f".{self.output_format}"

    output = self._get_output_function(df)
    if self.output_format == 'tsv':
      output(path, sep='\t', index=False)
    elif self.output_format == 'json':
      output(path, orient='split', index=False)
    else:
      output(path, index=False)
