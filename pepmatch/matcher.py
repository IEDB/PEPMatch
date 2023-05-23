import os
import argparse
import _pickle as pickle
import sqlite3
import pandas as pd
import numpy as np

from collections import Counter, defaultdict
from Levenshtein import hamming

from .helpers import parse_fasta, split_sequence, extract_metadata
from .preprocessor import Preprocessor


VALID_OUTPUT_FORMATS = ['dataframe', 'csv', 'xlsx', 'json', 'html']

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
  Supported formats are "csv", "xlsx", "json", and "html".
  """
  def __init__(self,
               query,
               proteome_file,
               max_mismatches=-1,
               k=0,
               preprocessed_files_path='.',
               best_match=False,
               output_format='csv',
               output_name=''):

    if k < 0 or k == 1:
      raise ValueError('Invalid k value given. k cannot be negative or 1.')

    if output_format not in VALID_OUTPUT_FORMATS:
      raise ValueError('Invalid output format, please choose `dataframe`, `csv`, `xlsx`, `json`, or `html`.')

    # initialize query and output name based on input type
    self.query = self._initialize_query(query, proteome_file, output_name)
    self.proteome = proteome_file
    self.proteome_name = proteome_file.split('/')[-1].split('.')[0]

    # discontinuous epitopes and linear epitopes handling - store separately
    self.discontinuous_epitopes = self._find_discontinuous_epitopes()
    self.query = self._clean_query()
    assert self.query, 'Query is empty. Please check your input.'
    
    self.lengths = sorted(set(len(peptide) for peptide in self.query))
    self.max_mismatches = max_mismatches
    self.preprocessed_files_path = preprocessed_files_path
    self.best_match = best_match
    self.output_format = output_format

    # select format based on # of mismatches to pass to Preprocessor
    # SQLite for exact matching - pickle for mismatching
    self.preprocess_format = 'sql' if not max_mismatches else 'pickle'
    
    # initialize k and k_specified
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

  def _initialize_query(self, query, proteome_file, output_name):
    """Initialize query and output name based on input type."""
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
        self.output_name = f"{query.split('/')[-1].split('.')[0]}_to_{proteome_file.split('/')[-1].split('.')[0]}"
      
      return [seq.upper() for seq in parsed_query]

  def _find_discontinuous_epitopes(self):
    """Find discontinuous epitopes in query and store separately."""
    discontinuous_epitopes = []
    for peptide in self.query:
      try:
        discontinuous_epitope = [(x[0], int(x[1:])) for x in peptide.split(', ')]
        discontinuous_epitopes.append(discontinuous_epitope)
      except ValueError:
        continue

    return discontinuous_epitopes

  def _clean_query(self):
    """Remove discontinous epitopes from query."""
    query_without_discontinuous_epitopes = set(self.query) - set([', '.join([x[0] + str(x[1]) for x in self.discontinuous_epitopes[i]]) for i in range(len(self.discontinuous_epitopes))])
    return list(query_without_discontinuous_epitopes)
  
  def _initialize_k(self, k):
    """Initialize k and k_specified values based on k and max_mismatches input."""
    if k > 1:
      return k, True
    else: # use the length of the shortest peptide when k is not specified and exact matching is requested
      if not self.max_mismatches and not k:
        return self.lengths[0], False
      else:
        return 0, False

  def _batch_query(self):
    """
    Batch peptides together by ideal k so it can be used when searching.
    If k is specified, just return the query as a dictionary with k as the key.
    """
    if self.k_specified:
      return {self.k: self.query}
    else:
      batched_peptides = defaultdict(list)
      for seq in self.query:
        key = len(seq) // (self.max_mismatches + 1)
        batched_peptides[key].append(seq)
      
      return dict(batched_peptides)

  def _best_match_ks(self):
    """
    For the special case where mismatching is to be done, k is not 
    specified and best match is selected, we need to get all the k values
    to preprocess the proteome multiple times. Starting with the length 
    of the smallest peptide and then halving until we get to 2.
    """
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

  def match(self):
    """
    Overarching function that calls the appropriate search matching function
    based on the parameters.
    """
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

  def exact_match_search(self):
    """
    Using preprocessed data within a SQLite database and the query peptides,
    find the peptide matches within the proteome without any residue 
    substitutions. 
    """
    if not os.path.isfile(os.path.join(self.preprocessed_files_path, self.proteome_name + '.db')):
      Preprocessor(self.proteome).sql_proteome(self.k)

    kmers_table_name = f'{self.proteome_name}_{str(self.k)}mers'
    metadata_table_name = f'{self.proteome_name}_metadata'

    conn = sqlite3.connect(os.path.join(self.preprocessed_files_path, f'{self.proteome_name}.db'))
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
      sql_query = f'SELECT kmer, idx FROM "{kmers_table_name}" WHERE kmer IN ({sql_placeholders})'
      cursor.execute(sql_query, target_kmers)
      kmer_indexes = cursor.fetchall()

      kmer_hit_list = []
      for kmer, index in kmer_indexes:
        kmer_hit_list.append(index - all_kmers.index(kmer))

      peptide_matches = []
      sum_hits = Counter(kmer_hit_list)
      for hit, count in sum_hits.items():
        if count == len(target_kmers): # number of index recordings that agree with the number of kmers used
          peptide_matches.append(hit)

      processed_matches = self._process_exact_matches(peptide, peptide_matches, cursor, metadata_table_name)
      all_matches.extend(processed_matches)

    cursor.close()
    conn.close()

    return all_matches

  def _get_target_kmers(self, all_kmers):
    """Return the target kmers that overlap the entire peptide."""
    if len(all_kmers) == self.k:
        return all_kmers

    target_kmers = all_kmers[::self.k]
    if all_kmers[-1] != target_kmers[-1]:
        target_kmers.append(all_kmers[-1])
    return target_kmers

  def _process_exact_matches(self, peptide, peptide_matches, cursor, metadata_table_name):
    """Extract all metadata for the exact matches and return as a list of tuples."""
    all_matches = []
    if not peptide_matches:
      all_matches.append((peptide,) + (np.nan,) * 13)
    else:
      for match in peptide_matches:
        protein_number = (match - (match % 1000000)) // 1000000
        query = f'SELECT * FROM "{metadata_table_name}" WHERE protein_number = "{protein_number}"'
        cursor.execute(query)
        protein_data = cursor.fetchone()
        match_data = (
          peptide,                          # query peptide
          peptide,                          # matched peptide (same as query for exact)
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
          protein_data[8])                  # gene priority flag
        
        all_matches.append(match_data)

    return all_matches

  def mismatch_search(self):
    """
    Using preprocessed data within a serialized pickle files, the query
    peptides, and a maximum number of residue substitutions, find all peptide
    matches up to and including the maximum number of residue substitutions.

    This will utilize the pigeon hole principle to find all possible matches.
    We first search for any k-mer exact matches in the proteome and then
    check the left and right k-mer neighbors for mismatches. If the number of
    mismatches is less than or equal to the max, then it's a match.
    """
    all_matches = []
    for k, peptides in self.batched_peptides.items(): # iterate through each batch - k: [peptides]
      if not self.k_specified:
        self.k = k

      try: # read in the preprocessed pickle files - create reverse kmer dictionary for neighbor seaches - index: kmer
        kmer_dict, metadata_dict = self._read_pickle_files()
        rev_kmer_dict = {i: k for k, v in kmer_dict.items() for i in v}
      except FileNotFoundError: # do preprocessing if pickle files don't exist
        Preprocessor(self.proteome).pickle_proteome(self.k)
        kmer_dict, metadata_dict = self._read_pickle_files()
        rev_kmer_dict = {i: k for k, v in kmer_dict.items() for i in v}

      for peptide in peptides:

        all_kmers = split_sequence(peptide, self.k)

        # faster search if possible
        if len(peptide) % self.k == 0:
          peptide_matches = self._find_even_split_matches(all_kmers, kmer_dict, rev_kmer_dict, len(peptide))

        # slower search if necessary
        else:
          peptide_matches = self._find_uneven_split_matches(all_kmers, kmer_dict, rev_kmer_dict, len(peptide))

        processed_matches = self._process_mismatch_matches(peptide, peptide_matches, metadata_dict)
        all_matches.extend(processed_matches)
    
    return all_matches

  def _read_pickle_files(self):
    """
    Read in the already created pickle files for each dictionary in the
    preprocessing step.
    """
    with open(os.path.join(self.preprocessed_files_path, 
      f'{self.proteome_name}_{str(self.k)}mers.pickle'), 'rb') as f:
      kmer_dict = pickle.load(f)

    with open(os.path.join(self.preprocessed_files_path, 
      f'{self.proteome_name}_metadata.pickle'), 'rb') as f:
      names_dict = pickle.load(f)

    return kmer_dict, names_dict

  def _find_even_split_matches(self, kmers, kmer_dict, rev_kmer_dict, peptide_length):
    """
    If the peptide length is evenly divisible by k, perform faster search. 
    Check the associated k-mers for the left and right neighbors of any 
    exact matching k-mers in the proteome.
    
    EXAMPLE: "YLLDLHSYL", k=3 -> check only ["YLL", "DLH", "SYL"]. 
    If, say "DLH" is found somewhere in the proteome, check the neighboring
    k-mers to the left and right to see if "YLL" (left neighbor) and "SYL" 
    (right neighbor) have mismatches with the proteome k-mers, return any
    matches that are less than or equal to self.max_mismatches.
    """
    peptide_matches = set()
    for idx in range(0, len(kmers), self.k): # gets only the k-mers to check
      try: # try to find an exact match hit for each k-mer
        for kmer_hit in kmer_dict[kmers[idx]]:

          mismatches = 0
          mismatches = self._check_left_neighbors(kmers, idx, kmer_hit, rev_kmer_dict, mismatches)

          if not mismatches > self.max_mismatches:
            mismatches = self._check_right_neighbors(kmers, idx, kmer_hit, rev_kmer_dict, mismatches)
          else:
            continue
          
          if mismatches <= self.max_mismatches:
            matched_peptide = ''
            for i in range(0, peptide_length, self.k): # add each proteome k-mer to get the full peptide
              matched_peptide += rev_kmer_dict[kmer_hit - idx + i]

            peptide_matches.add((matched_peptide, mismatches, kmer_hit - idx))

            # in case match is 0 mismatches and self.best_match is True, return right away
            if self.best_match and not mismatches:
              return list(peptide_matches)
          
          else:
            continue
      
      except KeyError:
        continue
    
    return list(peptide_matches)

  def _find_uneven_split_matches(self, kmers, kmer_dict, rev_kmer_dict, peptide_length):
    """
    If the peptide length is NOT evenly divisible by k, perform slow search. 
    Check the associated residues for the left and right neighbors of any 
    exact matching k-mers in the proteome.
    
    EXAMPLE: "YLLDLHSYL", k=3 -> check all k-mers: 
    ["YLL", "LLD", "LDL", "DLH", "LHS", "HSY", "SYL"].
    
    If, say "DLH" is found somewhere in the proteome, check left "L" from "LDL",
    left "L" from "LLD", and "Y" from "YLL" to for any mismatches. Then do the 
    same for the residues of the right neighbors. Return any matches that are
    less than or equal to self.max_mismatches.
    """
    peptide_matches = set()
    for idx in range(0, len(kmers)): # every k-mer
      try: # try to find an exact match hit for each k-mer
        for kmer_hit in kmer_dict[kmers[idx]]:

          mismatches = 0
          mismatches = self._check_left_residues(kmers, idx, kmer_hit, rev_kmer_dict, mismatches)

          if not mismatches > self.max_mismatches:
            mismatches = self._check_right_residues(kmers, idx, kmer_hit, rev_kmer_dict, mismatches)
            
          else:
            continue
          
          if mismatches <= self.max_mismatches:
            matched_peptide = rev_kmer_dict[kmer_hit - idx] # add the first k-mer
            for i in range(self.k - 1, peptide_length):
              matched_peptide += rev_kmer_dict[kmer_hit - idx + i][-1] # add the last residue of each remaining k-mer

            matched_peptide = matched_peptide[0:peptide_length]
            peptide_matches.add((matched_peptide, mismatches, kmer_hit - idx))

            # in case match is 0 mismatches and self.best_match is True, return right away
            if self.best_match and not mismatches:
              return peptide_matches

      except KeyError:
        continue

    return peptide_matches

  def _check_left_neighbors(self, kmers, idx, kmer_hit, rev_kmer_dict, mismatches):
    """Get mismatches of left k-mer neighbors in proteome."""
    for i in range(0, idx, self.k):
      try: # get proteome k-mer from reverse dictionary and check associated query k-mer
        mismatches += hamming(rev_kmer_dict[kmer_hit + i - idx], kmers[i])
        if mismatches > self.max_mismatches:
          return 100
      except KeyError:
        return 100
    
    return mismatches

  def _check_right_neighbors(self, kmers, idx, kmer_hit, rev_kmer_dict, mismatches):
    """Get mismatches of right k-mer neighbors in proteome."""
    for i in range(self.k + idx, len(kmers), self.k):
      try: # get proteome k-mer from reverse dictionary and check associated query k-mer
        mismatches += hamming(rev_kmer_dict[kmer_hit + i - idx], kmers[i])
        if mismatches > self.max_mismatches:
          return 100
      except KeyError:
        return 100
    
    return mismatches

  def _check_left_residues(self, kmers, idx, kmer_hit, rev_kmer_dict, mismatches):
    """Get mismatches of left residues of left k-mer neighbors in proteome."""
    for i in range(0, idx):
      try: # get proteome residue from reverse dictionary and check associated query residue
        if rev_kmer_dict[kmer_hit + i - idx][0] != kmers[i][0]:
          mismatches += 1
        if mismatches > self.max_mismatches:
          return 100
      except KeyError:
        return 100
    
    return mismatches

  def _check_right_residues(self, kmers, idx, kmer_hit, rev_kmer_dict, mismatches):
    """Get mismatches of right residues of right k-mer neighbors in proteome."""
    for i in range(idx + 1, len(kmers)):
      try: # get proteome residue from reverse dictionary and check associated query residue
        if rev_kmer_dict[kmer_hit + i - idx][-1] != kmers[i][-1]:
          mismatches += 1
        if mismatches > self.max_mismatches:
          return 100
      except KeyError:
        return 100
    
    return mismatches

  def _process_mismatch_matches(self, peptide, peptide_matches, metadata_dict):
    """Extract all metadata for the mismatch matches and return as a list of tuples."""
    all_matches = []
    if not peptide_matches:
      all_matches.append((peptide,) + (np.nan,) * 13)
    else:
      for match in peptide_matches:
        match_data = (
          peptide,                                                         # query peptide
          match[0],                                                        # matched peptide
          metadata_dict[(match[2] - (match[2] % 1000000)) // 1000000][0],  # protein ID
          metadata_dict[(match[2] - (match[2] % 1000000)) // 1000000][1],  # protein name
          metadata_dict[(match[2] - (match[2] % 1000000)) // 1000000][2],  # species
          int(metadata_dict[(match[2] - (match[2] % 1000000)) // 1000000][3]), # taxon ID
          metadata_dict[(match[2] - (match[2] % 1000000)) // 1000000][4],  # gene symbol
          int(match[1]),                                                   # mismatches count
          [i+1 for i in range(len(peptide)) if peptide[i] != match[0][i]], # mutated positions
          int((match[2] % 1000000) + 1),                                   # index start
          int((match[2] % 1000000) + len(peptide)),                        # index end
          int(metadata_dict[(match[2] - (match[2] % 1000000)) // 1000000][5]), # protein existence level
          metadata_dict[(match[2] - (match[2] % 1000000)) // 1000000][6],  # sequence version
          metadata_dict[(match[2] - (match[2] % 1000000)) // 1000000][7])  # gene priority flag
        
        all_matches.append(match_data)

    return all_matches

  def best_match_search(self):
    """
    After calculating the splits we would need (starting with lowest peptide
    length in query), we can then call the mismatch_search function. We use
    the maximum # of mismatches calculated from the each length and split.
    """
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

  def discontinuous_search(self):
    """Find matches for discontinuous epitopes."""
    all_matches = []
    for discontinuous_epitope in self.discontinuous_epitopes:
      for protein_record in parse_fasta(self.proteome):
        try:
          residue_matches = sum([x[0] == protein_record.seq[x[1] - 1] for x in discontinuous_epitope])
          if residue_matches >= (len(discontinuous_epitope) - self.max_mismatches):
            metadata = extract_metadata(protein_record)
            match_data = (
              ', '.join(   # query epitope
                [x[0] + str(x[1]) for x in discontinuous_epitope]), 
              ', '.join(   # matched epitope
                [protein_record.seq[x[1] - 1] + str(x[1]) for x in discontinuous_epitope]), 
              metadata[0], # protein ID
              metadata[1], # protein name                  
              metadata[2], # species               
              metadata[3], # taxon ID              
              metadata[4], # gene symbol       
              len(discontinuous_epitope) - residue_matches, # mismatches count         
              [x[1] for x in discontinuous_epitope if x[0] != protein_record.seq[x[1] - 1]], # mutated positions
              discontinuous_epitope[0][1],  # index start           
              discontinuous_epitope[-1][1], # index end  
              metadata[5], # protein existence level         
              metadata[6], # sequence version       
              metadata[7]) # gene priority flag
            
            all_matches.append(match_data)
        
        except IndexError:
          continue
    
    return all_matches

  def _dataframe_matches(self, all_matches):
    """Return Pandas dataframe of the results."""
    df = pd.DataFrame(all_matches,
                      columns=['Query Sequence', 'Matched Sequence', 'Protein ID', 
                               'Protein Name', 'Species', 'Taxon ID', 'Gene',
                               'Mismatches', 'Mutated Positions','Index start',
                               'Index end', 'Protein Existence Level', 'Sequence Version',
                               'Gene Priority'])

    if self.best_match:
      def filter_fragments(group):
        """
        Takes out matches with 'Fragment' in the protein name if there are other
        matches without 'Fragment'.
        """
        no_fragments = group[~group['Protein Name'].str.contains('Fragment')]
        if len(no_fragments) > 0:
          return no_fragments
        else:
          return group

      # take matches with the least number of mismatches   
      idx = df.groupby(['Query Sequence'])['Mismatches'].transform(min) == df['Mismatches']
      df = df[idx]

      # take matches that are in the gene priority proteome and with the best protein
      # existence level - first fill NaN values with 0 - otherwise pandas will drop the rows
      df[['Gene Priority', 'Protein Existence Level']] = df[['Gene Priority','Protein Existence Level']].fillna(value=0)
      
      idx = df.groupby(['Query Sequence'])['Gene Priority'].transform(max) == df['Gene Priority']
      df = df[idx]

      idx = df.groupby(['Query Sequence'])['Protein Existence Level'].transform(min) == df['Protein Existence Level']
      df = df[idx]

      # filter fragments
      df = df.groupby('Query Sequence', group_keys=False).apply(filter_fragments).reset_index(drop=True)

      # sort values by protein ID and drop duplicates, guaranteeing same results 
      df.sort_values(by=['Query Sequence', 'Protein ID', 'Index start'], inplace=True)
      df.drop_duplicates(['Query Sequence'], inplace=True)

    # combine protein ID and sequence version
    df['Sequence Version'] = df['Sequence Version'].apply(lambda x: f'.{int(x)}' if not pd.isna(x) else '')   
    df['Protein ID'] = df['Protein ID'] + df['Sequence Version']
    
    # drop "Sequence Version" and "Gene Priority" columns
    df.drop(columns=['Sequence Version'], inplace=True)
    df.drop(columns=['Gene Priority'], inplace=True)

    # force integers on some columns
    int_cols = ['Taxon ID', 'Mismatches', 'Index start', 'Index end', 'Protein Existence Level']
    df[int_cols] = df[int_cols].astype('Int64')

    return df

  def _output_matches(self, df):
    """Write Pandas dataframe to format that is specified"""
    if self.output_format == 'csv':
      return df.to_csv(f'{self.preprocessed_files_path}/{self.output_name}.csv', index=False)
    elif self.output_format == 'xlsx':
      return df.to_excel(f'{self.preprocessed_files_path}/{self.output_name}.xlsx', index=False)
    elif self.output_format == 'json':
      return df.to_json(f'{self.preprocessed_files_path}/{self.output_name}.json', index=False)
    elif self.output_format == 'html':
      return df.to_html(f'{self.preprocessed_files_path}/{self.output_name}.html', index=False)



# run via command line

def parse_arguments():
  parser = argparse.ArgumentParser()

  parser.add_argument('-q', '--query', required=True)
  parser.add_argument('-p', '--proteome_file', required=True)
  parser.add_argument('-m', '--max_mismatches', type=int)
  parser.add_argument('-k', '--kmer_size', type=int)
  parser.add_argument('-P', '--preprocessed_files_path', default='.')
  parser.add_argument('-b', '--best_match', default=False, action='store_true')
  parser.add_argument('-f', '--output_format', default='csv')
  parser.add_argument('-o', '--output_name', default='')

  args = parser.parse_args()

  return args

def run():
  args = parse_arguments()

  Matcher(args.query, args.proteome_file, args.max_mismatches, 
          args.kmer_size, args.preprocessed_files_path, 
          args.best_match, args.output_format, args.output_name).match()
