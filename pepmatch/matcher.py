import os
import re
import argparse
import _pickle as pickle
import sqlite3
import pandas as pd
import numpy as np

from collections import Counter, defaultdict
from Levenshtein import hamming

from .helpers import parse_fasta, split_sequence
from .preprocessor import Preprocessor


VALID_OUTPUT_FORMATS = ['dataframe', 'csv', 'xlsx', 'json', 'html']

class Matcher(Preprocessor):
  """
  Object class that inherits from Preprocessor. This is so queries can be
  preprocessed if the stored objects don't already exist. The class takes in the
  query (either FASTA file or Python list), proteome file, and maximum # of
  mismatches. There are multiple methods that will do exact/mismatch/best match
  searching against the proteome.

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

    super().__init__(self.proteome, self.preprocessed_files_path)

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

  def _process_peptide_exact_matches(self, peptide, peptide_matches, cursor, metadata_table_name):
    all_matches = []
    if not peptide_matches:
      all_matches.append((peptide,) + (np.nan,) * 12)
    else:
      for match in peptide_matches:
        protein_number = (match - (match % 1000000)) // 1000000
        query = f'SELECT * FROM "{metadata_table_name}" WHERE protein_number = "{protein_number}"'
        cursor.execute(query)
        protein_data = cursor.fetchone()
        match_data = (peptide,                          # query peptide
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

  def exact_match_search(self):
    if not os.path.isfile(os.path.join(self.preprocessed_files_path, self.proteome_name + '.db')):
      self.preprocess(self.k)

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
      target_kmers = all_kmers if self.k == len(peptide) else [all_kmers[i] for i in range(0, len(all_kmers), self.k)] + [all_kmers[-1]]
      
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

      processed_matches = self._process_peptide_exact_matches(peptide, peptide_matches, cursor, metadata_table_name)
      all_matches.extend(processed_matches)

    cursor.close()
    conn.close()

    return all_matches

  def even_split_mismatching(self, kmers, kmer_dict, rev_kmer_dict, peptide_length):
    # record matches in a set so as to not duplicate matches
    matches = set()

    for i in range(0, len(kmers), self.k):

      # find each hit for each k-mer
      try:
        for hit in kmer_dict[kmers[i]]:

          mismatches = 0

          # if the k-mer is found in the middle or end, check the neighboring
          # k-mers to the left
          for j in range(0, i, self.k):
            
            # use reverse dictionary to retrive k-mers for Hamming distance
            try:
              mismatches += hamming(rev_kmer_dict[hit+j-i], kmers[j])
              
              # if mismatches ever reach threshold, break out of loop
              if mismatches >= self.max_mismatches + 1:
                break

            # if first k-mer finds nothing, set mismatches to 100 to disqualify this
            # peptide from matching with this area
            except KeyError:
              mismatches = 100

          # if the k-mer is found in the middle or end, check the neighbors
          # k-mers to the right
          for k in range(i+self.k, len(kmers), self.k):
            try:

              # use reverse dictionary to retrive k-mers for Hamming distance
              mismatches += hamming(rev_kmer_dict[hit+k-i], kmers[k])

              # if mismatches ever reach threshold, break out of loop
              if mismatches >= self.max_mismatches + 1:
                break

            # if last k-mer finds nothing, set mismatches to 100 to disqualify this
            # peptide from matching with this area
            except KeyError:
              mismatches = 100

          # if the mismatches that were calculated is less than threshold
          # for all neighbors, then it's a match
          if mismatches < self.max_mismatches + 1:
            matched_peptide = ''

            try:
              for s in range(0, peptide_length, self.k):
                matched_peptide += rev_kmer_dict[hit-i+s]
            except KeyError:
              continue

            matches.add((matched_peptide, mismatches, hit - i))

            if self.best_match and not mismatches:
              return matches

      # if nothing is found, you can check the next k-mer, since it can still be a match
      except KeyError:
        continue

    return matches

  def uneven_split_mismatching(self, kmers, kmer_dict, rev_kmer_dict, peptide_length):
    # record matches in a set so as to not duplicate matches
    matches = set()

    for i in range(0, len(kmers)):

      # find each hit for each k-mer
      try:
        for hit in kmer_dict[kmers[i]]:

          mismatches = 0

          # if the k-mer is found in the middle or end, check the neighbors
          # k-mers to the left
          for j in range(0, i,):
            try:
              # use reverse dictionary to retrive k-mers and just check the
              # very first letter since it's a rolling split
              if rev_kmer_dict[hit+j-i][0] != kmers[j][0]:
                mismatches += 1

              # if mismatches ever reach threshold, break out of loop
              if mismatches >= self.max_mismatches + 1:
                break

            # if first k-mer finds nothing, set mismatches to 100 to disqualify this
            # peptide from matching with this area
            except KeyError:
              mismatches = 100

          # if the k-mer is found in the middle or end, check the neighbors
          # k-mers to the right
          for k in range(i+1, len(kmers),):
            try:
              # use reverse dictionary to retrive k-mers and just check the
              # very last letter since it's a rolling split
              if rev_kmer_dict[hit+k-i][-1] != kmers[k][-1]:
                mismatches += 1

              # if mismatches ever reach threshold, break out of loop
              if mismatches >= self.max_mismatches + 1:
                break

            # if last k-mer finds nothing, set mismatches to 100 to disqualify this
            # peptide from matching with this area
            except KeyError:
              mismatches = 100

          if mismatches >= self.max_mismatches + 1:
            continue

          # if the mismatches that were calculated is less than threshold
          # for all neighbors, then it's a match
          if mismatches < self.max_mismatches + 1:
            matched_peptide = ''
            try:
              for s in range(0, peptide_length, self.k):
                matched_peptide += rev_kmer_dict[hit-i+s]

            except KeyError:
              for r in range(1, peptide_length % self.k + 1):
                matched_peptide += rev_kmer_dict[hit-i+s-(self.k - r)][-1]

            matched_peptide = matched_peptide[0:peptide_length]
            matches.add((matched_peptide, mismatches, hit - i))

            if self.best_match and not mismatches:
              return matches

      # if nothing is found, you can check the next k-mer, since there can still be a match
      except KeyError:
        continue

    return matches

  def _read_pickle_files(self):
    '''
    Read in the already created pickle files for each dictionary in the
    preprocessing step.
    '''
    with open(os.path.join(self.preprocessed_files_path, 
      f'{self.proteome_name}_{str(self.k)}mers.pickle'), 'rb') as f:
      kmer_dict = pickle.load(f)

    with open(os.path.join(self.preprocessed_files_path, 
      f'{self.proteome_name}_names.pickle'), 'rb') as f:
      names_dict = pickle.load(f)

    return kmer_dict, names_dict

  def mismatching_search(self):
    '''
    Searches a preprocessed proteome for all matches of a given query of
    peptides in FASTA format up to a number of specified mismatches that was
    initialized by the class.
    '''
    all_matches_dict = {}

    for k, peptides in self.batched_peptides.items():
      if not self.k_specified:
        self.k = k

      # try reading in the preprocessed pickle files, if they don't exist
      # we'll need to preprocess the proteome using self.k
      try:
        kmer_dict, names_dict = self._read_pickle_files()
        rev_kmer_dict = {i: k for k, v in kmer_dict.items() for i in v}
      except FileNotFoundError:
        self.preprocess(self.k)
        kmer_dict, names_dict = self._read_pickle_files()
        rev_kmer_dict = {i: k for k, v in kmer_dict.items() for i in v}

      for peptide in peptides:

        # split peptide into all possible k-mers with size k (self.k)
        kmers = split_sequence(peptide, self.k)

        # if the peptide length has an even split of k, perform faster search
        if len(peptide) % self.k == 0:
          matches = self.even_split_mismatching(kmers, kmer_dict, rev_kmer_dict, len(peptide))

        # if the peptide length does NOT HAVE an even split of k, perform slower search (rolling split)
        else:
          matches = self.uneven_split_mismatching(kmers, kmer_dict, rev_kmer_dict, len(peptide))

        all_matches_dict[peptide] = list(matches)


    all_matches = []

    # compile all matches to output in tuples in the following format:
    # (peptide, matched peptide, protein matched in, index, # of mismatches)
    for peptide, matches in all_matches_dict.items():
      if matches == []:
        all_matches.append((peptide, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan))
      else:
        for match in matches:
          all_matches.append((
            peptide,                                                         # query peptide
            match[0],                                                        # matched peptide
            names_dict[(match[2] - (match[2] % 1000000)) // 1000000][0],       # taxon ID
            names_dict[(match[2] - (match[2] % 1000000)) // 1000000][1],       # species name
            names_dict[(match[2] - (match[2] % 1000000)) // 1000000][2],       # gene
            names_dict[(match[2] - (match[2] % 1000000)) // 1000000][3],       # protein ID
            names_dict[(match[2] - (match[2] % 1000000)) // 1000000][4],       # protein name
            match[1],                                                        # mismatches count
            [i+1 for i in range(len(peptide)) if peptide[i] != match[0][i]], # mutated positions
            (match[2] % 1000000) + 1,                                         # index start
            (match[2] % 1000000) + len(peptide),                              # index end
            names_dict[(match[2] - (match[2] % 1000000)) // 1000000][5],       # protein existence level
            names_dict[(match[2] - (match[2] % 1000000)) // 1000000][6],))     # gene priority binary

    return all_matches

  def best_match_search(self):
    '''
    After calculating the splits we would need (starting with lowest peptide
    length in query), we can then call the mismatching function. We use
    the maximum # of mismatches calculated from the each length and split.
    '''

    all_matches = []
    # once we reach k = 2, we will just increase the # of mismatches until we find some
    # match for each peptide in the query
    for k in self.ks:
      self.k = k
      self.k_specified = True

      # calculate maximum possible # of mismatches that can guaranteed to be found
      # using the lowest length and k-split
      self.max_mismatches = (self.lengths[0] // self.k) - 1

      matches = self.mismatching_search()

      # separate out peptides that did not match into a new query
      self.query = []
      for match in matches:
        if match[1]:
          all_matches.append(match)
        else:
          self.query.append(match[0])

      if not self.query:
        return all_matches
      
      self.batched_peptides = {0: self.query}

    return all_matches

  def get_protein_metadata(self, protein_record):
    """
    Extracts protein metadata from a protein record in FASTA format.
    This is exclusively used for discontinuous search as the metadata
    for exact match and mismatch searching is already stored in the
    preprocessed files.
    """
    protein_id = protein_record.id.split('|')[1] if '|' in protein_record.id else protein_record.id
    
    # use regex to get all other data from the UniProt FASTA header
    try:
      taxon_id = int(re.search('OX=(.*?) ', protein_record.description).group(1))
    except AttributeError:
      taxon_id = None

    try:
      species = re.search('OS=(.*) OX=', protein_record.description).group(1)
    except AttributeError:
      species = None

    try:
      gene = re.search('GN=(.*?) ', protein_record.description).group(1)
    except AttributeError:
      gene = None

    try:
      protein_name = re.search(' (.*) OS', protein_record.description).group(1)
    except AttributeError:
      protein_name = None

    try:
      pe_level = int(re.search('PE=(.*?) ', protein_record.description).group(1))
    except AttributeError:
      pe_level = 0

    return (taxon_id, species, gene, protein_id, protein_name, pe_level)

  def discontinuous_search(self):
    """Find matches for discontinuous epitopes."""
    all_matches = []
    for discontinuous_epitope in self.discontinuous_epitopes:
      for protein_record in parse_fasta(self.proteome):
        try:
          residue_matches = sum([x[0] == protein_record.seq[x[1] - 1] for x in discontinuous_epitope])
          if residue_matches >= (len(discontinuous_epitope) - self.max_mismatches):
            protein_metadata = self.get_protein_metadata(protein_record)
            all_matches.append((', '.join([x[0] + str(x[1]) for x in discontinuous_epitope]), 
                                ', '.join([protein_record.seq[x[1] - 1] + str(x[1]) for x in discontinuous_epitope]),
                                protein_metadata[0],
                                protein_metadata[1],
                                protein_metadata[2],
                                protein_metadata[3],
                                protein_metadata[4],
                                len(discontinuous_epitope) - residue_matches,
                                [x[1] for x in discontinuous_epitope if x[0] != protein_record.seq[x[1] - 1]],
                                discontinuous_epitope[0][1],
                                discontinuous_epitope[-1][1],
                                protein_metadata[5],
                                np.nan))
        except IndexError:
          continue
    
    return all_matches

  def dataframe_matches(self, all_matches):
    '''Return Pandas dataframe of the results.'''
    df = pd.DataFrame(all_matches,
                      columns=['Query Sequence', 'Matched Sequence', 'Protein ID', 
                               'Protein Name', 'Species', 'Taxon ID', 'Gene',
                               'Mismatches', 'Mutated Positions','Index start',
                               'Index end', 'Protein Existence Level', 'Sequence Version',
                               'Gene Priority'])

    if self.best_match:
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

      # sort values by protein ID and drop duplicates, guaranteeing same results 
      df.sort_values(by='Protein ID', inplace=True)
      df.drop_duplicates(['Query Sequence'], inplace=True)

    # combine protein ID and sequence version
    df['Sequence Version'] = df['Sequence Version'].apply(lambda x: f'.{int(x)}' if not pd.isna(x) else '')
    df['Protein ID'] = df['Protein ID'] + df['Sequence Version']
    
    # drop "Sequence Version" and "Gene Priority" columns
    df.drop(columns=['Sequence Version'], inplace=True)
    df.drop(columns=['Gene Priority'], inplace=True)

    return df

  def output_matches(self, df):
    '''Write Pandas dataframe to format that is specified'''
    if self.output_format == 'csv':
      return df.to_csv(f'{self.preprocessed_files_path}/{self.output_name}.csv', index=False)
    elif self.output_format == 'xlsx':
      return df.to_excel(f'{self.preprocessed_files_path}/{self.output_name}.xlsx', index=False)
    elif self.output_format == 'json':
      return df.to_json(f'{self.preprocessed_files_path}/{self.output_name}.json', index=False)
    elif self.output_format == 'html':
      return df.to_html(f'{self.preprocessed_files_path}/{self.output_name}.html', index=False)

  def match(self):
    '''
    Overarching function that calls the appropriate search matching function
    based on the parameters.
    '''
    if self.query:
      if self.max_mismatches == -1:
        query_df = self.dataframe_matches(self.best_match_search())

      elif self.max_mismatches > 0:
        query_df = self.dataframe_matches(self.mismatching_search())
      
      else:
        query_df = self.dataframe_matches(self.exact_match_search())
    else:
      query_df = pd.DataFrame()

    # search for discontinuous epitopes if they exist
    if self.discontinuous_epitopes:
      discontinuous_df = self.dataframe_matches(self.discontinuous_search())
    else:
      discontinuous_df = pd.DataFrame()

    # combine the dataframes if both query and discontinuous epitopes exist
    df = pd.concat([query_df, discontinuous_df], ignore_index=True)

    # return a dataframe instead of outputting file if specified
    if self.output_format == 'dataframe':
      return df
    
    else:
      self.output_matches(df)


# run via command line

def parse_arguments():
  parser = argparse.ArgumentParser()

  parser.add_argument('-q', '--query', required=True)
  parser.add_argument('-p', '--proteome_file', required=True)
  parser.add_argument('-m', '--max_mismatches', type=int, required=True)
  parser.add_argument('-k', '--kmer_size', type=int, required=True)
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
