import os
import argparse
import _pickle as pickle
import sqlite3
import pandas as pd

from collections import Counter, defaultdict
from Levenshtein import hamming

from .parser import parse_fasta
from .preprocessor import Preprocessor


VALID_OUTPUT_FORMATS = ['dataframe', 'csv', 'xlsx', 'json', 'html']

class Matcher(Preprocessor):
  '''
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
  '''
  def __init__(self,
               query,
               proteome,
               max_mismatches=-1,
               k=0,
               preprocessed_files_path='.',
               best_match=False,
               output_format='csv',
               output_name='',
               versioned_ids=True):
    
    # passing a Python list is possible, call results generic name 
    # if none is specified
    if type(query) == list:
      self.query = query
      if output_name:
        self.output_name = output_name
      else:
        self.output_name = 'PEPMatch_results'
    
    # parse from FASTA if not Python list
    else:
      self.query = [str(sequence.seq) for sequence in parse_fasta(query)]
      if output_name:
        self.output_name = output_name
      else:
        # make output name using query name to proteome name
        self.output_name = f"{query.replace('.fasta', '')}_to_{proteome.split('/')[-1].split('.')[0]}"

    self.lengths = sorted(set([len(peptide) for peptide in self.query]))
    self.proteome = proteome
    self.proteome_name = proteome.split('/')[-1].split('.')[0]
    self.max_mismatches = max_mismatches
    self.preprocessed_files_path = preprocessed_files_path
    self.best_match = best_match
    self.output_format = output_format
    self.versioned_ids = versioned_ids

    # select format based on # of mismatches to pass to Preprocessor
    # SQLite for exact matching - pickle for mismatching
    self.preprocess_format = 'sql' if not max_mismatches else 'pickle'
    
    # use the k that is specified if it is given in the parameters
    if k > 1:
      self.k = k
      self.k_specified = True
    else:
      self.k_specified = False

    # for exact matching, if no k is specified, use smallest length in query as k
    if not max_mismatches and not k:
      self.k = self.lengths[0]

    # for mismatching, if no k is specified, batch the peptides by length
    # then later we will use the ideal k to search them
    if max_mismatches > 0:
      self.batched_peptides = self.batch_query()
      if not k:
        self.k = 0

    # best match where no mismatches is specified means we will preprocess
    # the proteome by the smallest length, try to find exact matches,
    # then preprocess by half the length and search for more and repeat
    # until every peptide in the query has a match
    if max_mismatches == -1:
      self.ks = self.best_match_ks()
      self.k = 0

    assert k >= 0, 'Invalid k value given.'

    # if max_mismatches == -1:
    #   assert best_match, 'Number of mismatches not specified.'

    if self.output_format not in VALID_OUTPUT_FORMATS:
      raise ValueError('Invalid output format, please choose dataframe, csv, xlsx, json, or html.')

    if not all([seq.isupper() for seq in self.query]):
      raise ValueError('A peptide in the query contains a lowercase letter.')

    super().__init__(self.proteome, self.preprocess_format, self.preprocessed_files_path, self.versioned_ids)

  def batch_query(self):
    '''
    Batch peptides together by ideal k so it can be used when searching.
    '''
    batched_peptides = defaultdict(list)
    for seq in self.query:
      batched_peptides[len(seq) // (self.max_mismatches + 1)].append(seq)
    
    return dict(batched_peptides)

  def best_match_ks(self):
    '''
    For the special case where mismatching is to be done, k is not 
    specified and best match is selected, we need to get all the k values
    to preprocess the proteome multiple times. Starting with the length 
    of the smallest peptide and then halving until we get to 2.
    '''
    k = self.lengths[0]
    ks = [k]

    # halve the length and add k values until we get to 2 
    # and make sure 2 is a k in the list
    while k > 2:
      k = k // 2
      if k == 1:
        ks.append(2)
      else:
        ks.append(k)
    # mismatching function uses batched peptides, for best match, batch into one
    self.batched_peptides = {0: self.query}

    return ks

  def split_peptide(self, seq, k):
    '''
    Splits a peptide into equal sized k-mers on a rolling basis.
    Ex: k = 4, NSLFLTDLY --> ['NSLF', 'SLFL', 'LFLT', 'FLTD', 'LTDL', 'TDLY']
    '''
    kmers = []
    for i in range(0, len(seq) - k + 1):
      kmers.append(seq[i:i+k])
    return kmers

  def read_pickle_files(self):
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

  def exact_match_search(self):
    '''
    Use the preprocessed SQLite DB to perform the exact search query.
    '''
    if not os.path.isfile(os.path.join(self.preprocessed_files_path, self.proteome_name + '.db')):
      self.preprocess(self.k)

    kmers_table_name = f'{self.proteome_name}_{str(self.k)}mers'
    names_table_name = f'{self.proteome_name}_names'

    conn = sqlite3.connect(os.path.join(self.preprocessed_files_path, f'{self.proteome_name}.db'))
    c = conn.cursor()

    peptides = self.query
    all_matches_dict = {}

    for peptide in peptides:

      # skip peptide if shorter than the actual k-mer size
      if len(peptide) < self.k:
        continue

      all_matches_dict[peptide] = []
      kmers = self.split_peptide(peptide, self.k)

      # lookup each k-mer in the preprocessed proteome and subtract the offset
      # of the k-mer position in the peptide and keep track in a list
      # if the k split evenly divides the length of the peptide, the number of
      # lookups can be limited
      hit_list = []
      if len(peptide) % self.k == 0:
        for i in range(0, len(kmers), self.k):

          get_positions = f'SELECT position FROM "{kmers_table_name}" WHERE kmer = "{kmers[i]}"'
          c.execute(get_positions)
          positions_fetch = c.fetchall()

          try:
            for hit in positions_fetch:
              hit_list.append(hit[0] - i)

          except:
            continue

        # count all positions in the list and if a position appears as often as
        # the number of k-mer lookups, this will be a match
        sum_hits = Counter(hit_list)
        for hit, count in sum_hits.items():
          if count == len(peptide) // self.k:
            all_matches_dict[peptide].append(hit)

      # if the k split does not evenly divide the peptide, then all k-mers
      # need to be looked up on a rolling basis
      else:
        i = 0
        while i < len(peptide):
          try:
            get_positions = f'SELECT position FROM "{kmers_table_name}" WHERE kmer = "{kmers[i]}"'
            c.execute(get_positions)
            positions_fetch = c.fetchall()

            try:
              for hit in positions_fetch:
                hit_list.append(hit[0] - i)
            except:
              continue
            i += self.k

          # if i + k k-mer is out of range of k-mers, just check the final k-mer
          except IndexError:
            get_positions = f'SELECT position FROM "{kmers_table_name}" WHERE kmer = "{kmers[-1]}"'
            c.execute(get_positions)
            positions_fetch = c.fetchall()
            try:
              for hit in positions_fetch:
                hit_list.append(hit[0] - (len(kmers) - 1))
            except:
              continue
            i += self.k

        # if the position shows up as many times as there are k-mers, that means
        # they all agree to the same location and that is a match
        sum_hits = Counter(hit_list)
        for hit, count in sum_hits.items():
          if count == len(peptide) // self.k + 1:
            all_matches_dict[peptide].append(hit)

    all_matches = []

    # compile all matches into format used for benchmarking (comma separated)
    for peptide, matches in all_matches_dict.items():
      if matches == []:
        all_matches.append((peptide, '', '', '', '', '', '', '', '', '', '', '', ''))
      for match in matches:
        # retrieve protein IDs from the other created table
        get_protein_data = f'SELECT * FROM "{names_table_name}" WHERE protein_number = "{(match - (match % 100000)) // 100000}"'
        c.execute(get_protein_data)
        protein_data = c.fetchall()

        all_matches.append((peptide,               # query peptide
                            peptide,               # matched peptide (same as query if exact)
                            protein_data[0][1],    # taxon ID
                            protein_data[0][2],    # species name
                            protein_data[0][3],    # gene
                            protein_data[0][4],    # protein ID
                            protein_data[0][5],    # protein name
                            0,                     # 0 mismatches for exact matches
                            [],                    # mutated positions (none)
                            (match % 100000) + 1,  # index start
                            (match % 100000) + len(peptide), # index end
                            protein_data[0][6],    # protein existence level
                            protein_data[0][7]))   # gene priority binary

    c.close()
    conn.close()

    return all_matches

  def even_split_mismatching(self, kmers, kmer_dict, rev_kmer_dict, peptide_length):
    '''
    '''
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
    '''
    '''
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
        kmer_dict, names_dict = self.read_pickle_files()
        rev_kmer_dict = {i: k for k, v in kmer_dict.items() for i in v}
      except FileNotFoundError:
        self.preprocess(self.k)
        kmer_dict, names_dict = self.read_pickle_files()
        rev_kmer_dict = {i: k for k, v in kmer_dict.items() for i in v}

      for peptide in peptides:

        # split peptide into all possible k-mers with size k (self.k)
        kmers = self.split_peptide(peptide, self.k)

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
        all_matches.append((peptide, '', '', '', '', '', '', '', '', '', '', '', ''))
      else:
        for match in matches:
          all_matches.append((
            peptide,                                                         # query peptide
            match[0],                                                        # matched peptide
            names_dict[(match[2] - (match[2] % 100000)) // 100000][0],       # taxon ID
            names_dict[(match[2] - (match[2] % 100000)) // 100000][1],       # species name
            names_dict[(match[2] - (match[2] % 100000)) // 100000][2],       # gene
            names_dict[(match[2] - (match[2] % 100000)) // 100000][3],       # protein ID
            names_dict[(match[2] - (match[2] % 100000)) // 100000][4],       # protein name
            match[1],                                                        # mismatches count
            [i+1 for i in range(len(peptide)) if peptide[i] != match[0][i]], # mutated positions
            (match[2] % 100000) + 1,                                         # index start
            (match[2] % 100000) + len(peptide),                              # index end
            names_dict[(match[2] - (match[2] % 100000)) // 100000][5],       # protein existence level
            names_dict[(match[2] - (match[2] % 100000)) // 100000][6],))     # gene priority binary

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
      
  def dataframe_matches(self, all_matches):
    '''Return Pandas dataframe of the results.'''
    df = pd.DataFrame(all_matches,
                      columns=['Query Sequence', 'Matched Sequence', 'Proteome',
                               'Species', 'Gene', 'Protein ID', 'Protein Name',
                               'Mismatches', 'Mutated Positions', 'Index start', 
                               'Index end', 'Protein Existence Level', 'Gene Priority'])

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
    
    else:
      # gene priority and protein existence levels not needed if not best match
      df.drop(columns=['Gene Priority', 'Protein Existence Level'], inplace=True)
    
    # drop any columns that are entirely empty (usually gene priority column)
    df.dropna(how='all', axis=1, inplace=True)

    return df


  def output_matches(self, df):
    '''Write Pandas dataframe to format that is specified'''
    if self.output_format == 'csv':
      return df.to_csv(f'{self.output_name}.csv', index=False)
    elif self.output_format == 'xlsx':
      return df.to_excel(f'{self.output_name}.xlsx', index=False)
    elif self.output_format == 'json':
      return df.to_json(f'{self.output_name}.json')
    elif self.output_format == 'html':
      return df.to_html()

  def match(self):
    '''
    Overarching function that calls the appropriate search matching function
    based on the parameters.
    '''
    if self.max_mismatches == -1:
      df = self.dataframe_matches(self.best_match_search())

    elif self.max_mismatches > 0:
      df = self.dataframe_matches(self.mismatching_search())
    
    else:
      df = self.dataframe_matches(self.exact_match_search())

    # return a dataframe instead of outputting file if specified
    if self.output_format == 'dataframe':
      return df
    
    else:
      self.output_matches(df)


# run via command line

def parse_arguments():
  parser = argparse.ArgumentParser()

  parser.add_argument('-q', '--query', required=True)
  parser.add_argument('-p', '--proteome', required=True)
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

  Matcher(args.query, args.proteome, args.max_mismatches, 
          args.kmer_size, args.preprocessed_files_path, 
          args.best_match, args.output_format, args.output_name).match()
