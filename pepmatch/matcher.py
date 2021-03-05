#!/usr/bin/env python3

from collections import Counter
from Levenshtein import hamming
import _pickle as pickle
import pandas as pd
import sqlite3
import random

from .parser import parse_fasta
from .preprocessor import Preprocessor


splits = []
VALID_OUTPUT_FORMATS = ['csv', 'xlsx', 'json', 'html']

class Matcher(Preprocessor):
  '''
  Object class that inherits from Preprocessor. This is so queries can be
  preprocessed if the stored objects don't already exist. The class takes in the
  query (either FASTA file or Python list), proteome file, and maximum # of
  mismatches. There are multiple methods that will do exact/mismatch/best match
  searching against the proteome.

  Optional: output and output_format arguments to write results to file.
  Supported formats are "csv", "xlsx", "json", and "html".
  '''
  def __init__(self,
               query,
               proteome,
               max_mismatches,
               split=0,
               database='',
               one_match=False,
               output_df=False,
               output_format=''):

    if type(query) == list:
      self.query = query
    else:
      self.query = [str(sequence.seq) for sequence in parse_fasta(query)]

    self.lengths = sorted(list(set([len(peptide) for peptide in self.query])))
    self.proteome = proteome
    self.max_mismatches = max_mismatches
    self.database = database
    self.one_match = one_match
    self.output_df = output_df
    self.output_format = output_format

    # use sql format for exact matches
    if split == 0:
      if max_mismatches == 0:
        self.split = min(self.lengths)
        self.preprocess_format = 'sql'

      # use pickle format for mismatching
      # calculate optimal k-splits by the query passed
      elif max_mismatches > 0:
        self.splits = self.mismatching_splits(self.lengths)
        self.split = self.splits[0]
        self.preprocess_format = 'pickle'
      elif max_mismatches == -1:
        self.splits = self.best_match_splits(min(self.lengths))
        self.split = self.splits[0]
        self.preprocess_format = 'pickle'
    else:
      self.split = split
      if self.database != '':
        self.preprocess_format = 'sql'
      else:
        self.preprocess_format = 'pickle'

    if self.output_format != '' and self.output_format not in VALID_OUTPUT_FORMATS:
      raise ValueError('Invalid output format, please choose csv, xlsx, json, or html.')

    super().__init__(self.proteome, self.split, self.preprocess_format, self.database, True)

  def split_peptide(self, seq, k):
    '''
    Splits a peptide into equal sized k-mers on a rolling basis.
    Ex: k = 4, NSLFLTDLY --> ['NSLF', 'SLFL', 'LFLT', 'FLTD', 'LTDL', 'TDLY']
    '''
    kmers = []
    for i in range(0, len(seq) - k + 1):
      kmer = seq[i:i+k]
      kmers.append(kmer)
    return kmers

  def read_pickle_files(self):
    '''
    Read in the already created pickle files for each dictionary in the
    preprocessing step.
    '''
    name = self.proteome.replace('.fa', '').replace('.fasta', '')
    with open(name + '_' + str(self.split) + 'mers' + '.pickle', 'rb') as f:
      kmer_dict = pickle.load(f)
    with open(name + '_names.pickle', 'rb') as f:
      names_dict = pickle.load(f)
    return kmer_dict, names_dict

  def sql_exact_match(self):
    '''
    Use the preprocessed SQLite DB to perform the exact search query.
    '''
    proteome_name = self.proteome.split('/')[-1].split('.')[0]
    kmers_table_name = proteome_name + '_' + str(self.split) + 'mers'
    names_table_name = proteome_name + '_names'

    conn = sqlite3.connect(self.database)
    c = conn.cursor()

    peptides = self.query
    all_matches_dict = {}

    for peptide in peptides:
      all_matches_dict[peptide] = []
      kmers = self.split_peptide(peptide, self.split)

      # lookup each k-mer in the preprocessed proteome and subtract the offset
      # of the k-mer position in the peptide and keep track in a list
      # if the k split evenly divides the length of the peptide, the number of
      # lookups can be limited
      hit_list = []
      if len(peptide) % self.split == 0:
        for i in range(0, len(kmers), self.split):
          get_positions = 'SELECT position FROM "{kmer_table}" WHERE kmer = "{actual_kmer}"'.format(kmer_table = kmers_table_name, actual_kmer = kmers[i])
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
          if count == len(peptide) // self.split:
            all_matches_dict[peptide].append(hit)

      # if the k split does not evenly divide the peptide, then all k-mers
      # need to be looked up on a rolling basis
      else:
        i = 0
        while i < len(peptide):
          try:
            get_positions = 'SELECT position FROM "{kmer_table}" WHERE kmer = "{actual_kmer}"'.format(kmer_table = kmers_table_name, actual_kmer = kmers[i])
            c.execute(get_positions)
            positions_fetch = c.fetchall()

            try:
              for hit in positions_fetch:
                hit_list.append(hit[0] - i)
            except:
              continue
            i += self.split

          # if i + k k-mer is out of range of k-mers, just check the final k-mer
          except IndexError:
            get_positions = 'SELECT position FROM "{kmer_table}" WHERE kmer = "{actual_kmer}"'.format(kmer_table = kmers_table_name, actual_kmer = kmers[-1])
            c.execute(get_positions)
            positions_fetch = c.fetchall()
            try:
              for hit in positions_fetch:
                hit_list.append(hit[0] - (len(kmers) - 1))
            except:
              continue
            i += self.split

        # if the position shows up as many times as there are k-mers, that means
        # they all agree to the same location and that is a match
        sum_hits = Counter(hit_list)
        for hit, count in sum_hits.items():
          if count == len(peptide) // self.split + 1:
            all_matches_dict[peptide].append(hit)

    all_matches = []

    # compile all matches into format used for benchmarking (comma separated)
    for peptide, matches in all_matches_dict.items():
      if matches == []:
        all_matches.append((peptide, '', '', '', '', '', ''))
      for match in matches:
        # retrieve protein IDs from the other created table
        get_protein_data = 'SELECT * FROM "{names_table}" WHERE protein_number = "{protein_number}"'.format(
                  names_table = names_table_name, protein_number = (match - (match % 100000)) // 100000)
        c.execute(get_protein_data)
        protein_data = c.fetchall()

        all_matches.append((peptide, 
                            peptide, 
                            protein_data[0][1],
                            (match % 100000),
                            (match % 100000) + len(peptide), 
                            protein_data[0][2], 
                            protein_data[0][3]))

    c.close()
    conn.close()

    return all_matches

  def mismatching_splits(self, lengths):
    '''
    There is an ideal k split for each peptide and max # of mismatches.
    Some queries have many lengths, so we can calculate all the optimal k-size
    splits we would need for optimal searching here.
    '''
    splits = set()
    for length in lengths:
      if length // (self.max_mismatches + 1) not in [0,1]:
        splits.add(length // (self.max_mismatches + 1))

    splits = sorted(list(splits), reverse = True)

    return splits

  def mismatching(self):
    '''
    Searches a preprocessed proteome for all matches of a given query of
    peptides in FASTA format up to a number of specified mismatches that was
    initialized by the class.
    '''
    all_matches_dict = {}

    peptides = self.query
    try:
      kmer_dict, names_dict = self.read_pickle_files()
      rev_kmer_dict = {i: k for k, v in kmer_dict.items() for i in v}
    except FileNotFoundError:
      self.preprocess()
      kmer_dict, names_dict = self.read_pickle_files()
      rev_kmer_dict = {i: k for k, v in kmer_dict.items() for i in v}

    for peptide in peptides:
      # record matches in a set so as to not duplicate matches
      matches = set()

      # split peptide into k-mers (rolling basis)
      kmers = self.split_peptide(peptide, self.split)

      # if the peptide length has an even split of k, perform faster search
      if len(peptide) % self.split == 0:
        for i in range(0, len(kmers), self.split):
          try:

            # find each hit for each k-mer
            for hit in kmer_dict[kmers[i]]:

              mismatches = 0

              # if the k-mer is found in the middle or end, check the neighboring
              # k-mers to the left
              for j in range(0, i, self.split):
                try:

                  # use reverse dictionary to retrive k-mers for Hamming distance
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
              for k in range(i+self.split, len(kmers), self.split):
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
                  for s in range(0, len(peptide), self.split):
                    matched_peptide += rev_kmer_dict[hit-i+s]
                except KeyError:
                  continue

                matches.add((matched_peptide, mismatches, hit - i))

          # if nothing is found, you can check the next k-mer, since it can still be a match
          except KeyError:
            continue

      # if the peptide length does NOT HAVE an even split of k, perform slower search (rolling split)
      else:
        for i in range(0, len(kmers)):

          try:
            # find each hit for each k-mer
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
                  for s in range(0, len(peptide), self.split):
                    matched_peptide += rev_kmer_dict[hit-i+s]

                except KeyError:
                  for r in range(1, len(peptide) % self.split + 1):
                    matched_peptide += rev_kmer_dict[hit-i+s-(self.split - r)][-1]

                matched_peptide = matched_peptide[0:len(peptide)]
                matches.add((matched_peptide, mismatches, hit - i))

          # if nothing is found, you can check the next k-mer, since it can still be a match
          except KeyError:
            continue

      all_matches_dict[peptide] = list(matches)

    all_matches = []

    # compile all matches to output in tuples in the following format:
    # (peptide, matched peptide, protein matched in, index, # of mismatches)
    for peptide, matches in all_matches_dict.items():
      if matches == []:
        all_matches.append((peptide, '', '', '', '', '', ''))
      else:
        for match in matches:
          all_matches.append((
            peptide,
            match[0],
            names_dict[(match[2] - (match[2] % 100000)) // 100000][0],
            match[1],
            [i+1 for i in range(len(peptide)) if peptide[i] != match[0][i]],
            match[2] % 100000,
            (match[2] % 100000) + len(peptide)))

    return all_matches

  def best_match_splits(self, length):
    '''
    For best matching, we start with the lowest length then divide by half.
    We continue to divide subsequent k by half until we reach k = 2. Then we
    use these splits for multiple lookups until we find a match for every
    peptide.
    '''
    splits.append(length)
    if length > 3:
      return self.best_match_splits(length // 2)
    elif length == 2:
      return splits
    else:
      splits.append(2)
      return splits

  def best_match(self):
    '''
    After calculating the splits we would need (starting with lowest peptide
    length in query), we can then call the mismatching function. We use
    the maximum # of mismatches calculated from the each length and split.
    '''
    all_matches = []

    for split in self.splits:
      self.split = split

      # calculate maximum possible # of mismatches that can guaranteed to be found
      # using the lowest length and k-split
      self.max_mismatches = (min(self.lengths) // self.split)

      matches = self.mismatching()

      # separate out peptides that did not match into a new query
      self.query = []
      for match in matches:
        if match[1] != '':
          all_matches.append(match)
        else:
          self.query.append(match[0])

    # once we reach k = 2, we will just increase the # of mismatches until we find some
    # match for each peptide in the query
    while self.query != []:
      self.max_mismatches += 1

      matches = self.mismatching()

      # separate out peptides that did not match into a new query
      self.query = []
      for match in matches:
        if match[1] != '':
          all_matches.append(match)
        else:
          self.query.append(match[0])

    return all_matches

  def match(self):
    '''
    Overarching function that calls the appropriate search matching function
    based on the parameters.
    '''
    if self.max_mismatches == 0:
      all_matches = self.sql_exact_match()

    # if we specify the # of mismatches, we can the use the optimal k-split
    # given the lengths and # of mismatches
    elif self.max_mismatches > 0:
      if self.split:
        all_matches = self.mismatching()
      else:
        all_matches = []
        query = self.query
        for split in self.splits:
          self.split = split
          self.query = [peptide for peptide in query if (len(peptide) // (self.max_mismatches + 1)) == split]
          all_matches.extend(self.mismatching())

    elif self.max_mismatches == -1:
      all_matches = self.best_match()

    else:
      raise ValueError('Invalid input of mismatches.')

    # option to output results to a specified format, default is CSV
    if self.output_df:
      if self.max_mismatches == 0:
        df = self.dataframe_exact_matches(all_matches)
      else:
        df = self.dataframe_mismatch_matches(all_matches)
      
      if self.output_format == '':
        return df
      else:
        self.output_matches(df)

    return all_matches

  def dataframe_exact_matches(self, all_matches):
    '''Return Pandas dataframe of the results.'''
    df = pd.DataFrame(all_matches,
                      columns=['Peptide Sequence',
                               'Matched Peptide',
                               'Protein ID',
                               'Index start',
                               'Index end',
                               'Protein Evidence Level',
                               'Gene Priority'])
    
    if self.one_match:
      try:
        idx = df.groupby(['Peptide Sequence'])['Protein Evidence Level'].transform('min') == df['Protein Evidence Level']
        df = df[idx]

        idx = df.groupby(['Peptide Sequence'])['Gene Priority'].transform('max') == df['Gene Priority']
        df = df[idx]

        df.drop_duplicates(['Peptide Sequence'], inplace=True)

      except KeyError:
        df.drop_duplicates(['Peptide Sequence'], inplace=True)

    return df 

  def dataframe_mismatch_matches(self, all_matches):
    '''Return Pandas dataframe of the results.'''
    df = pd.DataFrame(all_matches,
                      columns=['Peptide Sequence',
                               'Matched Peptide',
                               'Protein ID',
                               'Mismatches',
                               'Mutated Positions',
                               'Index start',
                               'Index end'])

    if self.one_match:
      df.drop_duplicates(['Peptide Sequence'], inplace=True)

    return df 

  def output_matches(self, df):
    '''Write Pandas dataframe to format that is specified'''
    results_id = ''.join(random.choice('0123456789ABCDEFabcdef') for i in range(6))
    if self.output_format == 'csv':
      return df.to_csv('PEPMatch_results_' + results_id + '.csv')
    elif self.output_format == 'xlsx':
      return df.to_excel('PEPMatch_results_' + results_id + '.xlsx')
    elif self.output_format == 'json':
      return df.to_json('PEPMatch_results_' + results_id + '.json')
    elif self.output_format == 'html':
      return df.to_html()