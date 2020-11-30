#!/usr/bin/env python3

from Bio import SeqIO
from collections import Counter
from Levenshtein import hamming
import _pickle as pickle
import pandas as pd
import os, glob, sqlite3

splits = []

def parse_fasta(file):
  '''Return a parsed Biopython SeqRecord object from a FASTA file.'''
  return SeqIO.parse(file, 'fasta')

class Preprocessor(object):
  '''
  Object class that takes in a proteome FASTA file, k-size split, and format which
  then preprocesses the proteome. 

  With the preprocess method, it will break the proteome into equal size k-mers and 
  map them to locations within the individual proteins. The mapped keys and values 
  will then be stored in either pickle files or a SQLite database.

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
    name = self.proteome.split('.')[0]
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
    name = self.proteome.split('.')[0].split('/')[-1]
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

    c.execute('CREATE INDEX IF NOT EXISTS "{id}" ON "{k}"(kmer)'.format(id = kmers_table + '_id', k = kmers_table))
    c.execute('CREATE INDEX IF NOT EXISTS "{id}" ON "{n}"(protein_number)'.format(id = names_table + '_id', n = names_table))

    conn.commit()
    c.close()
    conn.close()

  def preprocess(self):
    '''
    Processing the given proteome, by splitting each protein into k-mers and assigning
    a unique index to each unique k-mer within each protein. This is done by assigning 
    a number to each protein and for each k-mer, multiplying the protein number by 100,000
    and adding the index position of the index within the protein. This guarantees a unique
    index for each and every possible k-mer. Also, each protein # assigned is also mapped
    to the protein ID to be read back later after searching.
    '''
    proteome = parse_fasta(self.proteome)
    kmer_dict = {}
    names_dict = {}
    protein_count = 1

    for protein in proteome:

      # split protein into k-mers
      kmers = self.split_protein(str(protein.seq), self.split)

      for i in range(len(kmers)):

        # if k-mer is already found in the dict, add position to value list
        if kmers[i] in kmer_dict.keys():
          kmer_dict[kmers[i]].append(protein_count * 100000 + i)
        
        # if k-mer is not in dict, add it along with first list entry
        else:
          kmer_dict[kmers[i]] = [protein_count * 100000 + i]

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

    # send the created dicts to be dumped into pickle files
    if self.preprocess_format == 'pickle':
      self.pickle_proteome(kmer_dict, names_dict)
    
    # or send dicts to function to be logged into SQLite DB
    elif self.preprocess_format == 'sql':
      self.sql_proteome(kmer_dict, names_dict)
    
    else:
      raise AssertionError('Unexpected value of preprocessing format', self.preprocess_format)

    return kmer_dict, names_dict


VALID_OUTPUT_FORMATS = ['csv', 'xlsx', 'json', 'html']

class Matcher(object):
  '''
  Object class that inherits from Preprocessor. This is so queries can be preprocessed if
  the stored objects don't already exist. The class takes in the query (either FASTA file 
  or Python list), proteome file, and maximum # of mismatches. There are multiple methods
  that will do exact/mismatch/best match searching against the proteome.
  Optional: output and output_format arguments to write results to file. Supported formats
  are "csv", "xlsx", "json", and "html".
  '''
  def __init__(self, query, proteome, max_mismatches, database='', output = False, output_format = 'csv'):
    if type(query) == list:
      self.query = query
    else:
      self.query = [str(sequence.seq) for sequence in parse_fasta(query)]

    self.proteome = proteome
    self.lengths = sorted(list(set([len(peptide) for peptide in self.query])))
    
    self.max_mismatches = max_mismatches
    self.database = database

    self.output = output
    self.output_format = output_format

    # use sql format for exact matches
    if max_mismatches == 0:
      self.split = min(self.lengths)
      self.split = 5
      self.preprocess_format = 'sql'
    
    # use pickle format for mismatching and calculate optimal k-splits by the query passed
    elif max_mismatches > 0:
      self.splits = self.mismatching_splits(self.lengths)
      self.split = self.splits[0]
      self.preprocess_format = 'pickle'
    elif max_mismatches == -1:
      self.splits = self.best_match_splits(min(self.lengths))
      self.split = self.splits[0]
      self.preprocess_format = 'pickle'
    
    if self.output_format not in VALID_OUTPUT_FORMATS:
      raise ValueError('Invalid output format, please choose csv, xlsx, json, or html.')

    #super().__init__(self.proteome, self.split, self.preprocess_format)

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
    name = self.proteome.split('.')[0]

    # read in k-mer dictionary pickle file
    with open(name + '_kmers' + '_' + str(self.split) + '.pickle', 'rb') as f:
      kmer_dict = pickle.load(f)

    # read in protein names dictionary pickle file
    with open(name + '_names.pickle', 'rb') as f:
      names_dict = pickle.load(f)

    return kmer_dict, names_dict

  def sql_exact_match(self):
    '''
    Use the preprocessed SQLite DB to perform the exact search query.
    '''
    proteome_name = self.proteome.replace('.fa', '').split('/')[-1]

    # make table names based on proteome and split size
    kmers_table_name = proteome_name + '_kmers' + '_' + str(self.split)
    names_table_name = proteome_name + '_names'

    conn = sqlite3.connect(self.database)
    c = conn.cursor()

    peptides = self.query

    all_matches_dict = {}

    for peptide in peptides:

      # keep track of all matches and split peptide into k-mers
      all_matches_dict[peptide] = []
      kmers = self.split_peptide(peptide, self.split)
      
      hit_list = []

      if len(peptide) % self.split == 0:

        # lookup each k-kmer and add positions with offset of k-mer position in peptide
        for i in range(0, len(kmers), self.split):
          get_positions = 'SELECT position FROM "{kmer_table}" WHERE kmer = "{actual_kmer}"'.format(kmer_table = kmers_table_name, actual_kmer = kmers[i])
          c.execute(get_positions)
          positions_fetch = c.fetchall()
          try:
            for hit in positions_fetch:
              hit_list.append(hit[0] - i)
          except:
            continue

        # find where all k-mer positions agree on 1 location from their corresponding offset
        # this will be a match if they all agree
        sum_hits = Counter(hit_list)
        for hit, count in sum_hits.items():
          if count == len(peptide) // self.split:
            all_matches_dict[peptide].append(hit)

      else:
        # start with first k-mer (0) look up then check i + k k-mer from there
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
        
        # find where all k-mer positions agree on 1 location from their corresponding offset
        # this will be a match if they all agree
        sum_hits = Counter(hit_list)
        for hit, count in sum_hits.items():
          if count == len(peptide) // self.split + 1:
            all_matches_dict[peptide].append(hit)

    all_matches = []

    # compile all matches into format used for benchmarking (comma separated)
    for peptide, matches in all_matches_dict.items():
      if matches == []:
        all_matches.append((peptide, '', '', '', ''))
      for match in matches:
        # retrieve protein IDs from the other created table
        get_protein_name = 'SELECT protein_ID FROM "{names_table}" WHERE protein_number = "{protein_number}"'.format(
                  names_table = names_table_name, protein_number = (match - (match % 100000)) // 100000)
        c.execute(get_protein_name)
        protein_ID = c.fetchall()

        all_matches.append((peptide, peptide, protein_ID[0][0], 0, match % 100000))

    c.close()
    conn.close()

    return all_matches

  def exact_match(self):
    '''
    Searches a preprocessed proteome for exact matches of a given
    query file in FASTA format. 
    '''

    # either load the pickle files or run the query via a SQLite DB
    if self.preprocess_format == 'pickle':
      try:
        kmer_dict, names_dict = self.read_pickle_files()
      except FileNotFoundError:
        self.preprocess()
    elif self.preprocess_format == 'sql':
      return self.sql_exact_match()

    # load in query peptides and pickle files of preprocessed proteome
    peptides = self.query

    # keep track of all the matches found within a dictionary
    all_matches_dict = {}

    for peptide in peptides:
      
      # create list of matches entry for each peptide 
      all_matches_dict[str(peptide.seq)] = []
      
      # split peptide into k-mers (rolling basis)
      kmers = self.split_peptide(str(peptide.seq), self.split)

      hit_list = []

      for offset, kmer in enumerate(kmers):
        try:

          # find each hit for each k-mer
          hits = kmer_dict[kmer]
          for hit in hits:
            
            # keep track of each index position for each hit
            # using offset to find k-mers that agree with start
            # index of the 1st k-mer
            hit_list.append(hit - offset)
        
        except KeyError:
          break

      # use Counter function to count all k-mers that agree with 1st position
      sum_hits = Counter(hit_list)
      
      for hit, count in sum_hits.items():

        # if all k-mers agree with 1st position, then it's a match
        if count == len(peptide.seq) - self.split + 1:
          all_matches_dict[str(peptide.seq)].append(hit) # add to dict entry

    all_matches = []
    
    # compile all matches to output in tuples in the following format:
    # (peptide, protein matched in, index, # of mismatches)
    for peptide, matches in all_matches_dict.items():
      if matches == []:
        all_matches.append((peptide, '', '', '', ''))
      for match in matches:
        all_matches.append((peptide, peptide, names_dict[(match - (match % 100000)) // 100000], 0, match % 100000))

    return all_matches

  def mismatching_splits(self, lengths):
    '''
    There is an ideal k split for each peptide and max # of mismatches.
    Some queries have many lengths, so we can calculate all the optimal k-size splits
    we would need for optimal searching here.
    '''
    splits = set()
    for length in lengths:
      if length // (self.max_mismatches + 1) not in [0,1]:
        splits.add(length // (self.max_mismatches + 1))
    
    splits = sorted(list(splits), reverse = True)

    return splits

  def mismatching(self):
    '''
    Searches a preprocessed proteome for all matches of a given query of peptides in 
    FASTA format up to a number of specified mismatches that was initialized by the class.
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
        all_matches.append((peptide, '', '', '', ''))
      else:
        for match in matches:
          all_matches.append((
            peptide,
            match[0],
            names_dict[(match[2] - (match[2] % 100000)) // 100000],
            match[1],
            match[2] % 100000,
            ))

    return all_matches

  def best_match_splits(self, length):
    '''
    For best matching, we start with the lowest length then divide by half.
    We continue to divide subsequent k by half until we reach k = 2. Then we 
    use these splits for multiple lookups until we find a match for every peptide.
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
    After calculating the splits we would need (starting with lowest peptide length 
    in query), we can then call the mismatching function. We the maximum # of mismatches
    calculated from the each length and split.
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
      all_matches = self.exact_match()
    
    # if we specify the # of mismatches, we can the use the optimal k-split
    # given the lengths and # of mismatches
    elif self.max_mismatches > 0:
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
    if self.output:
      df = self.dataframe_matches(all_matches)
      self.output_matches(df)

    return all_matches

  def dataframe_matches(self, all_matches):
    '''Return Pandas dataframe of the results.'''  
    return pd.DataFrame(all_matches, 
      columns=['Peptide Sequence', 'Matched Peptide', 'Protein ID', 'Mismatches', 'Index start'])

  def output_matches(self, df):
    '''Write Pandas dataframe to format that is specified'''
    if self.output_format == 'csv':
      return df.to_csv('PEPMatch_results.csv')
    elif self.output_format == 'xlsx':
      return df.to_excel('PEPMatch_results.xlsx')
    elif self.output_format == 'json':
      return df.to_json('PEPMatch_results.json')
    elif self.output_format == 'html':
      return df.to_html()


class Benchmarker(Matcher):
  '''
  Object used for benchmarking the PEPMatch code for the various applications.
  It uses three methods: 2 for preprocessing and 1 for searching as the benchmarking
  code is structured. PEPMatch does not do any query preprocessing so it raises
  a TypeError which is excepted in the benchmarking code.
  Inherits from Matcher object.
  '''
  def __init__(self, query, proteome, lengths, max_mismatches, algorithm_parameters):
    self.query = query
    self.proteome = proteome
    self.max_mismatches = max_mismatches
    self.lengths = lengths
    self.algorithm_parameters = algorithm_parameters
    
    super().__init__(query, proteome, max_mismatches)

  def __str__(self):
    return 'PEPMatch'

  def preprocess_query(self):
    '''No query preprocessing, raise TypeError'''
    raise TypeError(self.__str__() + ' does not preprocess queries.\n')

  def preprocess_proteome(self):
    '''Preprocess proteome once or multiple times for each split calculated.'''
    if self.max_mismatches == 0:
      self.preprocess()
    else:
      for split in self.splits:
        self.split = split
        self.preprocess()

  def search(self):
    '''
    Call overarching match function. Then convert results into the standard format
    needed to calculate accuracy.
    '''
    matches = self.match()

    all_matches = []
    for match in matches:
      if match[1] == '':
        continue
      match_string = ''
      for i in range(0, len(match)):

        if i == len(match) - 1:
          match_string += str(match[i])
        else:
          match_string += str(match[i]) + ','

      all_matches.append(match_string)

    # remove files after benchmarking as the only purpose is to extract benchmarks
    try:
      os.remove(self.database)
    except FileNotFoundError:
      for file in glob.glob(os.path.dirname(self.proteome) + '/*.pickle'):
        os.remove(file)

    return all_matches
