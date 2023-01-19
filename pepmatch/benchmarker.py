from .matcher import Matcher

import os, glob

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
    
    super().__init__(query, proteome, max_mismatches, output_format=algorithm_parameters['output_format'], versioned_ids=False)

  def __str__(self):
    return 'PEPMatch'

  def preprocess_query(self):
    '''No query preprocessing, raise TypeError'''
    raise TypeError(self.__str__() + ' does not preprocess queries.\n')

  def preprocess_proteome(self):
    '''Preprocess proteome once or multiple times for each split calculated.'''
    if self.max_mismatches == -1:
      for k in self.best_match_ks():
        self.preprocess(k)
    else:
      for k in self.batch_query().keys():
        self.preprocess(k)

  def search(self):
    '''
    Call overarching match function. Then convert results into the standard format
    needed to calculate accuracy.
    '''
    if self.max_mismatches == -1:
      self.k_specified = True
      self.k = 2
      self.max_mismatches = 7

    matches = self.match()

    all_matches = []
    for i, match in matches.iterrows():
      match_string = ''
      match_string += match['Query Sequence'] + ','
      match_string += match['Matched Sequence'] + ','
      match_string += match['Protein ID'] + ','
      match_string += str(match['Mismatches']) + ','
      try:
        match_string += str(match['Index start'] - 1)
      except TypeError:
        match_string += ''
        
      all_matches.append(match_string)
    
    return all_matches
