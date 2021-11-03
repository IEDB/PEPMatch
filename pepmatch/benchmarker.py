from .matcher import Matcher

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
