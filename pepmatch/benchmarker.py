from .preprocessor import Preprocessor
from .parallel_match import ParallelMatcher


class Benchmarker(ParallelMatcher):
  """Benchmarker class for PEPMatch. This will be the class that is called by
  benchmarking.py that calls the necessary methods to run the benchmarking."""

  def __init__(
    self, benchmark: str, query: str, proteome: str, lengths: list, max_mismatches: int,
    algorithm_parameters: dict
  ):
    self.benchmark = benchmark
    self.query = query
    self.proteome = proteome
    self.max_mismatches = max_mismatches
    self.lengths = lengths
    self.algorithm_parameters = algorithm_parameters
    
    super().__init__(
      query=query,
      proteome_file=proteome,
      max_mismatches=max_mismatches,
      best_match=True if max_mismatches == -1 else False,
      output_format='dataframe',
      n_jobs=5
    )

  def __str__(self):
    return 'PEPMatch'

  def preprocess_query(self) -> None:
    """No query preprocessing, raise TypeError"""
    raise TypeError(self.__str__() + ' does not preprocess queries.\n')

  def preprocess_proteome(self) -> None:
    """Preprocess proteome once or multiple times for each split calculated."""

    preprocessor = Preprocessor(self.proteome)

    if self.benchmark == 'mhc_ligands':
      preprocessor.sql_proteome(k=9)
    elif self.benchmark == 'coronavirus':
      for i in range(2, 6):
        preprocessor.pickle_proteome(k=i)
    elif self.benchmark == 'neoepitopes':
      preprocessor.pickle_proteome(k=3)
    elif self.benchmark == 'milk':
      preprocessor.sql_proteome(k=15)
      for i in [2, 3, 7]:
        preprocessor.pickle_proteome(k=i)
    

  def search(self) -> list:
    """Call overarching match function. Then convert results into the standard format
    needed to calculate accuracy.
    """
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
