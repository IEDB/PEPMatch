from .preprocessor import Preprocessor
from .parallel_match import ParallelMatcher


class Benchmarker(ParallelMatcher):
  """Benchmarker class for PEPMatch. This will be the class that is called by
  benchmarking.py that calls the necessary methods to run the benchmarking."""

  def __init__(
    self, benchmark: str, query: str, proteome: str, lengths: list, max_mismatches: int,
    method_parameters: dict
  ):
    self.benchmark = benchmark
    self.query = query
    self.proteome = proteome
    self.max_mismatches = max_mismatches
    self.lengths = lengths
    self.method_parameters = method_parameters
    
    super().__init__(
      query=query,
      proteome_file=proteome,
      max_mismatches=max_mismatches,
      best_match=True if max_mismatches == -1 else False,
      output_format='dataframe',
      sequence_version=False,
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
    """Call overarching match function which returns a dataframe of matches."""
    return self.match()
