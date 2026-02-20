from .preprocessor import Preprocessor
from .matcher import Matcher

class Benchmarker:
  def __init__(self, benchmark, query, proteome, lengths, max_mismatches, method_parameters):
    self.query = str(query)
    self.proteome = str(proteome)
    self.lengths = lengths
    self.max_mismatches = max_mismatches

    if max_mismatches == -1:
      self.best_match = True
      self.mm = 0
      self.k = 0
    elif max_mismatches == 0:
      self.best_match = False
      self.mm = 0
      self.k = min(lengths)
    else:
      self.best_match = False
      self.mm = max_mismatches
      if benchmark == 'coronavirus':
        self.k = 2
      elif benchmark == 'neoepitopes':
        self.k = 3
      else:
        self.k = min(lengths) // (max_mismatches + 1)
        if self.k < 2:
          self.k = 2

  def preprocess_proteome(self):
    if self.best_match:
      for k in [15, 7, 3, 2]:
        Preprocessor(self.proteome).preprocess(k=k)
    else:
      Preprocessor(self.proteome).preprocess(k=self.k)

  def preprocess_query(self):
    raise TypeError('PEPMatch does not preprocess queries.')

  def search(self):
    df = Matcher(
      query=self.query,
      proteome_file=self.proteome,
      max_mismatches=self.mm,
      k=self.k,
      best_match=self.best_match,
      output_format='dataframe',
      sequence_version=False
    ).match()

    return df.select(
      ['Query Sequence', 'Matched Sequence', 'Protein ID', 'Index start']
    ).to_pandas()
