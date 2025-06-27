from __future__ import annotations 
import polars as pl
import multiprocessing as mp
from .helpers import parse_fasta, output_matches
from .matcher import Matcher

class ParallelMatcher(Matcher):
  def __init__(self, n_jobs=1, **kwargs):
    self.kwargs = kwargs
    self.n_jobs = n_jobs
    self.query = kwargs.get('query', [])
    self.output_format = kwargs.get('output_format', 'csv')
    self.preprocessed_files_path = kwargs.get('preprocessed_files_path', '.')
    self.output_name = kwargs.get('output_name', '')
    if self.output_name == '':
      self.output_name = 'PEPMatch_results'
    assert self.output_format in ['csv', 'tsv', 'xlsx', 'json', 'html', 'dataframe']

  def _split_query(self):
    query = self.query
    if not isinstance(query, list):
      query = [str(record.seq) for record in parse_fasta(self.query)]
    
    # Do not create more jobs than there are queries
    if not query:
      self.n_jobs = 0
      return []
    self.n_jobs = min(self.n_jobs, len(query))
    
    query_chunks = []
    chunk_size = max(1, len(query) // self.n_jobs)
    for i in range(self.n_jobs):
      start_idx = i * chunk_size
      end_idx = start_idx + chunk_size
      if i == self.n_jobs - 1:
        end_idx = len(query)
      if start_idx < end_idx: # Ensure we don't create empty chunks
        query_chunks.append(query[start_idx:end_idx])
    return query_chunks

  @staticmethod
  def _run_match(matcher_instance: Matcher) -> pl.DataFrame:
    """Helper function for multiprocessing to call the instance method."""
    return matcher_instance.match()

  def match(self) -> pl.DataFrame:
    query_chunks = self._split_query()
    if not query_chunks:
      dummy_matcher = Matcher(query=[], proteome_file=self.kwargs.get('proteome_file'))
      return dummy_matcher.match()

    kwargs_without_query = {k: v for k, v in self.kwargs.items() if k != 'query'}
    kwargs_without_query['output_format'] = 'dataframe'

    matcher_instances = [
      Matcher(query=query, **kwargs_without_query) for query in query_chunks
    ]

    with mp.Pool(self.n_jobs) as pool:
      dfs = pool.map(ParallelMatcher._run_match, matcher_instances)

    df = pl.concat(dfs) 

    if self.output_format == 'dataframe':
      return df
    else:
      output_matches(df, self.output_format, self.output_name)

  def _output_matches(self, df: pl.DataFrame) -> None:
    path = self.output_name.__str__()
    if not path.lower().endswith(f".{self.output_format}"):
      path += f".{self.output_format}"
        
    if self.output_format == 'csv':
      df.write_csv(path)
    elif self.output_format == 'tsv':
      df.write_csv(path, separator='\t')
    elif self.output_format == 'xlsx':
      df.write_excel(path)
    elif self.output_format == 'json':
      df.write_json(path)
