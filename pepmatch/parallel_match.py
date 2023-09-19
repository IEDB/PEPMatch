import pandas as pd
import multiprocessing as mp

from .helpers import parse_fasta
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

    # set # of jobs to # of queries if # of jobs is greater than # of queries
    self.n_jobs = len(query) if self.n_jobs > len(query) else self.n_jobs

    query_chunks = []
    chunk_size = max(1, len(query) // self.n_jobs)

    for i in range(self.n_jobs):
      start_idx = i * chunk_size
      end_idx = start_idx + chunk_size

      if i == self.n_jobs - 1:
        end_idx = len(query)
      
      query_chunks.append(query[start_idx:end_idx])

    return query_chunks


  def match(self):
    query_chunks = self._split_query()

    kwargs_without_query = {k: v for k, v in self.kwargs.items() if k != 'query'}
    kwargs_without_query['output_format'] = 'dataframe'

    with mp.Pool(self.n_jobs) as pool:
      dfs = pool.map(
        Matcher.match,
        [Matcher(query=query, **kwargs_without_query) for query in query_chunks]
      )

    df = pd.concat(dfs, ignore_index=True)

    if self.output_format == 'dataframe':
      return df
    
    else:
      self._output_matches(df)


  def _output_matches(self, df: pd.DataFrame) -> None:
    path = f'{self.preprocessed_files_path}/{self.output_name}'
    if self.output_format == 'csv':
      return df.to_csv(f'{path}.csv', index=False)
    elif self.output_format == 'tsv':
      return df.to_csv(f'{path}.tsv', sep='\t', index=False)
    elif self.output_format == 'xlsx':
      return df.to_excel(f'{path}.xlsx', index=False)
    elif self.output_format == 'json':
      return df.to_json(f'{path}.json', index=False)
    elif self.output_format == 'html':
      return df.to_html(f'{path}.html', index=False)