import pandas as pd
import multiprocessing as mp
from Bio import SeqIO

from .matcher import Matcher


class ParallelMatcher(Matcher):
  def __init__(self, n_jobs=1, **kwargs):
    self.kwargs = kwargs
    self.n_jobs = n_jobs

    self.query = kwargs['query']
    self.preprocessed_files_path = kwargs['preprocessed_files_path']
    self.output_format = kwargs['output_format']

    self.output_name = kwargs['output_name']
    if self.output_name == '':
      self.output_name = 'PEPMatch_results'


  def _split_query(self):
    query = self.query
    if not isinstance(query, list):
      query = [str(record.seq) for record in list(SeqIO.parse(self.query, 'fasta'))]

    query_chunks = []
    chunk_size = len(query) // self.n_jobs

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
    elif self.output_format == 'xlsx':
      return df.to_excel(f'{path}.xlsx', index=False)
    elif self.output_format == 'json':
      return df.to_json(f'{path}.json', index=False)
    elif self.output_format == 'html':
      return df.to_html(f'{path}.html', index=False)