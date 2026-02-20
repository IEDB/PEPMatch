import os
import polars as pl
from pathlib import Path
from Bio import SeqIO
from ._rs import rs_preprocess, rs_match
from .helpers import output_matches

VALID_OUTPUT_FORMATS = ['dataframe', 'csv', 'tsv', 'xlsx', 'json']
FASTA_EXTENSIONS = {
  '.fasta', '.fas', '.fa', '.fna', '.ffn', '.faa', '.mpfa', '.frn'
}

class Matcher:
  """Searches query peptides against a preprocessed proteome index
  and returns matches as a Polars DataFrame or output file."""

  def __init__(
    self,
    query,
    proteome_file,
    max_mismatches=0,
    k=5,
    preprocessed_files_path='.',
    best_match=False,
    output_format='dataframe',
    output_name='',
    sequence_version=True,
  ):
    if k < 2:
      raise ValueError('k must be >= 2.')
    if output_format not in VALID_OUTPUT_FORMATS:
      raise ValueError(
        f'Invalid output format. Choose from: {VALID_OUTPUT_FORMATS}'
      )

    self.proteome_file = str(proteome_file)
    self.proteome_name = str(proteome_file).split('/')[-1].split('.')[0]
    self.max_mismatches = max_mismatches
    self.k = k
    self.preprocessed_files_path = preprocessed_files_path
    self.best_match = best_match
    self.output_format = output_format
    self.sequence_version = sequence_version
    self.output_name = output_name or 'PEPMatch_results'

    self.query = self._parse_query(query)
    if not self.query:
      raise ValueError('Query is empty.')

  def _parse_query(self, query):
    if isinstance(query, list):
      return [(str(i + 1), seq.upper()) for i, seq in enumerate(query)]

    path = Path(query)
    if not path.is_file():
      raise FileNotFoundError(f'Query file not found: {query}')

    if path.suffix.lower() in FASTA_EXTENSIONS:
      return [
        (record.id, str(record.seq).upper())
        for record in SeqIO.parse(str(path), 'fasta')
      ]
    elif path.suffix.lower() == '.txt':
      with open(path) as f:
        lines = [line.strip() for line in f if line.strip()]
      return [(str(i + 1), seq.upper()) for i, seq in enumerate(lines)]

    raise ValueError(
      f'Unsupported query format: {path.suffix}. '
      f'Use a Python list, .txt file, or FASTA file for the query.'
    )

  def match(self):
    pepidx_path = os.path.join(
      self.preprocessed_files_path, f'{self.proteome_name}_{self.k}mers.pepidx'
    )
    if not os.path.isfile(pepidx_path):
      rs_preprocess(self.proteome_file, self.k, pepidx_path)

    results = rs_match(pepidx_path, self.query, self.k, self.max_mismatches)
    df = self._to_dataframe(results)

    if self.output_format == 'dataframe':
      return df
    else:
      output_matches(df, self.output_format, self.output_name)

  def _to_dataframe(self, results):
    schema = {
      'Query ID': pl.Utf8,
      'Query Sequence': pl.Utf8,
      'Matched Sequence': pl.Utf8,
      'Protein ID': pl.Utf8,
      'Protein Name': pl.Utf8,
      'Species': pl.Utf8,
      'Taxon ID': pl.Utf8,
      'Gene': pl.Utf8,
      'Mismatches': pl.Utf8,
      'Mutated Positions': pl.Utf8,
      'Index start': pl.Utf8,
      'Index end': pl.Utf8,
      'Protein Existence Level': pl.Utf8,
      'Sequence Version': pl.Utf8,
      'Gene Priority': pl.Utf8,
      'SwissProt Reviewed': pl.Utf8,
    }

    if not results:
      return pl.DataFrame(schema=schema)

    df = pl.DataFrame(results, schema=schema, orient='row')

    df = df.with_columns(
      pl.when(pl.col(col) == '').then(None).otherwise(pl.col(col)).alias(col)
      for col in df.columns
    )

    df = df.with_columns([
      pl.col('Mismatches').cast(pl.Int64),
      pl.col('Index start').cast(pl.Int64),
      pl.col('Index end').cast(pl.Int64),
      pl.col('Protein Existence Level').cast(pl.Int64),
      pl.col('Sequence Version').cast(pl.Int64),
      pl.col('Gene Priority').cast(pl.Int64),
      pl.col('SwissProt Reviewed').cast(pl.Int64).cast(pl.Boolean),
    ])

    if self.sequence_version:
      df = df.with_columns(
        pl.when(pl.col('Sequence Version').is_not_null())
          .then(pl.col('Protein ID') + '.' + pl.col('Sequence Version').cast(pl.Utf8))
          .otherwise(pl.col('Protein ID'))
          .alias('Protein ID')
      )

    df = df.drop('Sequence Version')
    return df
