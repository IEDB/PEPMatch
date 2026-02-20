import os
import polars as pl
from pathlib import Path
from Bio import SeqIO
from ._rs import rs_preprocess, rs_match, rs_discontinuous

VALID_OUTPUT_FORMATS = ['dataframe', 'csv', 'tsv', 'xlsx', 'json']
FASTA_EXTENSIONS = {
  '.fasta', '.fas', '.fa', '.fna', '.ffn', '.faa', '.mpfa', '.frn'
}

def output_matches(df: pl.DataFrame, output_format: str, output_name: str) -> None:
  path = output_name.__str__()
  if not path.lower().endswith(f".{output_format}"):
    path += f".{output_format}"
  if output_format == 'csv':
    df.write_csv(path)
  elif output_format == 'tsv':
    df.write_csv(path, separator='\t')
  elif output_format == 'xlsx':
    df.write_excel(path)
  elif output_format == 'json':
    df.write_json(path)

class Matcher:
  """Searches query peptides against a preprocessed proteome index
  and returns matches as a Polars DataFrame or output file."""

  def __init__(
    self,
    query,
    proteome_file,
    max_mismatches=0,
    k=0,
    preprocessed_files_path='.',
    best_match=False,
    output_format='dataframe',
    output_name='',
    sequence_version=True,
  ):

    if best_match and k == 0:
      self.k = 0
      self.k_specified = False
    elif k == 0:
      self.k = 5
      self.k_specified = False
    else:
      self.k = k
      self.k_specified = True

    if k != 0 and k < 2:
      raise ValueError('k must be >= 2.')

    if output_format not in VALID_OUTPUT_FORMATS:
      raise ValueError(
        f'Invalid output format. Choose from: {VALID_OUTPUT_FORMATS}'
      )

    self.proteome_file = str(proteome_file)
    self.proteome_name = str(proteome_file).split('/')[-1].split('.')[0]
    self.max_mismatches = max_mismatches
    self.preprocessed_files_path = preprocessed_files_path
    self.best_match = best_match
    self.output_format = output_format
    self.sequence_version = sequence_version
    self.output_name = output_name or 'PEPMatch_results'

    self.query = self._parse_query(query)
    self.discontinuous_epitopes = self._find_discontinuous_epitopes()
    self.query = self._clean_query()
    if not self.query and not self.discontinuous_epitopes:
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

  def _find_discontinuous_epitopes(self):
    discontinuous_epitopes = {}
    for query_id, peptide in self.query:
      try:
        epitope = [(x[0], int(x[1:])) for x in peptide.split(', ')]
        discontinuous_epitopes[query_id] = epitope
      except (ValueError, IndexError):
        continue
    return discontinuous_epitopes

  def _clean_query(self):
    discontinuous_ids = set(self.discontinuous_epitopes.keys())
    return [
      (qid, seq) for qid, seq in self.query
      if qid not in discontinuous_ids
    ]

  def match(self):
    linear_df = pl.DataFrame()
    discontinuous_df = pl.DataFrame()

    if self.query:
      if self.best_match and self.k_specified:
        results = self._search(self.k, self.max_mismatches)
        linear_df = self._best_match_filter(self._to_dataframe(results))
      elif self.best_match:
        linear_df = self.best_match_search()
      else:
        results = self._search(self.k, self.max_mismatches)
        linear_df = self._to_dataframe(results)

    if self.discontinuous_epitopes:
      pepidx_path = self._pepidx_path(2)
      if not os.path.isfile(pepidx_path):
        rs_preprocess(self.proteome_file, 2, pepidx_path)
      epitopes = [
        (qid, residues) for qid, residues in self.discontinuous_epitopes.items()
      ]
      results = rs_discontinuous(pepidx_path, epitopes, self.max_mismatches)
      discontinuous_df = self._to_dataframe(results)

    dfs = [d for d in [linear_df, discontinuous_df] if d.height > 0]
    df = pl.concat(dfs, how="vertical") if dfs else linear_df

    if self.output_format == 'dataframe':
      return df
    else:
      self.output_matches(df, self.output_format, self.output_name)

  def _pepidx_path(self, k):
    return os.path.join(
      self.preprocessed_files_path, f'{self.proteome_name}_{k}mers.pepidx'
    )

  def _search(self, k, max_mismatches, peptides=None):
    pepidx_path = self._pepidx_path(k)
    if not os.path.isfile(pepidx_path):
      print(f"Preprocessing {self.proteome_name} with k={k}...")
      rs_preprocess(self.proteome_file, k, pepidx_path)
    query = peptides or self.query
    print(f"Searching {len(query)} peptides against {self.proteome_name} (k={k}, max_mismatches={max_mismatches})...")
    return rs_match(pepidx_path, query, k, max_mismatches)

  def best_match_search(self):
    peptides_remaining = self.query.copy()
    all_results = []

    initial_k = min(len(seq) for _, seq in peptides_remaining)
    ks = [initial_k]
    while initial_k > 2:
      initial_k //= 2
      ks.append(max(initial_k, 2))
    ks = sorted(set(ks), reverse=True)

    for k in ks:
      if not peptides_remaining:
        break

      min_len = min(len(seq) for _, seq in peptides_remaining)
      max_mm = (min_len // k) - 1
      if max_mm < 0:
        max_mm = 0

      results = self._search(k, max_mm, peptides_remaining)

      matched_ids = set()
      for row in results:
        if row[2] != '':
          all_results.append(row)
          matched_ids.add(row[0])

      peptides_remaining = [
        (qid, seq) for qid, seq in peptides_remaining if qid not in matched_ids
      ]
      print(f"  -> k={k}, mismatches<={max_mm}: {len(peptides_remaining)} remaining")

    if peptides_remaining:
      max_mm = min(len(seq) for _, seq in peptides_remaining) // 2

      while peptides_remaining:
        shortest_len = min(len(seq) for _, seq in peptides_remaining)
        if max_mm >= shortest_len:
          break

        max_mm += 1
        results = self._search(2, max_mm, peptides_remaining)

        matched_ids = set()
        for row in results:
          if row[2] != '':
            all_results.append(row)
            matched_ids.add(row[0])

        peptides_remaining = [
          (qid, seq) for qid, seq in peptides_remaining if qid not in matched_ids
        ]
        print(f"  -> k=2, mismatches<={max_mm}: {len(peptides_remaining)} remaining")

    if peptides_remaining:
      for qid, seq in peptides_remaining:
        all_results.append([qid, seq] + [''] * 14)

    df = self._to_dataframe(all_results)
    return self._best_match_filter(df)

  def _best_match_filter(self, df):
    matched_df = df.filter(pl.col("Matched Sequence").is_not_null())
    unmatched_df = df.filter(pl.col("Matched Sequence").is_null())

    if matched_df.height > 0:
      matched_df = matched_df.with_columns(
        pl.col("Mismatches").min().over("Query ID").alias("_min_mm")
      ).filter(
        pl.col("Mismatches") == pl.col("_min_mm")
      ).drop("_min_mm")

      matched_df = matched_df.with_columns(
        pl.col("Gene Priority").max().over("Query ID").alias("_max_gp")
      ).filter(
        pl.col("Gene Priority") == pl.col("_max_gp")
      ).drop("_max_gp")

      matched_df = (
        matched_df.with_columns(
          (~pl.col("Protein ID").str.contains("-")).cast(pl.Int8).alias("_is_canonical")
        ).with_columns(
          pl.col("_is_canonical").max().over("Query ID").alias("_max_canonical")
        ).filter(
          pl.col("_is_canonical") == pl.col("_max_canonical")
        ).drop(["_is_canonical", "_max_canonical"])
      )

      matched_df = matched_df.with_columns(
        pl.col("SwissProt Reviewed").cast(pl.Int8).max().over("Query ID").alias("_max_reviewed")
      ).filter(
        pl.col("SwissProt Reviewed").cast(pl.Int8) == pl.col("_max_reviewed")
      ).drop("_max_reviewed")

      matched_df = matched_df.with_columns(
        pl.col("Protein Existence Level").min().over("Query ID").alias("_min_pe")
      ).filter(
        pl.col("Protein Existence Level") == pl.col("_min_pe")
      ).drop("_min_pe")

      matched_df = (
        matched_df.with_columns(
          (~pl.col("Protein Name").str.contains("Fragment")).alias("_not_fragment")
        ).with_columns(
          pl.col("_not_fragment").any().over("Query ID").alias("_has_non_fragment")
        ).filter(
          (pl.col("_not_fragment") & pl.col("_has_non_fragment")) |
          (~pl.col("_has_non_fragment"))
        ).drop(["_not_fragment", "_has_non_fragment"])
      )

      matched_df = matched_df.unique(subset=["Query ID"], keep="first")

    return pl.concat([matched_df, unmatched_df], how="vertical")

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
