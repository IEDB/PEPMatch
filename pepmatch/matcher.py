import os
import polars as pl
from itertools import combinations
from pathlib import Path
from Bio import SeqIO
from ._rs import rs_preprocess, rs_match, rs_discontinuous, rs_metadata, rs_match_counts, rs_indel_match

VALID_OUTPUT_FORMATS = ['dataframe', 'polars', 'pandas', 'csv', 'tsv', 'xlsx', 'json']
DATAFRAME_OUTPUT_FORMATS = frozenset({'dataframe', 'polars', 'pandas'})
FASTA_EXTENSIONS = {
  '.fasta', '.fas', '.fa', '.fna', '.ffn', '.faa', '.mpfa', '.frn'
}

def output_matches(df: pl.DataFrame, output_format: str, output_name: str) -> None:
  """Write the results to disk in `output_format`, appending the extension if absent."""
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


def _return_dataframe(df: pl.DataFrame, output_format: str):
  """Return the finished frame as Polars ('dataframe' is an alias for it) or as pandas.

  The pandas conversion needs pyarrow, so pandas-only installs fall back to building
  the frame column by column."""
  if output_format in ('dataframe', 'polars'):
    return df
  try:
    import pandas as pd
  except ImportError as e:
    raise ImportError(
      "output_format='pandas' requires pandas. "
      "Install with: pip install 'pepmatch[pandas]'"
    ) from e
  try:
    return df.to_pandas()
  except ModuleNotFoundError:
    return pd.DataFrame({col: df[col].to_list() for col in df.columns})


def _indel_placements(query, matched):
  """Every valid placement of the edits between `query` and `matched`, each a tuple of
  (1-based query position, residue) per edit; a repeat makes placements equivalent.
  Query-terminal deletions and boundary insertions are barred, and an inserted
  residue's position counts the query residues preceding it."""
  L, M = len(query), len(matched)
  n = abs(M - L)
  if n == 0:
    return []
  out = []
  if M < L:
    for combo in combinations(range(1, L - 1), n):
      if ''.join(query[i] for i in range(L) if i not in combo) == matched:
        out.append(tuple((i + 1, query[i]) for i in combo))
  else:
    for combo in combinations(range(M), n):
      if ''.join(matched[i] for i in range(M) if i not in combo) != query:
        continue
      placement = []
      for k in combo:
        p = sum(1 for x in range(k) if x not in combo) + 1
        if p == 1 or p == L + 1:
          placement = None
          break
        placement.append((p, matched[k]))
      if placement:
        out.append(tuple(placement))
  return out


def _indel_edits(query, matched):
  """Entries for the Indel Positions column as (kind, residues, low, high).

  One edit gives one position range. Two edits at the same site collapse to a chunk,
  ranged only across starts that remove identical residues (in a periodic repeat the
  missing pair changes along the range); independent edits each get their own range."""
  placements = _indel_placements(query, matched)
  if not placements:
    return []
  kind = 'd' if len(matched) < len(query) else 'i'

  if len(placements[0]) == 1:
    positions = sorted(p for ((p, _),) in placements)
    residue = placements[0][0][1]
    return [(kind, residue, positions[0], positions[-1])]

  chunks = sorted({
    (p1, r1 + r2) for (p1, r1), (p2, r2) in placements
    if (p2 == p1 + 1 if kind == 'd' else p2 == p1)
  })
  if chunks:
    entries = []
    for start, residues in chunks:
      if entries and entries[-1][1] == residues and entries[-1][3] + 1 == start:
        k, r, low, _ = entries[-1]
        entries[-1] = (k, r, low, start)
      else:
        entries.append((kind, residues, start, start))
    return entries

  entries = []
  for slot in (0, 1):
    positions = sorted({pl[slot][0] for pl in placements})
    entries.append((kind, placements[0][slot][1], positions[0], positions[-1]))
  return entries


def format_indel_positions(query, matched):
  """Render the Indel Positions column, e.g. 'd: A[6]', 'd: AA[2,4]', 'd: C[3], d: E[5]',
  or '[]' for an exact match. A 1-based range collapses when the placement is unambiguous."""
  edits = _indel_edits(query, matched)
  if not edits:
    return '[]'
  return ', '.join(
    f'{kind}: {residues}[{low}]' if low == high else f'{kind}: {residues}[{low},{high}]'
    for kind, residues, low, high in edits
  )


class Matcher:
  """Searches query peptides against a preprocessed proteome index
  and returns matches as a DataFrame (Polars by default, or pandas) or output file."""

  def __init__(
    self,
    query,
    proteome_file,
    max_mismatches=0,
    max_indels=0,
    k=0,
    preprocessed_files_path='.',
    best_match=False,
    output_format='dataframe',
    output_name='',
    sequence_version=True,
    counts_only=False,
  ):

    if k == 0:
      self.k = 0
      self.k_specified = False
    else:
      self.k = k
      self.k_specified = True

    if k != 0 and k < 2:
      raise ValueError('k must be >= 2.')

    if max_indels > 2:
      raise ValueError('max_indels > 2 is not yet supported. Only max_indels<=2 has been validated.')

    if max_indels > 0 and max_mismatches > 0:
      raise ValueError('max_indels and max_mismatches are mutually exclusive.')

    if max_indels > 0 and best_match:
      raise ValueError('max_indels and best_match are not yet supported together.')

    if max_indels > 0 and counts_only:
      raise ValueError('max_indels and counts_only are not yet supported together.')

    self.max_indels = max_indels

    if output_format not in VALID_OUTPUT_FORMATS:
      raise ValueError(
        f'Invalid output format. Choose from: {VALID_OUTPUT_FORMATS}'
      )

    if counts_only and best_match:
      raise ValueError('counts_only is not supported with best_match.')

    self.proteome_file = str(proteome_file)
    self.proteome_name = Path(proteome_file).stem
    self.max_mismatches = max_mismatches
    self.preprocessed_files_path = preprocessed_files_path
    self.best_match = best_match
    self.output_format = output_format
    self.sequence_version = sequence_version
    self.counts_only = counts_only
    self.output_name = output_name or 'PEPMatch_results'

    self.query = self._parse_query(query)
    self.discontinuous_epitopes = self._find_discontinuous_epitopes()
    if self.max_indels > 0 and self.discontinuous_epitopes:
      raise ValueError('max_indels and discontinuous epitopes are not supported together.')
    self.query = self._clean_query()
    if not self.query and not self.discontinuous_epitopes:
      raise ValueError('Query is empty.')

    if self.k_specified and self.query:
      shortest = min(len(seq) for _, seq in self.query)
      if self.k > shortest:
        raise ValueError(
          f'k ({self.k}) cannot exceed the shortest query peptide length ({shortest}); '
          f'shorter peptides would be silently skipped.'
        )

  def _parse_query(self, query):
    """Normalize the query into (id, sequence) pairs from a list, .txt file, or FASTA."""
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
    """Pick out the queries written as residue/position lists, e.g. 'N1, A3, L5'."""
    discontinuous_epitopes = {}
    for query_id, peptide in self.query:
      try:
        epitope = [(x[0], int(x[1:])) for x in peptide.split(', ')]
        discontinuous_epitopes[query_id] = epitope
      except (ValueError, IndexError):
        continue
    return discontinuous_epitopes

  def _clean_query(self):
    """Drop the discontinuous epitopes from the linear query set."""
    discontinuous_ids = set(self.discontinuous_epitopes.keys())
    return [
      (qid, seq) for qid, seq in self.query
      if qid not in discontinuous_ids
    ]

  def match(self):
    """Run the search, returning a DataFrame or writing the results to `output_name`."""
    linear_df = pl.DataFrame()
    discontinuous_df = pl.DataFrame()

    if self.counts_only:
      k = self.k if self.k_specified else self._auto_k(self.max_mismatches)
      self._recall_warning(k, self.max_mismatches)
      df = self._counts_to_dataframe(self._search_counts(k, self.max_mismatches))
      if self.output_format in DATAFRAME_OUTPUT_FORMATS:
        return _return_dataframe(df, self.output_format)
      return output_matches(df, self.output_format, self.output_name)

    if self.query:
      if self.max_indels > 0:
        linear_df = self.indel_search()
      elif self.best_match and self.k_specified:
        results = self._search(self.k, self.max_mismatches)
        linear_df = self._best_match_filter(self._to_dataframe(results))
      elif self.best_match:
        linear_df = self.best_match_search()
      else:
        k = self.k if self.k_specified else self._auto_k(self.max_mismatches)
        self._recall_warning(k, self.max_mismatches)
        results = self._search(k, self.max_mismatches)
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
    df = pl.concat(dfs, how="vertical_relaxed") if dfs else linear_df

    if self.output_format in DATAFRAME_OUTPUT_FORMATS:
      return _return_dataframe(df, self.output_format)
    output_matches(df, self.output_format, self.output_name)

  def indel_search(self):
    """Indel search, partitioned by the k each group of queries needs.

    max_indels=2 is homogeneous: two insertions or two deletions, never one of each.
    Queries too short for k=3 run against a 2-mer table on their own rather than
    dragging the whole batch down to k=2, where 2-mer seeds explode. Every group still
    runs at k <= its own pigeonhole optimum, so the match set equals a single-k run.
    Groups are merged as columns, not frames, so an all-miss group's null column cannot
    clash with another group's String column."""
    n = self.max_indels
    short = [(q, s) for q, s in self.query if len(s) // (n + 1) < 3]
    rest = [(q, s) for q, s in self.query if len(s) // (n + 1) >= 3]

    groups = {}
    if short:
      groups.setdefault(self._clamp_k(2), []).extend(short)
    if rest:
      rest_optimal = max(2, min(len(s) for _, s in rest) // (n + 1))
      groups.setdefault(self._clamp_k(rest_optimal), []).extend(rest)

    combined = [[] for _ in range(8)]
    for k, queries in sorted(groups.items()):
      self._recall_warning(k, n, queries)
      pepidx_path = self._pepidx_path(k)
      if not os.path.isfile(pepidx_path):
        print(f"No preprocessed file found for k={k}, building index now "
              f"(this may take a moment)...")
        rs_preprocess(self.proteome_file, k, pepidx_path)
      print(f"Searching {len(queries)} peptides against {self.proteome_name} "
            f"(k={k}, max_indels={n})...")
      for i, column in enumerate(rs_indel_match(pepidx_path, queries, n)):
        combined[i].extend(column)
    return self._to_dataframe(combined, is_indels=True)

  def _clamp_k(self, optimal):
    """Resolve k for one query group: an explicit k is honored up to `optimal`, above
    which too few disjoint seeds remain to guarantee recall, so it is clamped down."""
    if self.k_specified and self.k > optimal:
      print(f"Requested k={self.k} exceeds k={optimal}, the largest k that "
            f"guarantees complete recall for these queries (pigeonhole); "
            f"using k={optimal}.")
      return optimal
    return self.k if self.k_specified else optimal

  def _recall_warning(self, k, edits, queries=None):
    """Warn where k cannot guarantee complete recall. The seed guarantee needs
    len >= k * (edits + 1), so a shorter peptide may have a real match that every seed
    misses. Shared by the mismatch and indel paths, the latter passing its own group."""
    queries = self.query if queries is None else queries
    unguaranteed = sorted({len(s) for _, s in queries if len(s) < k * (edits + 1)})
    if unguaranteed:
      print(f"k={k}, edits={edits}: complete recall is not guaranteed for query "
            f"peptides shorter than {k * (edits + 1)} residues (lengths: {unguaranteed}); "
            f"a match may exist but be missed.")

  def _auto_k(self, edits):
    """Largest k that still guarantees recall (pigeonhole): floor(L / (edits + 1)) taken
    on the shortest query peptide, so every peptide splits into edits+1 disjoint seeds.
    Floored at 2, the preprocessor minimum."""
    min_len = min(len(seq) for _, seq in self.query)
    return max(2, min_len // (edits + 1))

  def _pepidx_path(self, k):
    return os.path.join(
      self.preprocessed_files_path, f'{self.proteome_name}_{k}mers.pepidx'
    )

  def _search(self, k, max_mismatches, peptides=None):
    """Search `peptides` (default: the whole query) at k, building the index if absent."""
    pepidx_path = self._pepidx_path(k)
    if not os.path.isfile(pepidx_path):
      print(f"Preprocessing {self.proteome_name} with k={k}...")
      rs_preprocess(self.proteome_file, k, pepidx_path)
    query = peptides or self.query
    print(f"Searching {len(query)} peptides against {self.proteome_name} (k={k}, max_mismatches={max_mismatches})...")
    return rs_match(pepidx_path, query, k, max_mismatches)

  def _search_counts(self, k, max_mismatches):
    """Count hits per (peptide, mismatch level) at k, building the index if absent."""
    pepidx_path = self._pepidx_path(k)
    if not os.path.isfile(pepidx_path):
      print(f"Preprocessing {self.proteome_name} with k={k}...")
      rs_preprocess(self.proteome_file, k, pepidx_path)
    print(f"Counting {len(self.query)} peptides against {self.proteome_name} (k={k}, max_mismatches={max_mismatches})...")
    return rs_match_counts(pepidx_path, self.query, k, max_mismatches)

  def _counts_to_dataframe(self, cols):
    """Aggregate counts: one row per (query peptide, mismatch level) with a hit.
    Memory is O(unique queries), independent of total hit count."""
    qid, qseq, mm, count = cols
    if not qid:
      return pl.DataFrame(schema={
        'Query ID': pl.Utf8, 'Query Sequence': pl.Utf8,
        'Mismatches': pl.Int64, 'Count': pl.UInt64,
      })
    return pl.DataFrame({
      'Query ID': qid,
      'Query Sequence': qseq,
      'Mismatches': pl.Series(mm, dtype=pl.Int64),
      'Count': pl.Series(count, dtype=pl.UInt64),
    })

  def best_match_search(self):
    """Closest match per peptide, walking k down and the mismatch budget up until
    nothing is left unmatched."""
    peptides_remaining = self.query.copy()
    acc = tuple([] for _ in range(8))   # 8 columnar accumulators

    def collect_matched(cols):
      matched_ids = set()
      matched_col = cols[2]
      for i in range(len(cols[0])):
        if matched_col[i] is not None:
          for j in range(8):
            acc[j].append(cols[j][i])
          matched_ids.add(cols[0][i])
      return matched_ids

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
      matched_ids = collect_matched(self._search(k, max_mm, peptides_remaining))
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
        matched_ids = collect_matched(self._search(2, max_mm, peptides_remaining))
        peptides_remaining = [
          (qid, seq) for qid, seq in peptides_remaining if qid not in matched_ids
        ]
        print(f"  -> k=2, mismatches<={max_mm}: {len(peptides_remaining)} remaining")

    for qid, seq in peptides_remaining:   # leftover unmatched
      for j, val in enumerate((qid, seq, None, None, None, "[]", None, None)):
        acc[j].append(val)

    df = self._to_dataframe(acc)
    return self._best_match_filter(df)

  def _best_match_filter(self, df):
    """Keep one row per query, breaking ties in order: fewest mismatches, gene priority,
    canonical (unhyphenated) Protein ID, SwissProt reviewed, protein existence level,
    non-fragment. Unmatched rows pass through untouched."""
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

    return pl.concat([matched_df, unmatched_df], how="vertical_relaxed")

  def _metadata_table(self) -> pl.DataFrame:
    """Per-protein metadata (built once from this proteome's index) for the edge join."""
    import glob
    pattern = os.path.join(self.preprocessed_files_path, f'{self.proteome_name}_*mers.pepidx')
    idx_files = glob.glob(pattern)
    if not idx_files:
      raise FileNotFoundError(f'No .pepidx for {self.proteome_name} in {self.preprocessed_files_path}')
    pnum, pid, name, species, taxon, gene, exist, seqver, geneprio, swiss = rs_metadata(idx_files[0])
    m = pl.DataFrame({
      'protein_num': pl.Series(pnum, dtype=pl.UInt32),
      'Protein ID': pid, 'Protein Name': name, 'Species': species, 'Taxon ID': taxon,
      'Gene': gene, 'Protein Existence Level': exist, 'Sequence Version': seqver,
      'Gene Priority': geneprio, 'SwissProt Reviewed': swiss,
    })
    str_cols = ['Protein ID','Protein Name','Species','Taxon ID','Gene',
                'Protein Existence Level','Sequence Version','Gene Priority','SwissProt Reviewed']
    m = m.with_columns(
      pl.when(pl.col(c) == '').then(None).otherwise(pl.col(c)).alias(c) for c in str_cols
    )
    return m.with_columns([
      pl.col('Protein Existence Level').cast(pl.Int64),
      pl.col('Sequence Version').cast(pl.Int64),
      pl.col('Gene Priority').cast(pl.Int64),
      pl.col('SwissProt Reviewed').cast(pl.Int64).cast(pl.Boolean),
    ])

  def _final_columns(self, is_indels):
    """Column order of the results. Indel search reports Indels / Indel Positions and
    every other mode Mismatches / Mutated Positions, never both -- an empty twin column
    would imply a search that never ran."""
    edit_col = 'Indels' if is_indels else 'Mismatches'
    pos_col = 'Indel Positions' if is_indels else 'Mutated Positions'
    return [
      'Query ID','Query Sequence','Matched Sequence','Protein ID','Protein Name','Species',
      'Taxon ID','Gene', edit_col, pos_col,'Index start','Index end',
      'Protein Existence Level','Gene Priority','SwissProt Reviewed',
    ]

  def _to_dataframe(self, cols, is_indels=False):
    """Build the results frame from the columnar Rust output, joining protein metadata
    once instead of cloning it into every hit row. rs_indel_match reuses the mismatch
    slot for its edit count, so only the column name differs by mode; indel positions
    are derived here from (query, matched) rather than in Rust, and miss rows stay null."""
    qid, qseq, matched, pnum, mm, mutated, istart, iend = cols

    edit_col = 'Indels' if is_indels else 'Mismatches'
    pos_col = 'Indel Positions' if is_indels else 'Mutated Positions'
    final_columns = self._final_columns(is_indels)

    if not qid:
      schema = {c: pl.Utf8 for c in final_columns}
      for c in (edit_col, 'Index start', 'Index end', 'Protein Existence Level', 'Gene Priority'):
        schema[c] = pl.Int64
      schema['SwissProt Reviewed'] = pl.Boolean
      return pl.DataFrame(schema=schema)

    if is_indels:
      positions = [format_indel_positions(q, m) if m is not None else None
                   for q, m in zip(qseq, matched)]
    else:
      positions = mutated

    base = pl.DataFrame({
      'Query ID': qid,
      'Query Sequence': qseq,
      'Matched Sequence': matched,
      'protein_num': pl.Series(pnum, dtype=pl.UInt32),
      edit_col: pl.Series(mm, dtype=pl.Int64),
      pos_col: positions,
      'Index start': pl.Series(istart, dtype=pl.Int64),
      'Index end': pl.Series(iend, dtype=pl.Int64),
    })

    df = base.join(self._metadata_table(), on='protein_num', how='left').drop('protein_num')

    if self.sequence_version:
      df = df.with_columns(
        pl.when(pl.col('Sequence Version').is_not_null())
          .then(pl.col('Protein ID') + '.' + pl.col('Sequence Version').cast(pl.Utf8))
          .otherwise(pl.col('Protein ID'))
          .alias('Protein ID')
      )

    return df.drop('Sequence Version').select(final_columns)
