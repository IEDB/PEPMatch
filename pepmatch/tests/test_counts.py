import pytest
import polars as pl
import polars.testing as plt
from pathlib import Path
from pepmatch import Matcher


@pytest.fixture
def proteome_path() -> Path:
  return Path(__file__).parent / 'data' / 'proteome.fasta'

@pytest.fixture
def mismatch_query() -> Path:
  return Path(__file__).parent / 'data' / 'mismatching_query.fasta'

@pytest.fixture
def exact_query() -> Path:
  return Path(__file__).parent / 'data' / 'exact_match_query.fasta'


def _from_full(df: pl.DataFrame) -> pl.DataFrame:
  return (
    df.filter(pl.col('Matched Sequence').is_not_null())
      .group_by(['Query Sequence', 'Mismatches'])
      .agg(pl.len().alias('Count'))
      .with_columns(pl.col('Count').cast(pl.Int64))
      .sort(['Query Sequence', 'Mismatches'])
  )

def _from_counts(df: pl.DataFrame) -> pl.DataFrame:
  return (
    df.group_by(['Query Sequence', 'Mismatches'])
      .agg(pl.col('Count').sum().alias('Count'))
      .with_columns(pl.col('Count').cast(pl.Int64))
      .sort(['Query Sequence', 'Mismatches'])
  )


def test_counts_parity_mismatch(proteome_path, mismatch_query):
  """counts_only must equal the full output grouped by (peptide, mismatch level)."""
  full = Matcher(query=mismatch_query, proteome_file=proteome_path,
                 max_mismatches=3, k=3).match()
  counts = Matcher(query=mismatch_query, proteome_file=proteome_path,
                   max_mismatches=3, k=3, counts_only=True).match()
  assert counts.columns == ['Query ID', 'Query Sequence', 'Mismatches', 'Count']
  plt.assert_frame_equal(_from_full(full), _from_counts(counts))


def test_counts_parity_exact(proteome_path, exact_query):
  full = Matcher(query=exact_query, proteome_file=proteome_path,
                 max_mismatches=0, k=5).match()
  counts = Matcher(query=exact_query, proteome_file=proteome_path,
                   max_mismatches=0, k=5, counts_only=True).match()
  plt.assert_frame_equal(_from_full(full), _from_counts(counts))


def test_counts_only_rejects_best_match(proteome_path, exact_query):
  with pytest.raises(ValueError):
    Matcher(query=exact_query, proteome_file=proteome_path,
            counts_only=True, best_match=True)
