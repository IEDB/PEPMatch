import pytest
import polars as pl
import polars.testing as plt
from pathlib import Path
from pepmatch import Matcher

@pytest.fixture
def proteome_path() -> Path:
  return Path(__file__).parent / 'data' / 'proteome.fasta'

@pytest.fixture
def query_path() -> Path:
  return Path(__file__).parent / 'data' / 'exact_match_query.fasta'

@pytest.fixture
def expected_path() -> Path:
  return Path(__file__).parent / 'data' / 'exact_match_expected.csv'

def test_exact_match_k5(proteome_path, query_path, expected_path):
  df = Matcher(
    query=query_path,
    proteome_file=proteome_path,
    max_mismatches=0,
    k=5,
    output_format='dataframe'
  ).match()

  df = df.sort('Query Sequence')
  expected_df = pl.read_csv(expected_path).sort('Query Sequence')
  plt.assert_series_equal(df['Protein ID'], expected_df['Protein ID'])

def test_exact_match_k9(proteome_path, query_path, expected_path):
  df = Matcher(
    query=query_path,
    proteome_file=proteome_path,
    max_mismatches=0,
    k=9,
    output_format='dataframe'
  ).match()

  df = df.sort('Query Sequence')
  expected_df = pl.read_csv(expected_path).sort('Query Sequence')
  plt.assert_series_equal(df['Protein ID'], expected_df['Protein ID'])
