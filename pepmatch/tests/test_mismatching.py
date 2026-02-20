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
  return Path(__file__).parent / 'data' / 'mismatching_query.fasta'

@pytest.fixture
def expected_path() -> Path:
  return Path(__file__).parent / 'data' / 'mismatching_expected.csv'

def test_mismatching(proteome_path, query_path, expected_path):
  df = Matcher(
    query=query_path,
    proteome_file=proteome_path,
    max_mismatches=3,
    k=3,
    output_format='dataframe'
  ).match()

  expected_df = pl.read_csv(expected_path)
  cols_to_compare = ['Query Sequence', 'Matched Sequence', 'Protein ID']
  sort_key = ['Query Sequence', 'Matched Sequence']
  df_sorted = df.sort(sort_key)
  expected_df_sorted = expected_df.sort(sort_key)

  plt.assert_frame_equal(
    df_sorted.select(cols_to_compare),
    expected_df_sorted.select(cols_to_compare)
  )
