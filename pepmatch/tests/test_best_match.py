import os
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
  return Path(__file__).parent / 'data' / 'best_match_query.fasta'

@pytest.fixture
def expected_path() -> Path:
  return Path(__file__).parent / 'data' / 'best_match_expected.csv'

def test_best_match(proteome_path, query_path, expected_path):
  """Test best match of query peptides to a proteome. The query is various peptides
  with different best matches searched in the Dugbe virus proteome (isolate ArD44313).
  The peptides are taken from the human proteome at random."""

  df = Matcher(
    query=query_path,
    proteome_file=proteome_path,
    best_match=True,
    output_format='dataframe'
  ).match()

  for file_path in ['proteome_2mers.pkl', 'proteome_3mers.pkl', 'proteome_6mers.pkl', 'proteome_metadata.pkl', 'proteome.db']:
      if os.path.exists(file_path):
          os.remove(file_path)

  expected_df = pl.read_csv(expected_path)

  cols_to_compare = ['Query Sequence', 'Matched Sequence', 'Protein ID']
  sort_key = ['Query Sequence', 'Matched Sequence']

  df_sorted = df.sort(sort_key)
  expected_df_sorted = expected_df.sort(sort_key)

  plt.assert_frame_equal(
    df_sorted.select(cols_to_compare),
    expected_df_sorted.select(cols_to_compare),
  )
