import pytest
import polars as pl
import polars.testing as plt
from pathlib import Path
from pepmatch import Matcher

@pytest.fixture
def proteome_path() -> Path:
  return Path(__file__).parent / 'data' / 'proteome.fasta'

@pytest.fixture
def expected_path() -> Path:
  return Path(__file__).parent / 'data' / 'discontinuous_expected.csv'

@pytest.fixture
def query() -> list:
  return [
    'L354, V420, G461, Q468, E486, K499, D501, M503, G509',
    'T819, L822, A835, F840, S871, Y884, P886',
    'S2760, V2763, E2773, D2805, T2819, S2831, E2844, R2852, L2863'
  ]

def test_discontinuous_search(proteome_path, query, expected_path):
  df = Matcher(
    query=query,
    proteome_file=proteome_path,
    max_mismatches=0,
    output_format='dataframe'
  ).match()

  df = df.sort('Query Sequence')
  expected_df = pl.read_csv(expected_path).sort('Query Sequence')
  plt.assert_series_equal(df['Protein ID'], expected_df['Protein ID'])
