#!/usr/bin/env python3

import pandas as pd
import pandas.testing as pdt
from common_features import *
from pepmatch import Matcher


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
  """Test searching discontinuous peptides in a proteome."""
  df = Matcher(
    query=query,
    proteome_file=proteome_path,
    max_mismatches=0,
    output_format='dataframe'
  ).match()
  df = df.sort_values(by=['Query Sequence']).reset_index(drop=True)

  expected_df = pd.read_csv(expected_path)
  expected_df= expected_df.sort_values(by=['Query Sequence']).reset_index(drop=True)
  pdt.assert_series_equal(df['Protein ID'], expected_df['Protein ID'])