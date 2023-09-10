#!/usr/bin/env python3

import os
import pytest
import pandas as pd
import pandas.testing as pdt
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
  df = df.sort_values(by=['Query Sequence']).reset_index(drop=True)

  os.remove('proteome_2mers.pickle')
  os.remove('proteome_3mers.pickle')
  os.remove('proteome_6mers.pickle')
  os.remove('proteome_metadata.pickle')
  os.remove('proteome.db')

  expected_df = pd.read_csv(expected_path)
  expected_df= expected_df.sort_values(by=['Query Sequence']).reset_index(drop=True)
  pdt.assert_series_equal(df['Protein ID'], expected_df['Protein ID'])

