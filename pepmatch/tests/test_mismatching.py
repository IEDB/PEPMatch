#!/usr/bin/env python3

import os
import pytest
import pandas as pd
import pandas.testing as pdt
from pathlib import Path

from pepmatch import Preprocessor, Matcher


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
  """Test mismatching of query peptides to a proteome. The query is various peptides
  with different mismatches searched in the Dugbe virus proteome (isolate ArD44313)."""

  Preprocessor(proteome_path).pickle_proteome(k=3)

  df = Matcher(
    query=query_path,
    proteome_file=proteome_path,
    max_mismatches=3,
    k=3,
    output_format='dataframe'
  ).match()
  df = df.sort_values(by=['Query Sequence']).reset_index(drop=True)

  os.remove('proteome_3mers.pickle')
  os.remove('proteome_metadata.pickle')

  expected_df = pd.read_csv(expected_path)
  expected_df= expected_df.sort_values(by=['Query Sequence']).reset_index(drop=True)
  pdt.assert_series_equal(df['Protein ID'], expected_df['Protein ID'])

