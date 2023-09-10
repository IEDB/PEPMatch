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
  return Path(__file__).parent / 'data' / 'exact_match_query.fasta'


@pytest.fixture
def expected_path() -> Path:
  return Path(__file__).parent / 'data' / 'exact_match_expected.csv'


def test_exact_match(proteome_path, query_path, expected_path):
  """Test exact matching of query peptides to a proteome. The query is various peptides
  searched in the Dugbe virus proteome (isolate ArD44313). Test for k=4 and k=9."""

  preprocessor = Preprocessor(proteome_path)
  preprocessor.sql_proteome(k=4)
  preprocessor.sql_proteome(k=9)

  # match using k=4
  df = Matcher(
    query=query_path,
    proteome_file=proteome_path,
    max_mismatches=0,
    k=4,
    output_format='dataframe'
  ).match()
  df = df.sort_values(by=['Query Sequence']).reset_index(drop=True)

  expected_df = pd.read_csv(expected_path)
  expected_df= expected_df.sort_values(by=['Query Sequence']).reset_index(drop=True)
  pdt.assert_series_equal(df['Protein ID'], expected_df['Protein ID'])

  # match using k=9
  df = Matcher(
    query=query_path,
    proteome_file=proteome_path,
    max_mismatches=0,
    k=9,
    output_format='dataframe'
  ).match()
  df = df.sort_values(by=['Query Sequence']).reset_index(drop=True)

  os.remove('proteome.db')

  pdt.assert_series_equal(df['Protein ID'], expected_df['Protein ID'])

