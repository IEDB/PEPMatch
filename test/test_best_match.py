#!/usr/bin/env python3

import os
import glob

import pandas as pd
import pandas.testing as pdt

from pepmatch import Preprocessor, Matcher
  
def test_mismatch():

  # paths
  test_script_dir = os.path.dirname(os.path.realpath(__file__))
  proteome_fasta = os.path.join(test_script_dir, '../benchmarking/proteomes/human.fasta')
  query_fasta = os.path.join(test_script_dir, '../benchmarking/queries/milk_peptides_test.fasta')
  expected_csv = os.path.join(test_script_dir, '../benchmarking/expected/milk_peptides_expected.csv')

  # match neoepitopes to human proteome
  df = Matcher(
    query=query_fasta,
    proteome_file=proteome_fasta,
    best_match=True,
    output_format='dataframe').match()

  # remove preprocessed files
  os.remove('human.db')
  os.remove('human_7mers.pickle')
  os.remove('human_3mers.pickle')
  os.remove('human_2mers.pickle')
  os.remove('human_metadata.pickle')

  # load the expected data
  expected_df = pd.read_csv(expected_csv)

  # select only the necessary columns to test for
  df = df[['Query Sequence', 'Matched Sequence', 'Protein ID', 'Index start']]
  expected_df = expected_df[['Query Sequence', 'Matched Sequence', 'Protein ID', 'Index start']]

  # convert Index start to int64 for comparison
  df['Index start'] = df['Index start'].astype('int64')
  expected_df['Index start'] = expected_df['Index start'].astype('int64')

  # sort dataframes by Query Sequence, Protein ID, and Index start
  df = df.sort_values(by=['Query Sequence', 'Protein ID', 'Index start']).reset_index(drop=True)
  expected_df = expected_df.sort_values(by=['Query Sequence', 'Protein ID', 'Index start']).reset_index(drop=True)

  # assert dataframes are equal
  pdt.assert_frame_equal(df, expected_df, check_dtype=False, check_exact=True)
