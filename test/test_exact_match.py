#!/usr/bin/env python3

import unittest
import os

import pandas as pd
import pandas.testing as pdt

from pepmatch import Preprocessor, Matcher

class TestExactMatch(unittest.TestCase):
  def test_exact_match(self):
    # preprocess human proteome
    Preprocessor('../benchmarking/proteomes/human.fasta').sql_proteome(9)

    # match MHC ligands (9-mers) to human proteome
    df = Matcher(
      query='../benchmarking/queries/mhc_ligands_test.fasta',
      proteome_file='../benchmarking/proteomes/human.fasta',
      max_mismatches=0,
      k=9,
      output_format='dataframe').match()

    # remove preprocessed file
    os.remove('human.db')

    # load the expected data
    expected_df = pd.read_csv('../benchmarking/expected/mhc_ligands_expected.csv')

    # select only the necessary columns to test for
    df = df[['Query Sequence', 'Matched Sequence', 'Protein ID', 'Index start']]
    expected_df = expected_df[['Query Sequence', 'Matched Sequence', 'Protein ID', 'Index start']]

    # sort dataframes by Query Sequence and Protein ID
    df = df.sort_values(by=['Query Sequence', 'Protein ID']).reset_index(drop=True)
    expected_df = expected_df.sort_values(by=['Query Sequence', 'Protein ID']).reset_index(drop=True)

    # assert dataframes are equal
    pdt.assert_frame_equal(df, expected_df, check_dtype=False, check_exact=True)

if __name__ == '__main__':
  unittest.main()