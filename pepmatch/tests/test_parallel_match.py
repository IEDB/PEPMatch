#!/usr/bin/env python3

import os
import pandas as pd
import pandas.testing as pdt
from pepmatch import Preprocessor, ParallelMatcher

  
def test_mismatching(proteome_path, query_path, expected_path):
  """Test mismatching of query peptides to a proteome. The query is various peptides
  with different mismatches searched in the Dugbe virus proteome (isolate ArD44313)."""

  Preprocessor(proteome_path).pickle_proteome(k=3)

  df = ParallelMatcher(
    query=query_path,
    proteome_file=proteome_path,
    max_mismatches=3,
    k=3,
    output_format='dataframe',
    n_jobs=2
  ).match()
  df = df.sort_values(by=['Query Sequence']).reset_index(drop=True)

  os.remove('proteome_3mers.pkl')
  os.remove('proteome_metadata.pkl')

  expected_df = pd.read_csv(expected_path)
  expected_df = expected_df.sort_values(
    by=['Query Sequence', 'Matched Sequence']
  ).reset_index(drop=True)
  
  pdt.assert_frame_equal(
    df[['Query Sequence', 'Matched Sequence', 'Protein ID']],
    expected_df[['Query Sequence', 'Matched Sequence', 'Protein ID']]
  )

