#!/usr/bin/env python3

import os


from pepmatch import Matcher


def test_exact_match():
  # paths
  test_script_dir = os.path.dirname(os.path.realpath(__file__))
  proteome_fasta = os.path.join(
    test_script_dir, '../benchmarking/proteomes/human.fasta')


  # match MHC ligands (9-mers) to human proteome
  df = Matcher(
    query=[
      'R377, Q408, Q432, H433, F436, V441, S442, S464, K467, K489, I491, S492, N497'],
    proteome_file=proteome_fasta,
    max_mismatches=0,
    output_format='dataframe').match()

  assert df['Protein ID'].iloc[0] == 'P00533.1'