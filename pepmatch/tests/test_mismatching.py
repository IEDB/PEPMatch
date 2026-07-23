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
  return Path(__file__).parent / 'data' / 'mismatching_query.fasta'

@pytest.fixture
def expected_path() -> Path:
  return Path(__file__).parent / 'data' / 'mismatching_expected.csv'

def test_mismatching(proteome_path, query_path, expected_path):
  df = Matcher(
    query=query_path,
    proteome_file=proteome_path,
    max_mismatches=3,
    k=3,
    output_format='dataframe'
  ).match()

  expected_df = pl.read_csv(expected_path)
  cols_to_compare = ['Query Sequence', 'Matched Sequence', 'Protein ID']
  sort_key = ['Query Sequence', 'Matched Sequence']
  df_sorted = df.sort(sort_key)
  expected_df_sorted = expected_df.sort(sort_key)

  plt.assert_frame_equal(
    df_sorted.select(cols_to_compare),
    expected_df_sorted.select(cols_to_compare)
  )


def test_mismatch_recall_warning_fires_for_short_query(tmp_path, capsys):
  # Mismatch shares the indel search's pigeonhole limit: recall is only guaranteed for
  # peptides with len >= k*(max_mismatches+1). A length-4 query at 2 mismatches derives
  # k=2 (floor), so the threshold is 6 and the query is sub-threshold -- the search must
  # WARN (same shared _recall_warning as the indel path) but still run.
  proteome_path = tmp_path / 'proteome.fasta'
  proteome_path.write_text('>P1\nZZABCDEFGHIKLMZZ\n')
  Matcher(query=['ABCE'], proteome_file=str(proteome_path), max_mismatches=2,
          preprocessed_files_path=str(tmp_path), output_format='dataframe').match()
  out = capsys.readouterr().out
  assert 'complete recall is not guaranteed' in out
  assert 'lengths: [4]' in out


def test_mismatch_recall_warning_silent_when_guaranteed(tmp_path, capsys):
  # A length-9 query at 2 mismatches derives k=3, threshold 9, so every query length is
  # covered -- no warning should be printed.
  proteome_path = tmp_path / 'proteome.fasta'
  proteome_path.write_text('>P1\nZZABCDEFGHIKLMNPQRSZZ\n')
  Matcher(query=['ABCDEFGHI'], proteome_file=str(proteome_path), max_mismatches=2,
          preprocessed_files_path=str(tmp_path), output_format='dataframe').match()
  out = capsys.readouterr().out
  assert 'complete recall is not guaranteed' not in out
