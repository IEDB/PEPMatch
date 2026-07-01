import pytest
import polars as pl
import polars.testing as plt
from pathlib import Path
from pepmatch import Matcher

# QNALVEATRFC is designed so the only match in the test proteome would require
# deleting Q at query position 0 (a terminal deletion). The guard must block it.

@pytest.fixture
def proteome_path() -> Path:
  return Path(__file__).parent / 'data' / 'proteome.fasta'

@pytest.fixture
def query_path() -> Path:
  return Path(__file__).parent / 'data' / 'indel_query.fasta'

@pytest.fixture
def expected_path() -> Path:
  return Path(__file__).parent / 'data' / 'indel_expected.csv'

def test_indel_search(proteome_path, query_path, expected_path):
  df = Matcher(
    query=query_path,
    proteome_file=proteome_path,
    max_indels=1,
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


def test_terminal_deletion_allowed_with_single_residue_buffer(tmp_path):
  # Regression test: p_idx == 0 and p_idx == protein_len - 1 are real,
  # existing residues, not out-of-bounds — a deletion's gap only needs 1
  # residue of protein context on either side to be a genuine indel, not 2.
  proteome_path = tmp_path / 'proteome.fasta'
  proteome_path.write_text('>Front\nXBCDE\n>Back\nABCDX\n')
  df = Matcher(
    query=['ABCDE'],
    proteome_file=str(proteome_path),
    max_indels=1,
    preprocessed_files_path=str(tmp_path),
    output_format='dataframe'
  ).match()
  hits = set(zip(df['Protein ID'].to_list(), df['Matched Sequence'].to_list()))
  assert ('Front.1', 'BCDE') in hits
  assert ('Back.1', 'ABCD') in hits


def test_query_terminal_deletion_found(proteome_path):
  # Deletion of Q at query position 0 hits NALVEATRFC at protein position 3
  # in Q8V336 — not a protein boundary, so the match is valid.
  df = Matcher(
    query=['QNALVEATRFC'],
    proteome_file=proteome_path,
    max_indels=1,
    output_format='dataframe'
  ).match()
  matches = df.filter(pl.col('Matched Sequence').is_not_null())
  assert 'NALVEATRFC' in matches['Matched Sequence'].to_list(), (
    'Expected query-terminal deletion match NALVEATRFC not found'
  )
