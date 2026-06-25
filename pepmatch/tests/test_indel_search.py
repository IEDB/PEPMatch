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


def test_terminal_deletion_blocked(proteome_path):
  df = Matcher(
    query=['QNALVEATRFC'],
    proteome_file=proteome_path,
    max_indels=1,
    output_format='dataframe'
  ).match()
  assert df['Matched Sequence'].is_null().all(), (
    'Terminal deletion should be blocked but got match(es): '
    + str(df.filter(pl.col('Matched Sequence').is_not_null())
            .select(['Query Sequence', 'Matched Sequence']))
  )
