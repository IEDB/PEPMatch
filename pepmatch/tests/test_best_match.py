import pytest
import polars as pl
from pathlib import Path
from pepmatch import Matcher

@pytest.fixture
def proteome_path() -> Path:
  return Path(__file__).parent / 'data' / 'proteome.fasta'

@pytest.fixture
def query_path() -> Path:
  return Path(__file__).parent / 'data' / 'best_match_query.fasta'

def test_best_match(proteome_path, query_path):
  df = Matcher(
    query=query_path,
    proteome_file=proteome_path,
    best_match=True,
    output_format='dataframe'
  ).match()

  assert df.height == 5
  assert df.filter(pl.col("Matched Sequence").is_not_null()).height >= 3
  assert df.filter(pl.col("Mismatches").is_not_null()).select("Mismatches").to_series().min() >= 0
