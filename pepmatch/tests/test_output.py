import polars as pl
import pytest
from pathlib import Path
from os import remove
from pepmatch import Matcher
from pepmatch.matcher import VALID_OUTPUT_FORMATS, output_matches

@pytest.fixture
def proteome_path() -> Path:
  return Path(__file__).parent / 'data' / 'proteome.fasta'

@pytest.fixture
def query_path() -> Path:
  return Path(__file__).parent / 'data' / 'exact_match_query.fasta'

@pytest.fixture
def simple_dataframe() -> pl.DataFrame:
  return pl.DataFrame({"col1": [1, 2], "col2": [3, 4]})

def _creation_calls(path: Path, df: pl.DataFrame):
  output_name_no_ext = str(path.with_suffix(''))
  output_format = path.suffix[1:]
  if path.exists():
    remove(str(path))
  output_matches(df, output_format, output_name_no_ext)

def test_local_creation(simple_dataframe):
  local_location = Path("local_result.csv")
  _creation_calls(local_location, simple_dataframe)
  assert local_location.exists()
  remove(str(local_location))

def test_relative_creation(simple_dataframe):
  base_path = Path(__file__).parent
  relative_location = base_path / "relative.csv"
  _creation_calls(relative_location, simple_dataframe)
  assert relative_location.exists()
  remove(str(relative_location))

def test_absolute_creation(simple_dataframe):
  absolute_location = Path("absolute.csv").absolute()
  _creation_calls(absolute_location, simple_dataframe)
  assert absolute_location.exists()
  remove(str(absolute_location))

def test_type_creation(simple_dataframe):
  for accepted_type in VALID_OUTPUT_FORMATS[1:]:
    location = Path(f"test_output.{accepted_type}")
    if location.exists():
      remove(str(location))
    _creation_calls(location, simple_dataframe)
    assert location.exists()
    remove(str(location))
