import polars as pl
import pytest
from pathlib import Path
from os import remove
from pepmatch import Matcher
from pepmatch.matcher import VALID_OUTPUT_FORMATS
from pepmatch.helpers import output_matches 

@pytest.fixture
def proteome_path() -> Path:
  return Path(__file__).parent / 'data' / 'proteome.fasta'

@pytest.fixture
def query_path() -> Path:
  return Path(__file__).parent / 'data' / 'exact_match_query.fasta'

@pytest.fixture
def match(proteome_path, query_path):
  return Matcher(
    query=query_path,
    proteome_file=proteome_path,
    max_mismatches=0,
    k=0,
    output_format='csv'
  )

@pytest.fixture
def simple_dataframe() -> pl.DataFrame:
  return pl.DataFrame({"col1": [1, 2], "col2": [3, 4]})

def _creation_calls(path: Path, df: pl.DataFrame):
  """Helper function to test file creation."""
  
  output_name_no_ext = str(path.with_suffix(''))
  output_format = path.suffix[1:]

  if path.exists():
    remove(str(path))
  
  output_matches(df, output_format, output_name_no_ext)

def test_local_creation(simple_dataframe):
  """Test creation given the local path."""
  local_location = Path("local_result.csv")
  _creation_calls(local_location, simple_dataframe)
  assert local_location.exists()
  remove(str(local_location))

def test_relative_creation(simple_dataframe):
  """Tests for creation given a relative path."""
  base_path = Path(__file__).parent
  relative_location = base_path / "relative.csv"
  _creation_calls(relative_location, simple_dataframe)
  assert relative_location.exists()
  remove(str(relative_location))

def test_absolute_creation(simple_dataframe):
  """Tests creation given an absolute path."""
  absolute_location = Path("absolute.csv").absolute()
  _creation_calls(absolute_location, simple_dataframe)
  assert absolute_location.exists()
  remove(str(absolute_location))

def test_type_creation(simple_dataframe):
  """Tests Matcher's ability to convert all supported types
  into the expected file's output and format.
  """
  # skip 'dataframe' as it doesn't create a file
  for accepted_type in VALID_OUTPUT_FORMATS[1:]:
    location = Path(f"test_output.{accepted_type}")
    if location.exists():
      remove(str(location))
    
    _creation_calls(location, simple_dataframe)
    
    assert location.exists()
    remove(str(location))
