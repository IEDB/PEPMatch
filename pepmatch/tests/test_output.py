#!/usr/bin/env python3

import pandas as pd
import pytest
from pathlib import Path
from os import remove
from pepmatch import Matcher
from pepmatch.matcher import VALID_OUTPUT_FORMATS


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
def simple_dataframe():
  return pd.DataFrame(data={"col1" : [1, 2], "col2": [3, 4]})


def _creation_calls(_matcher, _path, _df):
  if _path.exists():
    remove(_path.__str__())
  _matcher.output_name = _path
  _matcher._output_matches(_df)


def test_local_creation(match, simple_dataframe):
  """Test creation given the local path"""

  local_location = Path("local_result.csv")
  _creation_calls(match, local_location, simple_dataframe)

  assert local_location.exists()
  remove(local_location.__str__())


def test_relative_creation(match, simple_dataframe):
  """Tests for creation given a relative path"""
  
  base_path = Path(__file__).parent
  relative_location = base_path / "relative.csv"
  _creation_calls(match, relative_location, simple_dataframe)

  assert relative_location.exists()
  remove(relative_location.__str__())


def test_absolute_creation(match, simple_dataframe):
  """Tests creation given an absolute path"""

  absolute_location = Path("absolute.csv").absolute()
  _creation_calls(match, absolute_location, simple_dataframe)

  assert absolute_location.exists()
  remove(absolute_location.__str__())


def test_type_creation(match, simple_dataframe):
  """Tests Matcher's ability to convert all supported types
  into the expected file's output and format
  """
  for accepted_type in VALID_OUTPUT_FORMATS[1:]:
      location = Path(f"test_output.{accepted_type}")
      if location.exists():
          remove(location.__str__())
      match.output_format = accepted_type
      _creation_calls(match, location, simple_dataframe)
      assert location.exists()
      remove(location.__str__())
