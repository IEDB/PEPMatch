#!/usr/bin/env python3
"""
Containing common fixtures used by most test files
"""

import pytest
from pathlib import Path


@pytest.fixture
def proteome_path() -> Path:
  return Path(__file__).parent / 'data' / 'proteome.fasta'


@pytest.fixture
def query_path() -> Path:
  return Path(__file__).parent / 'data' / 'mismatching_query.fasta'


@pytest.fixture
def expected_path() -> Path:
  return Path(__file__).parent / 'data' / 'mismatching_expected.csv'