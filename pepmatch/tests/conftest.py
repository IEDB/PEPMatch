import pytest
import glob
import os

@pytest.fixture(autouse=True, scope="session")
def cleanup_pepidx():
  yield
  for f in glob.glob("*.pepidx"):
    os.remove(f)
  for f in glob.glob("PEPMatch_results.*"):
    os.remove(f)
