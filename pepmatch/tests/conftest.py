import pytest
import multiprocessing as mp

@pytest.fixture(scope="session", autouse=True)
def set_multiprocessing_start_method():
  """Set the multiprocessing start method to 'spawn' for all tests.

  This is necessary to avoid deadlocks when forking from a multi-threaded
  pytest process. 'spawn' is safer and more cross-platform compatible."""

  try:
    mp.set_start_method("spawn")
  except RuntimeError:
    if mp.get_start_method() != "spawn":
      print(
        f"Warning: multiprocessing start method is '{mp.get_start_method()}'"
        " and could not be set to 'spawn'. This may cause issues."
      )
    pass
