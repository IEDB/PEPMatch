import sys

import polars as pl
import polars.testing as plt
import pytest

from pathlib import Path

from pepmatch import Matcher

# A row-for-row record of indel search output: 53 queries (exact windows, one and two edits of
# each kind, windows inside homopolymer and short repeats, at protein ends, a motif that occurs
# twice, a mixed edit that is never found, short peptides, decoys) against twelve proteins built
# around repeat tracts, two of them isoform-like (a stretch shared at the same position, so hits
# on both have the same coordinates and only the protein tells them apart). The brute-force
# property tests compare sets of matches, so they cannot see a change that reorders, duplicates
# or relabels rows without changing the set; this does. The last six queries were picked because
# an engine that pairs budget splits in another order, or dedups hits without the protein, fails
# on them.
# Recorded from the engine as of v1.17.2. Regenerate only for a deliberate output change:
#   python pepmatch/tests/test_indel_golden.py --write
DATA = Path(__file__).parent / 'data'
EXPECTED = DATA / 'indel_golden_expected.csv'
RUNS = {
  '1 indel': dict(max_indels=1),
  '2 indels': dict(max_indels=2),
  '1 indel, k=2': dict(max_indels=1, k=2),
  '2 indels, k=2': dict(max_indels=2, k=2),
}


def search(run, preprocessed_files_path):
  df = Matcher(
    query=DATA / 'indel_golden_query.fasta',
    proteome_file=DATA / 'indel_golden_proteome.fasta',
    preprocessed_files_path=str(preprocessed_files_path),
    output_format='dataframe',
    **RUNS[run],
  ).match()
  return df.with_columns(pl.all().cast(pl.Utf8))   # compared as text, so the CSV round-trips exactly


@pytest.mark.parametrize('run', RUNS)
def test_indel_search_output_is_unchanged_row_for_row(tmp_path, run):
  expected = pl.read_csv(EXPECTED, infer_schema_length=0).filter(pl.col('Run') == run).drop('Run')
  got = search(run, tmp_path)
  assert got['Matched Sequence'].is_not_null().sum() > 0
  plt.assert_frame_equal(got, expected)


if __name__ == '__main__' and sys.argv[1:] == ['--write']:
  import tempfile
  with tempfile.TemporaryDirectory() as tmp:
    frames = [search(run, tmp).insert_column(0, pl.lit(run).alias('Run')) for run in RUNS]
  pl.concat(frames).write_csv(EXPECTED)
  print(f'wrote {sum(f.height for f in frames)} rows to {EXPECTED}')
