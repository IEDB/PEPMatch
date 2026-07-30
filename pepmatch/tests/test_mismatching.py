import pytest
import polars as pl
import polars.testing as plt
from pathlib import Path
from pepmatch import Matcher

@pytest.fixture
def proteome_path() -> Path:
  return Path(__file__).parent / 'data' / 'proteome.fasta'

@pytest.fixture
def query_path() -> Path:
  return Path(__file__).parent / 'data' / 'mismatching_query.fasta'

@pytest.fixture
def expected_path() -> Path:
  return Path(__file__).parent / 'data' / 'mismatching_expected.csv'

def test_mismatching(proteome_path, query_path, expected_path):
  df = Matcher(
    query=query_path,
    proteome_file=proteome_path,
    max_mismatches=3,
    k=3,
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


def test_mismatch_recall_warning_fires_for_short_query(tmp_path, capsys):
  # Mismatch shares the indel search's pigeonhole limit: recall is only guaranteed for
  # peptides with len >= k*(max_mismatches+1). A length-4 query at 2 mismatches derives
  # k=2 (floor), so the threshold is 6 and the query is sub-threshold -- the search must
  # WARN (same shared _recall_warning as the indel path) but still run.
  proteome_path = tmp_path / 'proteome.fasta'
  proteome_path.write_text('>P1\nZZABCDEFGHIKLMZZ\n')
  Matcher(query=['ABCE'], proteome_file=str(proteome_path), max_mismatches=2,
          preprocessed_files_path=str(tmp_path), output_format='dataframe').match()
  out = capsys.readouterr().out
  assert 'complete recall is not guaranteed' in out
  assert 'lengths: [4]' in out


def test_mismatch_recall_warning_silent_when_guaranteed(tmp_path, capsys):
  # A length-9 query at 2 mismatches derives k=3, threshold 9, so every query length is
  # covered -- no warning should be printed.
  proteome_path = tmp_path / 'proteome.fasta'
  proteome_path.write_text('>P1\nZZABCDEFGHIKLMNPQRSZZ\n')
  Matcher(query=['ABCDEFGHI'], proteome_file=str(proteome_path), max_mismatches=2,
          preprocessed_files_path=str(tmp_path), output_format='dataframe').match()
  out = capsys.readouterr().out
  assert 'complete recall is not guaranteed' not in out


# --- Hamming oracle machinery -------------------------------------------------
#
# The mismatch engine has TWO code paths: when `len(query) % k == 0` it verifies the
# flanks as whole k-mer blocks, otherwise it verifies them residue by residue. Both
# must agree with a plain Hamming scan, so the proteome/query set below deliberately
# contains queries of length 9, 10, 11 and 12 -- at k=3 and k=4 each regime is
# populated (asserted in the test, so a future query edit cannot silently drop one).
#
# Every record is padded with a run of W at both ends, and no query contains a W.
# The engine's index is a flat concatenation with no per-protein end guard (see
# test_mismatch_does_not_cross_protein_boundary), so a window may run off a record
# into the next one; with >= len(query) W's on both sides of every boundary, any
# such window is all W and therefore far outside every threshold. That keeps this
# test about the two in-protein paths rather than the boundary bug.

MISMATCH_PAD = 'W' * 12

ORACLE_CORES = {
  'P1': 'ACDEFGHIK' + 'GG' + 'ACDEFGHIL' + 'GG' + 'ACDEFGHIK',
  'P2': 'LMNPQRSTVY' + 'PP' + 'LMNPQRSSVY' + 'PP' + 'LMNAQRSSVY',
  'P3': 'FGHIKLMNPQR' + 'DD' + 'FGHIKLMNPQK' + 'DD' + 'FGHIKLMNPKK',
  'P4': 'ACDEFGHIKLMN' + 'EE' + 'ACDEFGHIKLMY' + 'EE' + 'ACDEFGHIKLNY',
}

ORACLE_PROTEINS = {pid: MISMATCH_PAD + core + MISMATCH_PAD for pid, core in ORACLE_CORES.items()}

ORACLE_QUERIES = ['ACDEFGHIK', 'LMNPQRSTVY', 'FGHIKLMNPQR', 'ACDEFGHIKLMN']


def write_fasta(path, proteins):
  """SwissProt-style records so Protein ID is the accession between the first two pipes."""
  with open(path, 'w') as f:
    for i, (pid, seq) in enumerate(proteins.items(), 1):
      f.write(
        f'>sp|{pid}|T{i}_HUMAN Test protein {i} OS=Homo sapiens OX=9606 GN=G{i} PE=1 SV=1\n'
        f'{seq}\n'
      )


def hamming(a: str, b: str) -> int:
  return sum(x != y for x, y in zip(a, b))


def brute_force_hits(proteins, queries, max_mismatches):
  """Independent oracle: every window of len(query) in every protein whose Hamming
  distance is within the threshold, as (query, protein id, window, 1-based start)."""
  hits = set()
  for pid, seq in proteins.items():
    for query in queries:
      for start in range(len(seq) - len(query) + 1):
        window = seq[start:start + len(query)]
        if hamming(query, window) <= max_mismatches:
          hits.add((query, pid, window, start + 1))
  return hits


def reported_hits(df):
  return set(zip(
    df['Query Sequence'], df['Protein ID'], df['Matched Sequence'], df['Index start']
  ))


def assert_rows_sound(df, max_mismatches):
  """Every reported row must be internally consistent, independent of any oracle."""
  matched_rows = df.filter(pl.col('Matched Sequence').is_not_null())
  assert matched_rows.height > 0
  for row in matched_rows.iter_rows(named=True):
    query, matched = row['Query Sequence'], row['Matched Sequence']
    assert len(matched) == len(query), (query, matched)
    assert row['Mismatches'] == hamming(query, matched), (query, matched)
    assert row['Mismatches'] <= max_mismatches, (query, matched)
    expected_positions = [i + 1 for i, (a, b) in enumerate(zip(query, matched)) if a != b]
    assert row['Mutated Positions'] == f'[{", ".join(str(p) for p in expected_positions)}]'
    assert len(expected_positions) == row['Mismatches']
    assert row['Index end'] == row['Index start'] + len(query) - 1


@pytest.fixture
def oracle_proteome(tmp_path) -> Path:
  path = tmp_path / 'oracle.fasta'
  write_fasta(path, ORACLE_PROTEINS)
  return path


# Each (k, max_mismatches) pair satisfies len >= k*(max_mismatches+1) for every query
# above, so the pigeonhole seeding guarantees complete recall and the engine's result
# set must equal the oracle exactly -- no missing hits, no extras.
@pytest.mark.parametrize('k,max_mismatches', [(3, 1), (3, 2), (4, 1)])
def test_mismatch_matches_hamming_oracle(oracle_proteome, tmp_path, k, max_mismatches):
  # Both engine path regimes must be exercised by this query set at this k.
  regimes = {len(q) % k == 0 for q in ORACLE_QUERIES}
  assert regimes == {True, False}, f'k={k} does not exercise both mismatch paths'
  assert all(len(q) >= k * (max_mismatches + 1) for q in ORACLE_QUERIES)

  df = Matcher(
    query=ORACLE_QUERIES,
    proteome_file=str(oracle_proteome),
    max_mismatches=max_mismatches,
    k=k,
    preprocessed_files_path=str(tmp_path),
    sequence_version=False,
    output_format='dataframe'
  ).match()

  expected = brute_force_hits(ORACLE_PROTEINS, ORACLE_QUERIES, max_mismatches)
  assert reported_hits(df) == expected

  # Guard against a vacuous oracle: the fixture must produce hits in every protein and
  # at every mismatch count from 0 up to the threshold.
  assert len(expected) >= 12
  assert set(df['Protein ID'].to_list()) == set(ORACLE_PROTEINS)
  assert set(df['Mismatches'].to_list()) == set(range(max_mismatches + 1))
  assert_rows_sound(df, max_mismatches)


@pytest.mark.parametrize('k,max_mismatches', [(3, 1), (3, 2), (4, 1), (3, 3)])
def test_mismatch_reported_rows_are_self_consistent(oracle_proteome, tmp_path, k, max_mismatches):
  # Soundness holds regardless of recall: (3, 3) is deliberately below the pigeonhole
  # threshold for the shortest query, so hits may be missing -- but whatever IS reported
  # must still line up with the query residue for residue.
  df = Matcher(
    query=ORACLE_QUERIES,
    proteome_file=str(oracle_proteome),
    max_mismatches=max_mismatches,
    k=k,
    preprocessed_files_path=str(tmp_path),
    sequence_version=False,
    output_format='dataframe'
  ).match()
  assert_rows_sound(df, max_mismatches)


def test_mismatch_rows_sound_on_reference_proteome(proteome_path, query_path, tmp_path):
  # Same soundness contract on the real reference data behind the golden CSV above:
  # Mismatches and Mutated Positions must be recomputable from the two sequences.
  df = Matcher(
    query=query_path,
    proteome_file=proteome_path,
    max_mismatches=3,
    k=3,
    preprocessed_files_path=str(tmp_path),
    output_format='dataframe'
  ).match()
  assert_rows_sound(df, 3)


def test_mismatch_includes_hits_below_the_threshold(tmp_path):
  # max_mismatches is an upper bound, not an exact count: a query with both an exact
  # occurrence and a 1-off occurrence must report BOTH when searched at 2 mismatches.
  proteome_path = tmp_path / 'proteome.fasta'
  core = 'ACDEFGHIK' + 'GG' + 'ACDEFGHIL'
  proteins = {'P1': MISMATCH_PAD + core + MISMATCH_PAD}
  write_fasta(proteome_path, proteins)

  df = Matcher(
    query=['ACDEFGHIK'],
    proteome_file=str(proteome_path),
    max_mismatches=2,
    k=3,
    preprocessed_files_path=str(tmp_path),
    sequence_version=False,
    output_format='dataframe'
  ).match()

  assert reported_hits(df) == brute_force_hits(proteins, ['ACDEFGHIK'], 2)
  by_matched = {row['Matched Sequence']: row for row in df.iter_rows(named=True)}
  assert by_matched['ACDEFGHIK']['Mismatches'] == 0
  assert by_matched['ACDEFGHIK']['Mutated Positions'] == '[]'
  assert by_matched['ACDEFGHIL']['Mismatches'] == 1
  assert by_matched['ACDEFGHIL']['Mutated Positions'] == '[9]'


@pytest.mark.parametrize('query', ['ACDACDACD', 'ACDACDACDA'])
def test_mismatch_reports_each_location_once(tmp_path, query):
  # A tandem repeat makes every seed k-mer of the query hit every repeat unit, so one
  # window is reachable from several seeds. Each (Protein ID, Index start) must still be
  # reported exactly once. Both lengths run at k=3, covering both path regimes (9 % 3 == 0,
  # 10 % 3 != 0) -- the block path and the residue path dedupe independently.
  proteome_path = tmp_path / 'proteome.fasta'
  proteins = {'P1': MISMATCH_PAD + 'ACDACDACDACD' + MISMATCH_PAD}
  write_fasta(proteome_path, proteins)

  df = Matcher(
    query=[query],
    proteome_file=str(proteome_path),
    max_mismatches=1,
    k=3,
    preprocessed_files_path=str(tmp_path),
    sequence_version=False,
    output_format='dataframe'
  ).match()

  locations = list(zip(df['Protein ID'], df['Index start']))
  assert len(locations) == len(set(locations))
  assert reported_hits(df) == brute_force_hits(proteins, [query], 1)
  assert df.height > 1  # more than one genuine location, so dedupe is actually load-bearing


def test_mismatch_miss_row_shape(tmp_path):
  # A query with no window inside the threshold yields exactly one row with the match
  # columns nulled out -- and Mutated Positions as the literal '[]', not null. The
  # matching query in the same batch proves the miss is not a batch-wide failure.
  proteome_path = tmp_path / 'proteome.fasta'
  write_fasta(proteome_path, {'P1': MISMATCH_PAD + 'ACDEFGHIK' + MISMATCH_PAD})

  df = Matcher(
    query=['ACDEFGHIK', 'YTYTYTYTY'],
    proteome_file=str(proteome_path),
    max_mismatches=2,
    k=3,
    preprocessed_files_path=str(tmp_path),
    sequence_version=False,
    output_format='dataframe'
  ).match()

  miss = df.filter(pl.col('Query Sequence') == 'YTYTYTYTY')
  assert miss.height == 1
  row = miss.to_dicts()[0]
  for column in ('Matched Sequence', 'Protein ID', 'Mismatches', 'Index start', 'Index end'):
    assert row[column] is None, column
  assert row['Mutated Positions'] == '[]'
  assert df.filter(pl.col('Query Sequence') == 'ACDEFGHIK')['Mismatches'].to_list() == [0]


def test_mismatches_and_indels_mutually_exclusive_raises():
  # Mismatch and indel search are separate engines with different edit-count columns;
  # asking for both is undefined and must fail in the constructor.
  with pytest.raises(ValueError, match='mutually exclusive'):
    Matcher(query=['ACDEFGHIK'], proteome_file='unused.fasta', max_mismatches=1, max_indels=1)


def test_counts_only_with_best_match_raises():
  # best_match collapses to one row per query by tie-break; counts_only replaces the rows
  # with per-protein tallies. There is no coherent combined output, so it must raise.
  with pytest.raises(ValueError, match='counts_only is not supported with best_match'):
    Matcher(query=['ACDEFGHIK'], proteome_file='unused.fasta', counts_only=True, best_match=True)


@pytest.mark.xfail(reason="engine walks past protein boundary into the next record", strict=False)
def test_mismatch_does_not_cross_protein_boundary(tmp_path):
  # DEFGHI exists only in the CONCATENATION of P1 (ABCDEF) and P2 (GHIKLM); neither
  # protein contains it, nor anything within 1 mismatch of it. The mismatch path resolves
  # k-mer positions against the flat sequence with no per-protein end guard, so it reports
  # a bogus 0-mismatch hit in P1 at 4..9 -- past P1's 6 residues.
  proteome_path = tmp_path / 'proteome.fasta'
  proteins = {'P1': 'ABCDEF', 'P2': 'GHIKLM'}
  write_fasta(proteome_path, proteins)

  df = Matcher(
    query=['DEFGHI'],
    proteome_file=str(proteome_path),
    max_mismatches=1,
    k=3,
    preprocessed_files_path=str(tmp_path),
    sequence_version=False,
    output_format='dataframe'
  ).match()

  assert brute_force_hits(proteins, ['DEFGHI'], 1) == set()  # the oracle sees no hit
  for row in df.filter(pl.col('Matched Sequence').is_not_null()).iter_rows(named=True):
    assert row['Index end'] <= len(proteins[row['Protein ID']]), row
  assert df['Matched Sequence'].is_null().all()
