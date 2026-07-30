import pytest
import polars as pl
import polars.testing as plt
from pathlib import Path
from pepmatch import Matcher

@pytest.fixture
def proteome_path() -> Path:
  return Path(__file__).parent / 'data' / 'proteome.fasta'

@pytest.fixture
def expected_path() -> Path:
  return Path(__file__).parent / 'data' / 'discontinuous_expected.csv'

@pytest.fixture
def query() -> list:
  return [
    'L354, V420, G461, Q468, E486, K499, D501, M503, G509',
    'T819, L822, A835, F840, S871, Y884, P886',
    'S2760, V2763, E2773, D2805, T2819, S2831, E2844, R2852, L2863'
  ]

def test_discontinuous_search(proteome_path, query, expected_path):
  df = Matcher(
    query=query,
    proteome_file=proteome_path,
    max_mismatches=0,
    output_format='dataframe'
  ).match()

  df = df.sort('Query Sequence')
  expected_df = pl.read_csv(expected_path).sort('Query Sequence')
  plt.assert_series_equal(df['Protein ID'], expected_df['Protein ID'])

def test_mixed_unmatched_linear_and_matched_discontinuous(proteome_path):
  unmatched_linear = 'WWWWWWWWW'
  matched_discontinuous = 'L354, V420, G461, Q468, E486, K499, D501, M503, G509'

  df = Matcher(
    query=[unmatched_linear, matched_discontinuous],
    proteome_file=proteome_path,
    max_mismatches=0,
    output_format='dataframe'
  ).match()

  assert df.schema['Matched Sequence'] == pl.String
  assert df.filter(pl.col('Query Sequence') == unmatched_linear)['Matched Sequence'].item() is None
  assert (
    df.filter(pl.col('Query Sequence') == matched_discontinuous)['Matched Sequence'].item()
    == matched_discontinuous
  )


def _pairs(epitope: str) -> list:
  """Parse 'N1, A3, L5' the way the Matcher does: (residue, 1-based position) pairs."""
  return [(x[0], int(x[1:])) for x in epitope.split(', ')]


def write_proteome(path: Path, proteins: dict) -> None:
  """Write a minimal SwissProt-style FASTA. SV is absent, so accession ACC reports as 'ACC.1'."""
  path.write_text(''.join(
    f'>sp|{acc}|{acc}_TEST Test protein OS=Homo sapiens OX=9606 GN={acc} PE=1 SV=1\n{seq}\n'
    for acc, seq in proteins.items()
  ))


def search(tmp_path: Path, proteins: dict, query, **kwargs):
  """Search a throwaway proteome, keeping the .pepidx index inside tmp_path."""
  proteome_path = tmp_path / 'proteome.fasta'
  write_proteome(proteome_path, proteins)
  return Matcher(
    query=query,
    proteome_file=str(proteome_path),
    preprocessed_files_path=str(tmp_path),
    output_format='dataframe',
    **kwargs
  ).match()


def oracle(proteins: dict, epitope: str, max_mismatches: int) -> dict:
  """Brute-force reference: a protein satisfies a discontinuous epitope when every listed
  position is within it and at most `max_mismatches` of them hold the wrong residue.
  Returns {'ACC.1': (mismatches, matched sequence, mutated positions)}."""
  pairs = _pairs(epitope)
  hits = {}
  for acc, seq in proteins.items():
    if any(pos < 1 or pos > len(seq) for _, pos in pairs):
      continue
    wrong = [pos for residue, pos in pairs if seq[pos - 1] != residue]
    if len(wrong) > max_mismatches:
      continue
    hits[f'{acc}.1'] = (
      len(wrong),
      ', '.join(f'{seq[pos - 1]}{pos}' for _, pos in pairs),
      f'[{", ".join(str(pos) for pos in wrong)}]',
    )
  return hits


def rows_for(df, query: str) -> dict:
  """The engine's hits for one query, keyed like `oracle` for direct comparison."""
  return {
    row['Protein ID']: (row['Mismatches'], row['Matched Sequence'], row['Mutated Positions'])
    for row in df.filter(pl.col('Query Sequence') == query).iter_rows(named=True)
  }


def assert_miss(df, query: str) -> None:
  """A miss is exactly one row per query with the null hit fields and '[]' mutated positions."""
  rows = df.filter(pl.col('Query Sequence') == query)
  assert rows.height == 1, f'{query!r} should produce exactly one row'
  row = rows.row(0, named=True)
  for column in ('Matched Sequence', 'Protein ID', 'Mismatches', 'Index start', 'Index end'):
    assert row[column] is None, f'{column} should be null for missed query {query!r}'
  assert row['Mutated Positions'] == '[]'


def test_discontinuous_exact_hit_reports_listed_positions(tmp_path):
  # Q2/L4/M7 are the residues of P1 at exactly those 1-based positions; P2 has none of them.
  # Index start/end are the first and last listed positions, so an off-by-one or a
  # contiguous-window interpretation would not report 2..7.
  proteins = {'P1': 'NQALLKMWY', 'P2': 'WWWWWWWWW'}
  epitope = 'Q2, L4, M7'
  df = search(tmp_path, proteins, [epitope])

  assert rows_for(df, epitope) == oracle(proteins, epitope, 0)
  assert df.height == 1
  row = df.row(0, named=True)
  assert row['Protein ID'] == 'P1.1'
  assert row['Mismatches'] == 0
  assert row['Matched Sequence'] == epitope
  assert row['Mutated Positions'] == '[]'
  assert (row['Index start'], row['Index end']) == (2, 7)


def test_discontinuous_tolerated_mismatch_echoes_protein_residue(tmp_path):
  # P1 holds C, not A, at listed position 3: reportable within a budget of one wrong
  # residue, and Matched Sequence must echo what the protein actually has there.
  proteins = {'P1': 'NQCLLKMWY'}
  epitope = 'N1, A3, L5'
  df = search(tmp_path, proteins, [epitope], max_mismatches=1)

  assert rows_for(df, epitope) == oracle(proteins, epitope, 1)
  row = df.row(0, named=True)
  assert row['Protein ID'] == 'P1.1'
  assert row['Mismatches'] == 1
  assert row['Matched Sequence'] == 'N1, C3, L5'
  assert row['Mutated Positions'] == '[3]'
  assert (row['Index start'], row['Index end']) == (1, 5)

  strict = search(tmp_path, proteins, [epitope], max_mismatches=0)
  assert oracle(proteins, epitope, 0) == {}
  assert_miss(strict, epitope)


def test_discontinuous_mismatch_budget_counts_every_listed_position(tmp_path):
  # P1 is wrong at two listed positions: over budget at 1, reported at 2 with both
  # protein positions named, in listed order.
  proteins = {'P1': 'NQCLYKMWY'}
  epitope = 'N1, A3, L5'

  assert oracle(proteins, epitope, 1) == {}
  assert_miss(search(tmp_path, proteins, [epitope], max_mismatches=1), epitope)

  df = search(tmp_path, proteins, [epitope], max_mismatches=2)
  assert rows_for(df, epitope) == oracle(proteins, epitope, 2)
  row = df.row(0, named=True)
  assert row['Mismatches'] == 2
  assert row['Matched Sequence'] == 'N1, C3, Y5'
  assert row['Mutated Positions'] == '[3, 5]'


def test_discontinuous_reports_every_satisfying_protein(tmp_path):
  # Three of the four proteins satisfy the epitope, including one whose length stops
  # exactly at the last listed position: the scan is not first-match-only.
  proteins = {
    'P1': 'NQALLKMWY',
    'P2': 'NYAYLYYYY',
    'P3': 'WWWWWWWWW',
    'P4': 'NGAGL',
  }
  epitope = 'N1, A3, L5'
  df = search(tmp_path, proteins, [epitope])

  expected = oracle(proteins, epitope, 0)
  assert set(expected) == {'P1.1', 'P2.1', 'P4.1'}
  assert rows_for(df, epitope) == expected
  assert df.height == len(expected)


def test_discontinuous_out_of_range_positions_miss_cleanly(tmp_path):
  # 50 is past every protein's end, and 0 is not a 1-based position. Both disqualify a
  # protein rather than raising or being skipped -- skipping the 0 would match both proteins.
  proteins = {'P1': 'NQALLKMWY', 'P2': 'NQALL'}
  past_end = 'N1, A3, L50'
  zero = 'N0, A3, L5'
  df = search(tmp_path, proteins, [past_end, zero])

  assert oracle(proteins, past_end, 0) == {}
  assert oracle(proteins, zero, 0) == {}
  assert df.height == 2
  assert_miss(df, past_end)
  assert_miss(df, zero)


def test_discontinuous_position_past_protein_end_is_not_read_from_next_protein(tmp_path):
  # P1 is 4 residues long, so listed position 5 lands on P2's first residue in the
  # concatenated index. Bounds are per protein: a miss, not a 'C5' hit inside P1.
  proteins = {'P1': 'AAAA', 'P2': 'CCCCGGGG'}
  epitope = 'A1, C5'
  df = search(tmp_path, proteins, [epitope])

  assert oracle(proteins, epitope, 0) == {}
  assert_miss(df, epitope)


def test_discontinuous_listed_position_order_is_preserved(tmp_path):
  # Listed out of ascending order, Index start/end are the FIRST and LAST listed
  # positions -- 5 and 3 -- not the min and max, and Matched Sequence keeps that order.
  proteins = {'P1': 'NQALLKMWY'}
  epitope = 'L5, N1, A3'
  df = search(tmp_path, proteins, [epitope])

  assert rows_for(df, epitope) == oracle(proteins, epitope, 0)
  row = df.row(0, named=True)
  assert row['Matched Sequence'] == epitope
  assert (row['Index start'], row['Index end']) == (5, 3)
  assert row['Index start'] > row['Index end']


def test_peptide_and_residue_position_queries_keep_their_own_semantics(tmp_path):
  # 'NQALL' is letters only: a linear peptide, matched as a substring at offset 3.
  # 'G1, N3, A5' parses as residue/position pairs and is anchored to those positions.
  # 'N1, A3, L5' lists the peptide's residues at the positions they do NOT occupy, so
  # it must miss -- position semantics, not a shifted substring scan.
  proteins = {'P1': 'GGNQALLKMWY'}
  peptide = 'NQALL'
  anchored = 'G1, N3, A5'
  shifted = 'N1, A3, L5'
  df = search(tmp_path, proteins, [peptide, anchored, shifted])

  assert df.schema['Matched Sequence'] == pl.String
  assert df.height == 3

  linear = df.filter(pl.col('Query Sequence') == peptide).row(0, named=True)
  start = proteins['P1'].index(peptide) + 1
  assert linear['Matched Sequence'] == peptide
  assert (linear['Index start'], linear['Index end']) == (start, start + len(peptide) - 1)

  assert rows_for(df, anchored) == oracle(proteins, anchored, 0)
  assert rows_for(df, anchored) == {'P1.1': (0, anchored, '[]')}
  assert_miss(df, shifted)


def test_discontinuous_with_max_indels_raises():
  # A position-anchored epitope has no contiguous window for the indel engine to seed.
  with pytest.raises(ValueError, match='discontinuous'):
    Matcher(query=['N1, A3, L5'], proteome_file='unused.fasta', max_indels=1)


def test_single_residue_position_pair_is_a_one_residue_epitope(tmp_path):
  # 'N1' is one residue/position pair, so it is anchored, not a peptide: it hits P1,
  # whose first residue is N, and not P2, which has an N at position 2 instead.
  proteins = {'P1': 'NQALLKMWY', 'P2': 'WNQALL'}
  df = search(tmp_path, proteins, ['N1'])

  assert oracle(proteins, 'N1', 0) == {'P1.1': (0, 'N1', '[]')}
  assert rows_for(df, 'N1') == oracle(proteins, 'N1', 0)
  assert df.height == 1
  row = df.row(0, named=True)
  assert row['Matched Sequence'] == 'N1'
  assert row['Index start'] == row['Index end'] == 1
