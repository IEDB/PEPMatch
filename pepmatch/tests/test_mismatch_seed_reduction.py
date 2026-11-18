"""Recall guard for the reduced-seed mismatch search.

When a query holds n+1 disjoint k-mers, mismatch search looks up only the n+1 rarest of its
tiles instead of every tile (or every offset, when k does not divide the length). Complete
recall does not depend on which disjoint seeds are chosen: n mismatches damage at most n of
the n+1, so one always survives inside any match and anchors it. It DOES depend on the search
verifying the whole peptide from whatever single seed survives, including an interior one
chosen for rarity rather than position.

The other mismatch tests use one query length at its optimal k, where there is nothing to
drop. This test forces the reduction: mixed-length queries searched at a small k so long ones
carry many tiles, mismatches placed anywhere so the surviving seed is often interior, and
repeat-rich proteins so the rarest-seed choice routes around common tiles. It then asserts
the engine's matches equal an independent brute-force scan, row for row, in both the search
and the counts-only path.

A query too short for the guarantee keeps the full seeding (every tile when k divides its
length, every offset otherwise), so it must still find exactly the matches that contain one
of those seeds unchanged. The brute-force scan computes that set directly, which pins the
best-effort behaviour as well as the guaranteed one.
"""
import random
from collections import Counter

import pytest
from pepmatch import Matcher

AMINO_ACIDS = 'ACDEFGHIKLMNPQRSTVWY'
MISMATCHES = (1, 2, 3)
KS = (0, 2, 3, 4)   # 0 lets Matcher choose k from the shortest query


def _random_sequence(rng, length):
  return ''.join(rng.choice(AMINO_ACIDS) for _ in range(length))


def _repeat_rich(rng):
  """Homopolymer and periodic tracts: where tiles are least unique, so the rarest-seed
  choice has the most to do."""
  blocks = ['AAAAAAAA', 'EKEKEKEK', 'QQQQQQ', 'LSLSLSLS', 'GGGGGG', 'PAPAPAPA']
  rng.shuffle(blocks)
  return ''.join(blocks) + _random_sequence(rng, rng.randint(10, 20))


def _proteome(rng):
  """(id, sequence) pairs: repeat-rich proteins, ordinary ones, and two shorter than most
  queries so windows running past a protein's end are exercised."""
  sequences = [_repeat_rich(rng) for _ in range(5)]
  sequences += [_random_sequence(rng, rng.randint(40, 120)) for _ in range(12)]
  sequences += [_random_sequence(rng, 9), _random_sequence(rng, 12)]
  return [(f'p{i + 1}', seq) for i, seq in enumerate(sequences)]


def _queries(rng, proteins):
  """Proteome windows of mixed length carrying 0-3 substitutions anywhere, plus a few
  random peptides that should mostly miss."""
  queries = []
  for _ in range(40):
    length = rng.randint(8, 25)
    _, source = rng.choice([p for p in proteins if len(p[1]) >= length])
    start = rng.randrange(len(source) - length + 1)
    window = list(source[start:start + length])
    for position in rng.sample(range(length), rng.randint(0, 3)):
      window[position] = rng.choice(AMINO_ACIDS)
    queries.append(''.join(window))
  queries += [_random_sequence(rng, rng.randint(8, 25)) for _ in range(6)]
  return queries


def _brute_force(peptide, proteins, max_mismatches):
  """Every proteome window within the budget, by the definition of a mismatch match: no
  seeds, no k-mers. Rows are (protein id, 1-based start, window, mismatches, positions)."""
  length = len(peptide)
  hits = set()
  for protein_id, seq in proteins:
    for i in range(len(seq) - length + 1):
      window = seq[i:i + length]
      positions = [j + 1 for j in range(length) if peptide[j] != window[j]]
      if len(positions) <= max_mismatches:
        hits.add((protein_id, i + 1, window, len(positions), str(positions)))
  return hits


def _expected(peptide, proteins, max_mismatches, k):
  """What the search must report at k: every match when the peptide holds
  max_mismatches + 1 disjoint k-mers, otherwise the matches containing a seed of the full
  seeding unchanged."""
  length = len(peptide)
  hits = _brute_force(peptide, proteins, max_mismatches)
  if length < k:
    return set()
  offsets = range(0, length - k + 1, k if length % k == 0 else 1)
  reachable = {
    hit for hit in hits
    if any(hit[2][j:j + k] == peptide[j:j + k] for j in offsets)
  }
  if length >= k * (max_mismatches + 1):
    assert reachable == hits   # the pigeonhole guarantee itself
  return reachable


def _engine_rows(df):
  """Reported hits per Query ID, as lists so a row emitted twice is caught."""
  rows = {}
  for row in df.iter_rows(named=True):
    if row['Matched Sequence'] is not None:
      rows.setdefault(row['Query ID'], []).append((
        row['Protein ID'], row['Index start'], row['Matched Sequence'],
        row['Mismatches'], row['Mutated Positions'],
      ))
  return rows


def _setup(seed, tmp_path):
  rng = random.Random(20261114 + seed)
  proteins = _proteome(rng)
  fasta = tmp_path / 'proteome.fasta'
  fasta.write_text(''.join(f'>{protein_id}\n{seq}\n' for protein_id, seq in proteins))
  return proteins, fasta, _queries(rng, proteins)


def _effective_k(k, queries, max_mismatches):
  return k or max(2, min(len(q) for q in queries) // (max_mismatches + 1))


@pytest.mark.parametrize('seed', range(8))
def test_mismatch_search_equals_brute_force(seed, tmp_path):
  proteins, fasta, queries = _setup(seed, tmp_path)
  compared = dropped = unguaranteed = 0
  for m in MISMATCHES:
    for k in KS:
      df = Matcher(
        query=queries, proteome_file=str(fasta), max_mismatches=m, k=k,
        preprocessed_files_path=str(tmp_path), sequence_version=False,
      ).match()
      rows = _engine_rows(df)
      used_k = _effective_k(k, queries, m)
      for i, peptide in enumerate(queries):
        actual = rows.get(str(i + 1), [])
        expected = _expected(peptide, proteins, m, used_k)
        assert len(actual) == len(set(actual)), f'duplicate rows for {peptide} (m={m}, k={used_k})'
        assert set(actual) == expected, f'{peptide} (m={m}, k={used_k})'
        compared += len(expected)
        if len(peptide) >= used_k * (m + 1):
          dropped += len(peptide) // used_k > m + 1
        else:
          unguaranteed += 1
  # The corpus must contain real matches and exercise both the reduction and the full seeding.
  assert compared > 100 and dropped > 100 and unguaranteed > 20, (compared, dropped, unguaranteed)


@pytest.mark.parametrize('seed', range(4))
def test_counts_only_equals_brute_force(seed, tmp_path):
  proteins, fasta, queries = _setup(seed, tmp_path)
  for m in MISMATCHES:
    for k in KS:
      df = Matcher(
        query=queries, proteome_file=str(fasta), max_mismatches=m, k=k,
        preprocessed_files_path=str(tmp_path), counts_only=True,
      ).match()
      actual = {(row['Query ID'], row['Mismatches']): row['Count'] for row in df.iter_rows(named=True)}
      used_k = _effective_k(k, queries, m)
      for i, peptide in enumerate(queries):
        expected = Counter(hit[3] for hit in _expected(peptide, proteins, m, used_k))
        got = {level: count for (query_id, level), count in actual.items() if query_id == str(i + 1)}
        assert got == dict(expected), f'{peptide} (m={m}, k={used_k})'


@pytest.mark.parametrize('m', MISMATCHES)
def test_a_match_never_crosses_into_the_neighbouring_protein(m, tmp_path):
  # Proteins are stored end to end. Each query below matches only if the search reads past
  # one protein's edge into its neighbour, so neither may be reported.
  proteins = [('p1', 'WWWWWMKTAYIAKQ'), ('p2', 'RQISFVKSHFSRQLEE'), ('p3', 'GGGGGGGG')]
  fasta = tmp_path / 'edges.fasta'
  fasta.write_text(''.join(f'>{protein_id}\n{seq}\n' for protein_id, seq in proteins))
  queries = [
    'MKTAYIAKQR',     # p1's last nine residues + p2's first: would run off p1's end
    'QRQISFVKSH',     # p1's last residue + p2's first nine: would start before p2
    'AKQRQISF',       # straddles the boundary in the middle
  ]
  for k in (2, 3):
    df = Matcher(
      query=queries, proteome_file=str(fasta), max_mismatches=m, k=k,
      preprocessed_files_path=str(tmp_path), sequence_version=False,
    ).match()
    rows = _engine_rows(df)
    for i, peptide in enumerate(queries):
      assert set(rows.get(str(i + 1), [])) == _expected(peptide, proteins, m, k), f'{peptide} (m={m}, k={k})'
