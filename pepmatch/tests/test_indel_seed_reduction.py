"""Recall guard for the reduced-seed indel search.

Indel search looks up only the n+1 rarest disjoint query tiles instead of the full tiling.
Complete recall does not depend on which disjoint seeds are chosen (n edits damage at most n
of the n+1, so one always survives inside any match and anchors it), but it DOES depend on
extend_bidirectional reconstructing the full alignment from whatever single seed survives --
including a purely interior one, and including the case where the surviving seed was chosen
for rarity rather than position.

These cases are under-exercised by the generic property tests (short queries, few tiles, no
reduction). This test forces them: mixed-length queries so long ones carry many tiles that
get dropped, interior edit placement so the surviving seed is often not at a boundary, and a
low-complexity protein so the rarest-seed choice actively routes around common tiles. It then
asserts the engine's match set equals the independent brute-force oracle -- if the reduction
ever dropped a match, actual != expected.
"""
import random

import pytest
from pepmatch import Matcher

from indel_brute_force import brute_force_search as bf1
from indel2_brute_force import brute_force_search as bf2

AMINO_ACIDS = 'ACDEFGHIKLMNPQRSTVWY'


def _random_sequence(rng, length):
  return ''.join(rng.choice(AMINO_ACIDS) for _ in range(length))


def _repeat_rich(rng):
  """Homopolymer / periodic tracts: where tiles are least unique, so the rarest-seed choice
  has the most to do and a common seed would explode into many candidates."""
  blocks = ['AAAAAAAA', 'EKEKEKEK', 'QQQQQQ', 'LSLSLSLS', 'GGGGGG', 'PAPAPAPA']
  rng.shuffle(blocks)
  return ''.join(blocks) + _random_sequence(rng, rng.randint(10, 20))


def _inject_interior(rng, seq, n, kind):
  """Apply n edits at INTERIOR positions (never the first/last residue, which the engine
  blocks as terminal), so a real <=n-indel match to `seq` exists and the boundary tiles of
  the resulting query are the ones most likely damaged."""
  s = seq
  for _ in range(n):
    if kind == 'deletion' and len(s) > 3:
      d = rng.randrange(1, len(s) - 1)          # interior delete
      s = s[:d] + s[d + 1:]
    elif kind == 'insertion':
      i = rng.randrange(1, len(s))              # interior insert
      s = s[:i] + rng.choice(AMINO_ACIDS) + s[i:]
  return s


def _oracle(qseq, pseq, n):
  return bf1(qseq, pseq) if n == 1 else bf2(qseq, pseq, n)


@pytest.mark.parametrize('n', [1, 2])
@pytest.mark.parametrize('seed', range(25))
def test_reduced_seed_indel_matches_oracle(tmp_path, n, seed):
  rng = random.Random((seed << 1) | (n - 1))

  proteins = {f'P{i}': _random_sequence(rng, rng.randint(45, 80)) for i in range(3)}
  proteins['R0'] = _repeat_rich(rng)            # low-complexity: exercises rarest choice
  proteome_path = tmp_path / 'proteome.fasta'
  proteome_path.write_text(''.join(f'>{pid}\n{seq}\n' for pid, seq in proteins.items()))

  # Mixed query lengths: a couple of short ones pin k small so the LONG queries carry many
  # tiles (most of which the reduction drops), which is the regime this test exists for.
  queries = []
  for i in range(10):
    pid = rng.choice(list(proteins))
    pseq = proteins[pid]
    qlen = rng.choice([8, 9, 10] + [18, 20, 22, 24, 25])
    if len(pseq) < qlen + 2:
      qlen = min(qlen, len(pseq) - 2)
    start = rng.randrange(0, len(pseq) - qlen + 1)
    base = pseq[start:start + qlen]
    mode = rng.choice(['deletion', 'insertion', 'random'])
    q = _random_sequence(rng, qlen) if mode == 'random' else _inject_interior(rng, base, n, mode)
    queries.append((f'q{i}', q))

  query_path = tmp_path / 'queries.fasta'
  query_path.write_text(''.join(f'>{qid}\n{seq}\n' for qid, seq in queries))

  df = Matcher(
    query=str(query_path),
    proteome_file=str(proteome_path),
    max_indels=n,
    preprocessed_files_path=str(tmp_path),
    output_format='dataframe',
  ).match()

  actual = set()
  for row in df.iter_rows(named=True):
    if row['Matched Sequence'] is not None:
      actual.add((row['Query Sequence'], row['Protein ID'], row['Matched Sequence']))

  expected = set()
  for _, qseq in queries:
    for pid, pseq in proteins.items():
      for _, matched in _oracle(qseq, pseq, n):
        expected.add((qseq, f'{pid}.1', matched))

  assert actual == expected
