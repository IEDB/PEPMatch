import itertools
import random

import pytest

from pepmatch import matcher
from pepmatch.matcher import _indel_placements, format_indel_positions

AMINO_ACIDS = 'ACDEFGHIKLMNPQRSTVWY'


def _reference_placements(query, matched):
  """The original implementation, kept as the reference: try every combination of edit
  positions, rebuild the string, compare. Slow (L^n per call) but obviously correct, so
  the fast version must agree with it exactly -- placements, and their order."""
  L, M = len(query), len(matched)
  n = abs(M - L)
  if n == 0:
    return []
  out = []
  if M < L:
    for combo in itertools.combinations(range(1, L - 1), n):
      if ''.join(query[i] for i in range(L) if i not in combo) == matched:
        out.append(tuple((i + 1, query[i]) for i in combo))
  else:
    for combo in itertools.combinations(range(M), n):
      if ''.join(matched[i] for i in range(M) if i not in combo) != query:
        continue
      placement = []
      for k in combo:
        p = sum(1 for x in range(k) if x not in combo) + 1
        if p == 1 or p == L + 1:
          placement = None
          break
        placement.append((p, matched[k]))
      if placement:
        out.append(tuple(placement))
  return out


def _edited_pairs(query, alphabet):
  """Every (query, matched) pair one or two deletions or insertions away from `query`,
  including terminal edits, which the annotation must reject identically."""
  L = len(query)
  for n in (1, 2):
    for pos in itertools.combinations(range(L), n):
      yield query, ''.join(c for i, c in enumerate(query) if i not in pos)
    for pos in itertools.combinations(range(L + n), n):
      for residues in itertools.product(alphabet, repeat=n):
        it, extra = iter(query), iter(residues)
        yield query, ''.join(next(extra) if i in pos else next(it) for i in range(L + n))


def _assert_same(query, matched, monkeypatch):
  assert _indel_placements(query, matched) == _reference_placements(query, matched), \
    (query, matched)
  fast = format_indel_positions(query, matched)
  monkeypatch.setattr(matcher, '_indel_placements', _reference_placements)
  slow = format_indel_positions(query, matched)
  monkeypatch.undo()
  assert fast == slow, (query, matched)


@pytest.mark.parametrize('alphabet,max_len', [('AB', 8), ('ABC', 6)])
def test_exhaustive_small_alphabet(alphabet, max_len, monkeypatch):
  """A two- or three-letter alphabet makes every query a homopolymer or periodic repeat
  somewhere, which is where placements are ambiguous and ranges/chunks get reported."""
  pairs = set()
  for L in range(3, max_len + 1):
    for query in map(''.join, itertools.product(alphabet, repeat=L)):
      pairs.update(_edited_pairs(query, alphabet))
  for query, matched in pairs:
    _assert_same(query, matched, monkeypatch)


@pytest.mark.parametrize('seed', range(5))
def test_random_peptides(seed, monkeypatch):
  """Realistic lengths (6-25) with a share of low-complexity queries."""
  rng = random.Random(seed)
  for _ in range(400):
    L = rng.randint(6, 25)
    pool = 'AAKKL' if rng.random() < 0.3 else AMINO_ACIDS
    query = ''.join(rng.choice(pool) for _ in range(L))
    n = rng.choice((1, 2))
    if rng.random() < 0.5:
      pos = rng.sample(range(L), n)
      matched = ''.join(c for i, c in enumerate(query) if i not in pos)
    else:
      chars = list(query)
      for _ in range(n):
        chars.insert(rng.randint(0, len(chars)), rng.choice(AMINO_ACIDS))
      matched = ''.join(chars)
    _assert_same(query, matched, monkeypatch)


def test_unrelated_and_short_pairs(monkeypatch):
  """Every pair over {A, B} one or two residues apart in length, related by indels or not,
  down to the empty string: both versions must agree on which have no valid placement."""
  for q_len in range(0, 7):
    for m_len in {q_len - 2, q_len - 1, q_len + 1, q_len + 2} - {-2, -1}:
      for query in map(''.join, itertools.product('AB', repeat=q_len)):
        for matched in map(''.join, itertools.product('AB', repeat=m_len)):
          _assert_same(query, matched, monkeypatch)


@pytest.mark.parametrize('query,matched,expected', [
  ('ABCDEF', 'ABCDEF', '[]'),            # exact match
  ('CAAAD', 'CAAD', 'd: A[2,4]'),        # homopolymer: the range of valid positions
  ('ABABAB', 'ABAB', 'd: BA[2], d: AB[3], d: BA[4]'),   # periodic: residues change along the range
])
def test_documented_examples(query, matched, expected):
  assert format_indel_positions(query, matched) == expected
