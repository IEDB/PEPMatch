import random
from pathlib import Path

import pytest
from Bio import SeqIO
from pepmatch import Matcher

# Guards the auto-k mismatch fix in matcher.py: when k is UNSPECIFIED, Matcher must
# derive the pigeonhole-optimal k = max(2, min_len // (m + 1)) so that a seed-based
# search returns the COMPLETE set of matches. The old code hardcoded k=5, which is
# too large for short peptides / high m and silently dropped real matches. These
# tests pin completeness (and soundness) against an independent brute-force oracle.

AMINO_ACIDS = 'ACDEFGHIKLMNPQRSTVWY'
LENGTHS = (8, 9, 10, 11, 12, 15, 20)
MISMATCHES = (1, 2, 3)


@pytest.fixture
def proteome_path() -> Path:
  return Path(__file__).parent / 'data' / 'proteome.fasta'


def _load_proteins(proteome_path):
  """Load (accession, sequence) pairs. The accession is exactly how PEPMatch labels
  a protein in its output ('<accession>.<version>' in the Protein ID column), so
  keying the oracle by accession lets us compare match LOCATIONS, not just sequences."""
  proteins = []
  for record in SeqIO.parse(str(proteome_path), 'fasta'):
    accession = record.id.split('|')[1] if '|' in record.id else record.id
    proteins.append((accession, str(record.seq).upper()))
  return proteins


def _hamming(a, b):
  return sum(1 for x, y in zip(a, b) if x != y)


def brute_force_mismatch(peptide, proteins, max_mismatches):
  """Independent oracle: every proteome window of len(peptide) whose Hamming
  distance to the peptide is <= max_mismatches. No seeding, no k-mers -- just the
  raw definition of a mismatch match, so it cannot share a blind spot with the
  seed-based implementation. Returns {(accession, start_0based, window)}."""
  length = len(peptide)
  hits = set()
  for accession, seq in proteins:
    for i in range(len(seq) - length + 1):
      window = seq[i:i + length]
      if _hamming(peptide, window) <= max_mismatches:
        hits.add((accession, i, window))
  return hits


def _matcher_rows(df):
  """Accepted (non-null) matches from Matcher output, canonicalized into the
  oracle's key space: strip the sequence-version suffix from Protein ID to recover
  the accession, and convert the 1-based Index start to a 0-based offset. Returns a
  list (not a set) so callers can detect duplicate emission of the same match."""
  rows = []
  for row in df.iter_rows(named=True):
    if row['Matched Sequence'] is None:
      continue
    accession = row['Protein ID'].rsplit('.', 1)[0]
    rows.append((accession, row['Index start'] - 1, row['Matched Sequence']))
  return rows


def _make_query(rng, proteins, length, m):
  """Take a random proteome window of the given length and inject EXACTLY m
  mismatches at distinct positions, each replaced by a different amino acid. The
  source window is therefore always a real match at Hamming distance m."""
  candidates = [(acc, seq) for acc, seq in proteins if len(seq) >= length]
  accession, seq = rng.choice(candidates)
  start = rng.randrange(0, len(seq) - length + 1)
  window = list(seq[start:start + length])
  for p in rng.sample(range(length), m):
    window[p] = rng.choice([a for a in AMINO_ACIDS if a != window[p]])
  return ''.join(window)


@pytest.mark.parametrize('m', MISMATCHES)
def test_auto_k_search_is_complete_vs_brute_force(proteome_path, tmp_path, m):
  # Seeded so a red run reproduces exactly; seed differs per m to vary the corpus.
  rng = random.Random(20240717 + m)
  proteins = _load_proteins(proteome_path)
  seq_by_acc = dict(proteins)

  with_matches = 0
  for length in LENGTHS:
    if length <= m:
      continue
    for _ in range(12):
      peptide = _make_query(rng, proteins, length, m)
      expected = brute_force_mismatch(peptide, proteins, m)
      if not expected:
        continue  # a query with no oracle match cannot exercise completeness

      df = Matcher(
        query=[peptide],
        proteome_file=str(proteome_path),
        max_mismatches=m,
        preprocessed_files_path=str(tmp_path),
      ).match()

      length = len(peptide)
      rows = _matcher_rows(df)
      # A reported hit whose Matched Sequence does not equal the protein's own
      # window at that offset is the pre-existing concatenated-index boundary
      # artifact (a window running off the end of one protein into the next). That
      # is a soundness quirk unrelated to k selection, so it is excluded from the
      # in-protein comparison below; it can only ever ADD matches, never drop one.
      in_protein = [
        t for t in rows if seq_by_acc[t[0]][t[1]:t[1] + length] == t[2]
      ]
      actual = set(in_protein)

      # COMPLETENESS (the contract the fix guards): every real proteome window
      # within m mismatches must be reported -- nothing silently dropped. Set
      # equality also pins SOUNDNESS: no genuine window is over-reported.
      assert actual == expected, (
        f'auto-k disagreed with brute force for {peptide!r} (m={m}); '
        f'dropped={sorted(expected - actual)} spurious={sorted(actual - expected)}'
      )
      # No duplicate emission of a genuine match (a bare set would absorb dupes).
      assert len(in_protein) == len(expected), (
        f'auto-k emitted {len(in_protein)} in-protein rows but oracle has '
        f'{len(expected)} distinct matches for {peptide!r} (m={m})'
      )
      with_matches += 1

  # The suite is only meaningful on peptides that actually have matches; refuse to
  # pass vacuously if the corpus somehow degenerated to all-misses.
  assert with_matches >= 50


def test_auto_k_recovers_match_that_hardcoded_k5_dropped(proteome_path, tmp_path):
  # Focused regression for the exact failure mode the matcher.py fix repairs.
  #
  # 'MSKRALIIA' is the proteome window 'MSKRVLIIA' (accession Q02004, index 0)
  # carrying a single mismatch at the centre. A length-9 / 1-mismatch query needs
  # two disjoint seeds -- i.e. k <= 4 -- for the pigeonhole guarantee to hold.
  # The OLD matcher hardcoded k=5 whenever k was unspecified: with a 5-mer seed no
  # sub-window of the query avoids the central mismatch, so the seed lookup finds
  # no candidate and this real match is SILENTLY DROPPED. The fix derives
  # k = max(2, 9 // (1 + 1)) = 4, which recovers it.
  proteins = _load_proteins(proteome_path)
  peptide = 'MSKRALIIA'
  source_match = ('Q02004', 0, 'MSKRVLIIA')

  expected = brute_force_mismatch(peptide, proteins, 1)
  assert source_match in expected  # the oracle confirms this is a genuine match

  auto_k = set(_matcher_rows(
    Matcher(
      query=[peptide],
      proteome_file=str(proteome_path),
      max_mismatches=1,
      preprocessed_files_path=str(tmp_path),
    ).match()
  ))
  hardcoded_k5 = set(_matcher_rows(
    Matcher(
      query=[peptide],
      proteome_file=str(proteome_path),
      max_mismatches=1,
      k=5,
      preprocessed_files_path=str(tmp_path),
    ).match()
  ))

  # auto-k is complete and finds the match; the old hardcoded k=5 drops it.
  assert source_match in auto_k
  assert auto_k == expected
  assert source_match not in hardcoded_k5
