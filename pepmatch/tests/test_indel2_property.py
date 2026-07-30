import random
import pytest
from pepmatch import Matcher
from indel2_brute_force import brute_force_search

AMINO_ACIDS = 'ACDEFGHIKLMNPQRSTVWY'


def _random_sequence(rng, length):
  return ''.join(rng.choice(AMINO_ACIDS) for _ in range(length))


def _repeat_rich_protein(rng):
  """A protein built from homopolymer and periodic-repeat tracts -- the MSI /
  frameshift regime 2-indel search exists to catch, and exactly where 2-mer seeds are
  least unique. A purely random proteome never generates these, so without injecting
  one the engine's behavior on repeats (many candidate seeds per lookup, ambiguous
  indel placement) is never exercised."""
  blocks = ['AAAAAAAA', 'EKEKEKEK', 'QQQQQ', 'LSLSLSLS', 'GGGGGG']
  rng.shuffle(blocks)
  return 'WW' + ''.join(blocks) + 'WW'


# Named for the edit applied to the query string, not the resulting Indels label — e.g.
# _with_deletions removes query residues, so the protein has "extra" content there
# relative to the query, which the search reports as an insertion match.
def _with_deletions(rng, seq, n):
  for _ in range(n):
    d = rng.randrange(len(seq))
    seq = seq[:d] + seq[d + 1:]
  return seq


def _with_insertions(rng, seq, n):
  for _ in range(n):
    i = rng.randrange(1, len(seq))
    seq = seq[:i] + rng.choice(AMINO_ACIDS) + seq[i:]
  return seq


def _random_query(rng, proteins):
  protein = proteins[rng.choice(list(proteins))]
  qlen = rng.randint(8, 14)  # min 8 so a 2-deletion still yields length >= 6 (k=2 floor)
  start = rng.randrange(0, len(protein) - qlen + 1)
  base = protein[start:start + qlen]
  mode = rng.choice(['exact', 'del1', 'ins1', 'del2', 'ins2', 'random'])
  if mode == 'del1':
    return _with_deletions(rng, base, 1)
  if mode == 'ins1':
    return _with_insertions(rng, base, 1)
  if mode == 'del2':
    return _with_deletions(rng, base, 2)
  if mode == 'ins2':
    return _with_insertions(rng, base, 2)
  if mode == 'random':
    return _random_sequence(rng, qlen)
  return base


@pytest.mark.parametrize('seed', range(15))
def test_indel2_search_matches_brute_force_oracle(tmp_path, seed):
  rng = random.Random(seed)

  proteins = {f'P{i}': _random_sequence(rng, rng.randint(20, 40)) for i in range(4)}
  proteins['REP'] = _repeat_rich_protein(rng)  # homopolymer / periodic-repeat coverage
  proteome_path = tmp_path / 'proteome.fasta'
  proteome_path.write_text(
    ''.join(f'>{pid}\n{seq}\n' for pid, seq in proteins.items())
  )

  queries = [(f'q{i}', _random_query(rng, proteins)) for i in range(8)]
  query_path = tmp_path / 'queries.fasta'
  query_path.write_text(
    ''.join(f'>{qid}\n{seq}\n' for qid, seq in queries)
  )

  df = Matcher(
    query=str(query_path),
    proteome_file=str(proteome_path),
    max_indels=2,
    preprocessed_files_path=str(tmp_path),
    output_format='dataframe'
  ).match()

  actual = set()
  for row in df.iter_rows(named=True):
    if row['Matched Sequence'] is not None:
      actual.add((row['Query Sequence'], row['Protein ID'], row['Matched Sequence']))

  expected = set()
  for _, qseq in queries:
    for pid, pseq in proteins.items():
      for _, matched in brute_force_search(qseq, pseq, 2):
        expected.add((qseq, f'{pid}.1', matched))

  assert actual == expected


def _oracle(qseq, proteins):
  hits = set()
  for pid, pseq in proteins.items():
    for _, matched in brute_force_search(qseq, pseq, 2):
      hits.add((qseq, f'{pid}.1', matched))
  return hits


@pytest.mark.parametrize('seed', range(12))
def test_indel2_short_queries_stay_sound_and_warn(tmp_path, seed, capsys):
  # Exercises the qlen 4-8 regime the batch oracle test deliberately avoids. For 2
  # indels the pigeonhole guarantee needs len >= k*(edits+1) = 6 at the k=2 floor:
  #   - len >= 6  -> recall guaranteed: the engine must find EVERY oracle match.
  #   - len < 6   -> recall NOT guaranteed: we require only SOUNDNESS (no spurious
  #                  hit) and that the engine WARNS; it may miss some oracle matches.
  # This is the regime that let the original bug slip through (old fuzzer used qlen>=8).
  rng = random.Random(1000 + seed)

  proteins = {f'P{i}': _random_sequence(rng, rng.randint(12, 24)) for i in range(3)}
  proteins['REP'] = _repeat_rich_protein(rng)
  proteome_path = tmp_path / 'proteome.fasta'
  proteome_path.write_text(
    ''.join(f'>{pid}\n{seq}\n' for pid, seq in proteins.items())
  )

  # Build queries whose FINAL length lands in 4-8, keeping length >= 3 so the search
  # always has a >=2-mer to seed on.
  queries = []
  for i in range(8):
    protein = proteins[rng.choice(list(proteins))]
    qlen = rng.randint(4, 8)
    start = rng.randrange(0, len(protein) - qlen + 1)
    base = protein[start:start + qlen]
    mode = rng.choice(['exact', 'del1', 'ins1', 'del2', 'ins2'])
    if mode == 'del1' and len(base) >= 4:
      q = _with_deletions(rng, base, 1)
    elif mode == 'del2' and len(base) >= 5:
      q = _with_deletions(rng, base, 2)
    elif mode == 'ins1':
      q = _with_insertions(rng, base, 1)
    elif mode == 'ins2':
      q = _with_insertions(rng, base, 2)
    else:
      q = base
    queries.append((f'q{i}', q))

  query_path = tmp_path / 'queries.fasta'
  query_path.write_text(''.join(f'>{qid}\n{seq}\n' for qid, seq in queries))

  df = Matcher(
    query=str(query_path),
    proteome_file=str(proteome_path),
    max_indels=2,
    preprocessed_files_path=str(tmp_path),
    output_format='dataframe'
  ).match()

  actual = set()
  for row in df.iter_rows(named=True):
    if row['Matched Sequence'] is not None:
      actual.add((row['Query Sequence'], row['Protein ID'], row['Matched Sequence']))

  for _, qseq in queries:
    oracle = _oracle(qseq, proteins)
    reported = {t for t in actual if t[0] == qseq}
    # SOUNDNESS holds in every regime: no reported match may be spurious.
    assert reported <= oracle, f'spurious hits for {qseq!r}: {sorted(reported - oracle)}'
    # COMPLETENESS only where the pigeonhole guarantee applies (len >= 6).
    if len(qseq) >= 6:
      assert reported == oracle, (
        f'guaranteed-regime miss for {qseq!r} (len {len(qseq)}): '
        f'dropped={sorted(oracle - reported)}'
      )

  # If any query is sub-threshold, the recall warning must fire.
  if any(len(qseq) < 6 for _, qseq in queries):
    assert 'complete recall is not guaranteed' in capsys.readouterr().out
