import random
import pytest

from pepmatch import Matcher, Preprocessor
from indel2_brute_force import brute_force_search

# The index packs each k-mer into 16 bytes. Auto-k is the shortest query's length divided
# by (edits + 1), so before the cap any batch of long enough queries asked for a 17-mer index
# and the preprocessor panicked; best_match started from the shortest query's full length.
PROTEIN = 'GGGGMKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQGGGG'
AMINO_ACIDS = 'ACDEFGHIKLMNPQRSTVWY'


def window(n):
  return PROTEIN[4:4 + n]


@pytest.fixture
def proteome(tmp_path):
  path = tmp_path / 'p.fasta'
  path.write_text(f'>sp|P|P P OS=Test OX=1 GN=P PE=1 SV=1\n{PROTEIN}\n')
  return path


def search(tmp_path, proteome, **kwargs):
  return Matcher(proteome_file=str(proteome), preprocessed_files_path=str(tmp_path),
                 output_format='dataframe', **kwargs).match()


@pytest.mark.parametrize('length, kwargs', [
  (17, {}),                        # exact: auto-k was the full 17
  (30, {}),
  (34, {'max_mismatches': 1}),     # 34 // 2 = 17
  (17, {'best_match': True}),      # best_match started at the shortest query's length
  (34, {'max_indels': 1}),         # 34 // 2 = 17
  (51, {'max_indels': 2}),         # 51 // 3 = 17
])
def test_long_query_is_found_with_k_capped_at_16(tmp_path, proteome, capsys, length, kwargs):
  df = search(tmp_path, proteome, query=[window(length)], **kwargs)
  assert window(length) in df['Matched Sequence'].to_list()
  assert 'k=16' in capsys.readouterr().out


def test_counts_only_long_query(tmp_path, proteome):
  df = search(tmp_path, proteome, query=[window(17)], counts_only=True)
  assert df.filter(df['Mismatches'] == 0)['Count'].to_list() == [1]


def test_k_below_the_cap_is_unchanged(tmp_path, proteome, capsys):
  search(tmp_path, proteome, query=[window(16)])
  assert '(k=16,' in capsys.readouterr().out


def test_explicit_k_above_16_is_rejected(tmp_path, proteome):
  with pytest.raises(ValueError, match='k must be <= 16'):
    search(tmp_path, proteome, query=[window(20)], k=17)


def test_preprocessor_rejects_k_above_16(tmp_path, proteome):
  with pytest.raises(ValueError, match='k must be <= 16'):
    Preprocessor(str(proteome), preprocessed_files_path=str(tmp_path)).preprocess(17)


def _substituted(rng, seq, n):
  for _ in range(n):
    i = rng.randrange(len(seq))
    seq = seq[:i] + rng.choice(AMINO_ACIDS) + seq[i + 1:]
  return seq


@pytest.mark.parametrize('mismatches', [0, 1, 2])
@pytest.mark.parametrize('seed', range(5))
def test_capped_k_keeps_complete_recall_for_mismatches(tmp_path, seed, mismatches):
  # Capping k only adds seeds, so long queries must still match a brute-force Hamming scan.
  rng = random.Random(seed * 3 + mismatches)
  proteins = {f'P{i}': ''.join(rng.choice(AMINO_ACIDS) for _ in range(rng.randint(60, 90))) for i in range(3)}
  queries = []
  for _ in range(6):
    p = proteins[rng.choice(list(proteins))]
    length = rng.randint(17 * (mismatches + 1), min(len(p), 17 * (mismatches + 1) + 8))
    start = rng.randrange(len(p) - length + 1)
    queries.append(_substituted(rng, p[start:start + length], rng.randint(0, mismatches)))
  fasta = tmp_path / 'proteins.fasta'
  fasta.write_text(''.join(f'>{pid}\n{seq}\n' for pid, seq in proteins.items()))

  df = search(tmp_path, fasta, query=queries, max_mismatches=mismatches)
  found = {(r['Query Sequence'], r['Protein ID'], r['Matched Sequence'])
           for r in df.iter_rows(named=True) if r['Matched Sequence'] is not None}
  expected = {
    (q, f'{pid}.1', p[s:s + len(q)])
    for q in queries for pid, p in proteins.items() for s in range(len(p) - len(q) + 1)
    if sum(a != b for a, b in zip(q, p[s:s + len(q)])) <= mismatches
  }
  assert found == expected


@pytest.mark.parametrize('indels', [1, 2])
@pytest.mark.parametrize('seed', range(3))
def test_capped_k_keeps_complete_recall_for_indels(tmp_path, seed, indels):
  rng = random.Random(100 + seed * 2 + indels)
  proteins = {f'P{i}': ''.join(rng.choice(AMINO_ACIDS) for _ in range(rng.randint(70, 90))) for i in range(2)}
  queries = []
  for _ in range(4):
    p = proteins[rng.choice(list(proteins))]
    length = rng.randint(17 * (indels + 1), 17 * (indels + 1) + 6)
    start = rng.randrange(len(p) - length + 1)
    q = p[start:start + length]
    for _ in range(rng.randint(0, indels)):   # interior deletions from the query
      i = rng.randrange(1, len(q) - 1)
      q = q[:i] + q[i + 1:]
    queries.append(q)
  fasta = tmp_path / 'proteins.fasta'
  fasta.write_text(''.join(f'>{pid}\n{seq}\n' for pid, seq in proteins.items()))

  df = search(tmp_path, fasta, query=queries, max_indels=indels)
  found = {(r['Query Sequence'], r['Protein ID'], r['Matched Sequence'])
           for r in df.iter_rows(named=True) if r['Matched Sequence'] is not None}
  expected = {(q, f'{pid}.1', w) for q in queries for pid, p in proteins.items()
              for _, w in brute_force_search(q, p, indels)}
  assert found == expected
