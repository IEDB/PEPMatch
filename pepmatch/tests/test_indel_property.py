import random
import pytest
from pepmatch import Matcher
from indel import brute_force_search

AMINO_ACIDS = 'ACDEFGHIKLMNPQRSTVWY'


def _random_sequence(rng, length):
  return ''.join(rng.choice(AMINO_ACIDS) for _ in range(length))


# Named for the edit applied to the query string, not the resulting Indels label —
# e.g. _with_deletion removes a query residue, so the protein has "extra" content
# there relative to the query, which the search reports as an insertion match.
def _with_deletion(rng, seq):
  d = rng.randrange(len(seq))
  return seq[:d] + seq[d + 1:]


def _with_insertion(rng, seq):
  i = rng.randrange(1, len(seq))
  return seq[:i] + rng.choice(AMINO_ACIDS) + seq[i:]


def _random_query(rng, proteins):
  protein_id = rng.choice(list(proteins))
  protein = proteins[protein_id]
  qlen = rng.randint(8, 14)
  start = rng.randrange(0, len(protein) - qlen + 1)
  base = protein[start:start + qlen]
  mode = rng.choice(['exact', 'deletion', 'insertion', 'random'])
  if mode == 'deletion':
    return _with_deletion(rng, base)
  if mode == 'insertion':
    return _with_insertion(rng, base)
  if mode == 'random':
    return _random_sequence(rng, qlen)
  return base


@pytest.mark.parametrize('seed', range(25))
def test_indel_search_matches_brute_force_oracle(tmp_path, seed):
  rng = random.Random(seed)

  proteins = {f'P{i}': _random_sequence(rng, rng.randint(20, 40)) for i in range(4)}
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
    max_indels=1,
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
      for _, matched in brute_force_search(qseq, pseq):
        expected.add((qseq, f'{pid}.1', matched))

  assert actual == expected
