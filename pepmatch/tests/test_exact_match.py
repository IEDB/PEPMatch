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
  return Path(__file__).parent / 'data' / 'exact_match_query.fasta'

@pytest.fixture
def expected_path() -> Path:
  return Path(__file__).parent / 'data' / 'exact_match_expected.csv'

def test_exact_match_k5(proteome_path, query_path, expected_path):
  df = Matcher(
    query=query_path,
    proteome_file=proteome_path,
    max_mismatches=0,
    k=5,
    output_format='dataframe'
  ).match()

  df = df.sort('Query Sequence')
  expected_df = pl.read_csv(expected_path).sort('Query Sequence')
  plt.assert_series_equal(df['Protein ID'], expected_df['Protein ID'])

def test_exact_match_k9(proteome_path, query_path, expected_path):
  df = Matcher(
    query=query_path,
    proteome_file=proteome_path,
    max_mismatches=0,
    k=9,
    output_format='dataframe'
  ).match()

  df = df.sort('Query Sequence')
  expected_df = pl.read_csv(expected_path).sort('Query Sequence')
  plt.assert_series_equal(df['Protein ID'], expected_df['Protein ID'])


# A tiny hand-built proteome. Every reported Protein ID carries its sequence version:
# SV=2 for P11111, SV=1 for P22222, and Q33333 has no SV= at all so it defaults to 1.
MINI_PROTEINS = {
  'P11111.2': 'MKVLATGGWQMKVLATRRDE',
  'P22222.1': 'YHFDSAMKVLATPCNQEI',
  'Q33333.1': 'WWCCTTGGAAFFLLNNDD',
}

MINI_FASTA = (
  '>sp|P11111|AAA_HUMAN Alpha protein OS=Homo sapiens OX=9606 GN=AAA PE=1 SV=2\n'
  'MKVLATGGWQMKVLATRRDE\n'
  '>sp|P22222|BBB_HUMAN Beta protein OS=Homo sapiens OX=9606 GN=BBB PE=2 SV=1\n'
  'YHFDSAMKVLATPCNQEI\n'
  '>tr|Q33333|CCC_HUMAN Gamma protein OS=Homo sapiens OX=9606 GN=CCC PE=4\n'
  'WWCCTTGGAAFFLLNNDD\n'
)

# MKVLAT hits twice in P11111 and once in P22222; GGWQMK/HFDSAM/FFLLNN hit once each;
# QQQQQQ is absent; RRDEYH straddles the P11111/P22222 record boundary. All length 6,
# so every k in 2..6 is a valid explicit k for this batch.
MINI_QUERIES = ['MKVLAT', 'GGWQMK', 'HFDSAM', 'FFLLNN', 'QQQQQQ', 'RRDEYH']


@pytest.fixture
def mini_proteome(tmp_path) -> Path:
  path = tmp_path / 'mini.fasta'
  path.write_text(MINI_FASTA)
  return path


def brute_force_hits(queries):
  """Oracle: every (query, protein id, matched, 1-based start) by plain substring scan."""
  hits = set()
  for query in queries:
    for protein_id, seq in MINI_PROTEINS.items():
      for i in range(len(seq) - len(query) + 1):
        if seq[i:i + len(query)] == query:
          hits.add((query, protein_id, query, i + 1))
  return hits


def reported_hits(df):
  """The same tuples, read off a result frame with the miss rows dropped."""
  rows = df.filter(pl.col('Matched Sequence').is_not_null())
  return set(zip(
    rows['Query Sequence'], rows['Protein ID'], rows['Matched Sequence'], rows['Index start']
  ))


def run_mini(mini_proteome, tmp_path, **kwargs):
  kwargs.setdefault('query', MINI_QUERIES)
  return Matcher(
    proteome_file=str(mini_proteome),
    max_mismatches=0,
    preprocessed_files_path=str(tmp_path),
    output_format='dataframe',
    **kwargs
  ).match()


def test_exact_match_equals_bruteforce_oracle(mini_proteome, tmp_path):
  # Exact search is definitionally a substring scan, so the engine must agree with one
  # residue for residue -- same proteins, same 1-based offsets, no extras, no drops.
  assert len(brute_force_hits(['MKVLAT'])) == 3  # the repeat is really in the fixture
  assert brute_force_hits(['QQQQQQ']) == set()   # and so is a genuine miss
  df = run_mini(mini_proteome, tmp_path)
  assert reported_hits(df) == brute_force_hits(MINI_QUERIES)


def test_exact_match_never_spans_a_protein_boundary(mini_proteome, tmp_path):
  # RRDEYH is P11111's last four residues followed by P22222's first two. The records
  # sit back to back in the index; a match here would mean the search walked off the
  # end of one protein into the next.
  df = run_mini(mini_proteome, tmp_path, query=['RRDEYH'])
  assert df.height == 1
  assert df['Matched Sequence'].item() is None


def test_exact_match_coordinates_locate_the_matched_sequence(mini_proteome, tmp_path):
  # Index start/end must be a usable 1-based inclusive slice of the protein, and an
  # exact hit is by definition the query itself with zero mismatches.
  df = run_mini(mini_proteome, tmp_path).filter(pl.col('Matched Sequence').is_not_null())
  assert df.height == len(brute_force_hits(MINI_QUERIES))
  for row in df.iter_rows(named=True):
    protein = MINI_PROTEINS[row['Protein ID']]
    start, end = row['Index start'], row['Index end']
    assert protein[start - 1:end] == row['Matched Sequence']
    assert end - start + 1 == len(row['Query Sequence'])
    assert row['Matched Sequence'] == row['Query Sequence']
    assert row['Mismatches'] == 0
    assert row['Mutated Positions'] == '[]'


def test_exact_match_miss_row_shape(mini_proteome, tmp_path):
  # A miss is reported, not omitted: one row, the query echoed back, every match and
  # metadata field null -- except Mutated Positions, which is the empty-list string.
  df = run_mini(mini_proteome, tmp_path, query=['QQQQQQ'])
  assert df.height == 1
  row = df.row(0, named=True)
  for column in (
    'Matched Sequence', 'Protein ID', 'Protein Name', 'Species', 'Taxon ID', 'Gene',
    'Mismatches', 'Index start', 'Index end', 'Protein Existence Level',
    'Gene Priority', 'SwissProt Reviewed',
  ):
    assert row[column] is None, column
  assert row['Query Sequence'] == 'QQQQQQ'
  assert row['Mutated Positions'] == '[]'


def test_every_query_represented_in_mixed_batch(mini_proteome, tmp_path):
  # No query may vanish from a batch of hits and misses, and the row count per query
  # is exactly its hit count -- or one, for a miss.
  df = run_mini(mini_proteome, tmp_path)
  assert set(df['Query Sequence']) == set(MINI_QUERIES)
  assert set(df['Query ID']) == {str(i + 1) for i in range(len(MINI_QUERIES))}
  counts = {query: 0 for query in MINI_QUERIES}
  for query in df['Query Sequence']:
    counts[query] += 1
  for query in MINI_QUERIES:
    assert counts[query] == max(1, len(brute_force_hits([query]))), query


def test_exact_match_invariant_to_k(mini_proteome, tmp_path):
  # k sizes the seed k-mer, nothing else. Every k valid for these length-6 queries --
  # and the auto-derived one -- must return an identical match set.
  expected = brute_force_hits(MINI_QUERIES)
  for k in (2, 3, 5, 6):
    assert reported_hits(run_mini(mini_proteome, tmp_path, k=k)) == expected, k
  assert reported_hits(run_mini(mini_proteome, tmp_path)) == expected


def test_query_input_forms_agree(mini_proteome, tmp_path):
  # A list, a one-per-line .txt and a FASTA are three spellings of one query set;
  # only the Query IDs may differ, never the matches.
  txt_path = tmp_path / 'queries.txt'
  txt_path.write_text('\n'.join(MINI_QUERIES) + '\n')
  fasta_path = tmp_path / 'queries.fasta'
  fasta_path.write_text(''.join(f'>q{i}\n{seq}\n' for i, seq in enumerate(MINI_QUERIES)))

  expected = brute_force_hits(MINI_QUERIES)
  assert reported_hits(run_mini(mini_proteome, tmp_path)) == expected
  assert reported_hits(run_mini(mini_proteome, tmp_path, query=str(txt_path))) == expected
  assert reported_hits(run_mini(mini_proteome, tmp_path, query=str(fasta_path))) == expected


def test_lowercase_query_is_upcased_before_searching(mini_proteome, tmp_path):
  # The proteome is upper case, so a lower-case query would miss everything unless it
  # is normalized on parse -- and the normalized form is what gets reported back.
  df = run_mini(mini_proteome, tmp_path, query=[query.lower() for query in MINI_QUERIES])
  assert reported_hits(df) == brute_force_hits(MINI_QUERIES)
  assert set(df['Query Sequence']) == set(MINI_QUERIES)


def test_k_below_two_raises(mini_proteome):
  # k=1 degenerates the index to single residues; rejected up front.
  with pytest.raises(ValueError, match='k must be >= 2'):
    Matcher(query=MINI_QUERIES, proteome_file=str(mini_proteome), k=1)


def test_k_above_shortest_query_raises(mini_proteome):
  # A peptide shorter than k has no k-mer seed and would be silently skipped, so the
  # constructor refuses. k == the shortest peptide is still fine.
  with pytest.raises(ValueError, match='cannot exceed the shortest query peptide'):
    Matcher(query=MINI_QUERIES, proteome_file=str(mini_proteome), k=7)
  Matcher(query=MINI_QUERIES, proteome_file=str(mini_proteome), k=6)


def test_missing_query_file_raises(mini_proteome, tmp_path):
  with pytest.raises(FileNotFoundError):
    Matcher(query=str(tmp_path / 'absent.txt'), proteome_file=str(mini_proteome))


def test_unsupported_query_extension_raises(mini_proteome, tmp_path):
  # The file exists, so this is the extension check firing and not FileNotFoundError.
  bad_path = tmp_path / 'queries.docx'
  bad_path.write_text('MKVLAT\n')
  with pytest.raises(ValueError, match='Unsupported query format'):
    Matcher(query=str(bad_path), proteome_file=str(mini_proteome))


def test_sequence_version_toggle_controls_accession_suffix(mini_proteome, tmp_path):
  # Default appends the header's SV (defaulting to 1 when absent); the flag off gives
  # the bare accession, which is what joins against version-less external tables.
  queries = ['MKVLAT', 'FFLLNN']
  versioned = run_mini(mini_proteome, tmp_path, query=queries)
  bare = run_mini(mini_proteome, tmp_path, query=queries, sequence_version=False)
  assert set(versioned['Protein ID']) == {'P11111.2', 'P22222.1', 'Q33333.1'}
  assert set(bare['Protein ID']) == {'P11111', 'P22222', 'Q33333'}
