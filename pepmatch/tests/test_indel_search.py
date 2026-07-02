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
  return Path(__file__).parent / 'data' / 'indel_query.fasta'

@pytest.fixture
def expected_path() -> Path:
  return Path(__file__).parent / 'data' / 'indel_expected.csv'

def test_indel_search(proteome_path, query_path, expected_path):
  df = Matcher(
    query=query_path,
    proteome_file=proteome_path,
    max_indels=1,
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


def test_terminal_deletion_allowed_with_single_residue_buffer(tmp_path):
  # Regression test: p_idx==0/protein_len-1 are real residues, not out-of-bounds —
  # 1 residue of protein buffer is enough, not 2. Uses interior query positions
  # (d=1,3 of 'ABCDE') to isolate this from the query-terminal guard below.
  proteome_path = tmp_path / 'proteome.fasta'
  proteome_path.write_text('>Front\nACDEFG\n>Back\nXYABCE\n')
  df = Matcher(
    query=['ABCDE'],
    proteome_file=str(proteome_path),
    max_indels=1,
    preprocessed_files_path=str(tmp_path),
    output_format='dataframe'
  ).match()
  hits = set(zip(df['Protein ID'].to_list(), df['Matched Sequence'].to_list()))
  assert ('Front.1', 'ACDE') in hits
  assert ('Back.1', 'ABCE') in hits


def test_query_terminal_deletion_blocked(proteome_path):
  # Deleting Q at query position 0 is the only way QNALVEATRFC could match
  # NALVEATRFC in Q8V336 — blocked unconditionally, a query-boundary deletion.
  df = Matcher(
    query=['QNALVEATRFC'],
    proteome_file=proteome_path,
    max_indels=1,
    output_format='dataframe'
  ).match()
  matches = df.filter(pl.col('Matched Sequence').is_not_null())
  assert 'NALVEATRFC' not in matches['Matched Sequence'].to_list(), (
    'Query-terminal deletion match NALVEATRFC should be blocked'
  )


def test_terminal_insertion_blocked(tmp_path):
  # A leading/trailing insertion pads the query with no context to confirm
  # it's genuine, not an arbitrary flanking residue — the exact match still hits.
  proteome_path = tmp_path / 'proteome.fasta'
  proteome_path.write_text(
    '>Leading\nXNALVEATRFC\n'
    '>Trailing\nNALVEATRFCX\n'
  )
  df = Matcher(
    query=['NALVEATRFC'],
    proteome_file=str(proteome_path),
    max_indels=1,
    preprocessed_files_path=str(tmp_path),
    output_format='dataframe'
  ).match()
  hits = set(zip(df['Protein ID'].to_list(), df['Matched Sequence'].to_list()))
  assert ('Leading.1', 'NALVEATRFC') in hits
  assert ('Trailing.1', 'NALVEATRFC') in hits
  assert ('Leading.1', 'XNALVEATRFC') not in hits
  assert ('Trailing.1', 'NALVEATRFCX') not in hits


def test_insertion_at_second_to_last_position_found(tmp_path):
  # The repeated S makes the inserted residue's exact position ambiguous, but
  # the second-to-last placement is verifiable (query's final S still matches
  # a real residue) and must be found — unlike a genuinely terminal insertion.
  proteome_path = tmp_path / 'proteome.fasta'
  proteome_path.write_text('>P\nLPDGVWEESS\n')
  df = Matcher(
    query=['LPDGVWEES'],
    proteome_file=str(proteome_path),
    max_indels=1,
    preprocessed_files_path=str(tmp_path),
    output_format='dataframe'
  ).match()
  hits = set(zip(df['Protein ID'].to_list(), df['Matched Sequence'].to_list()))
  assert ('P.1', 'LPDGVWEES') in hits
  assert ('P.1', 'LPDGVWEESS') in hits


def test_max_indels_greater_than_one_raises():
  with pytest.raises(ValueError, match='max_indels > 1 is not yet supported'):
    Matcher(query=['NALVEATRFC'], proteome_file='unused.fasta', max_indels=2)


def test_indels_and_mismatches_mutually_exclusive_raises():
  with pytest.raises(ValueError, match='mutually exclusive'):
    Matcher(
      query=['NALVEATRFC'], proteome_file='unused.fasta',
      max_indels=1, max_mismatches=1
    )


def test_max_indels_and_best_match_raises():
  with pytest.raises(ValueError, match='not yet supported together'):
    Matcher(
      query=['NALVEATRFC'], proteome_file='unused.fasta',
      max_indels=1, best_match=True
    )


def test_max_indels_and_counts_only_raises():
  with pytest.raises(ValueError, match='not yet supported together'):
    Matcher(
      query=['NALVEATRFC'], proteome_file='unused.fasta',
      max_indels=1, counts_only=True
    )


def test_indel_peptide_shorter_than_k(proteome_path):
  # A single-residue query forces k = max(2, min_len // 2) = 2 while
  # peptide_len=1 < k; must miss cleanly rather than erroring.
  df = Matcher(
    query=['A'],
    proteome_file=proteome_path,
    max_indels=1,
    output_format='dataframe'
  ).match()
  assert df['Matched Sequence'].is_null().all()


def test_indel_multi_hit_different_proteins(tmp_path):
  # NALVEATRFC (10 residues) with 1 indel matches two DIFFERENT proteins:
  # ProtDel is missing the second A (a deletion), ProtIns has an extra X
  # inserted after E (an insertion). Both hand-verified against the DFS's
  # seed ("NALVE") + bidirectional extension.
  proteome_path = tmp_path / 'proteome.fasta'
  proteome_path.write_text(
    '>ProtDel\nMKVNALVETRFCGHI\n'
    '>ProtIns\nMKVNALVEXATRFCGHI\n'
  )
  df = Matcher(
    query=['NALVEATRFC'],
    proteome_file=str(proteome_path),
    max_indels=1,
    preprocessed_files_path=str(tmp_path),
    output_format='dataframe'
  ).match()
  hits = set(zip(df['Protein ID'].to_list(), df['Matched Sequence'].to_list()))
  assert ('ProtDel.1', 'NALVETRFC') in hits
  assert ('ProtIns.1', 'NALVEXATRFC') in hits
