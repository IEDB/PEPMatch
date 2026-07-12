import pytest
import polars as pl
import polars.testing as plt
from pathlib import Path
from pepmatch import Matcher
from pepmatch.matcher import format_indel_positions

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


def test_max_indels_greater_than_two_raises():
  with pytest.raises(ValueError, match='max_indels > 2 is not yet supported'):
    Matcher(query=['NALVEATRFC'], proteome_file='unused.fasta', max_indels=3)


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


def test_max_indels_and_discontinuous_raises():
  # A discontinuous (position-anchored) epitope has no contiguous window for the
  # indel engine to seed/extend, so the combination is undefined — it must error
  # up front rather than crash later when the indel and discontinuous result
  # frames (now carrying different edit-count columns) fail to concat.
  with pytest.raises(ValueError, match='discontinuous'):
    Matcher(query=['N1, A3, L5'], proteome_file='unused.fasta', max_indels=1)


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


def test_indel_mode_emits_indels_column_only(proteome_path):
  # One edit-count and one edit-position column per mode: an indel search reports
  # Indels + Indel Positions and must NOT carry the mismatch-mode twins.
  df = Matcher(
    query=['NALVEATRFC'],
    proteome_file=proteome_path,
    max_indels=1,
    output_format='dataframe'
  ).match()
  assert 'Indels' in df.columns and 'Indel Positions' in df.columns
  assert 'Mismatches' not in df.columns and 'Mutated Positions' not in df.columns


def test_mismatch_mode_emits_mismatches_column_only(proteome_path):
  # The mirror: a non-indel search keeps Mismatches + Mutated Positions and must
  # NOT gain the indel-mode twins suggesting an indel search that never ran.
  df = Matcher(
    query=['NALVEATRFC'],
    proteome_file=proteome_path,
    max_mismatches=0,
    output_format='dataframe'
  ).match()
  assert 'Mismatches' in df.columns and 'Mutated Positions' in df.columns
  assert 'Indels' not in df.columns and 'Indel Positions' not in df.columns


def test_indel_positions_annotation_unit():
  # Hand-verified annotations, format `<d|i>: <residues>[<pos or range>]`, 1-based.
  # Deletion residue comes from the query, insertion residue from the protein. In a
  # repeat the exact position is ambiguous, so a range of all valid positions is
  # reported (collapsing to a single number when unambiguous).
  assert format_indel_positions('ABCDEF', 'ABCDEF') == '[]'           # exact
  assert format_indel_positions('YYADGY', 'YADGY') == 'd: Y[2]'       # only pos 2 (1 is terminal)
  assert format_indel_positions('NALVEATRFC', 'NALVETRFC') == 'd: A[6]'  # the 2nd A, unambiguous
  assert format_indel_positions('ABCDEF', 'ABXCDEF') == 'i: X[3]'     # X inserted before C
  assert format_indel_positions('AAAAA', 'AAAA') == 'd: A[2,4]'       # deletable at 2, 3 or 4
  assert format_indel_positions('AAAAAA', 'AAAAA') == 'd: A[2,5]'
  assert format_indel_positions('AAAAAA', 'AAAAAAA') == 'i: A[2,6]'   # insertable across the run


def test_indel_positions_deletion_end_to_end(tmp_path):
  # Query NALVEATRFC vs a protein missing the 2nd A -> matched NALVETRFC,
  # annotated as a deletion of A at query position 6.
  proteome_path = tmp_path / 'proteome.fasta'
  proteome_path.write_text('>P\nMKVNALVETRFCGHI\n')
  df = Matcher(
    query=['NALVEATRFC'],
    proteome_file=str(proteome_path),
    max_indels=1,
    preprocessed_files_path=str(tmp_path),
    output_format='dataframe'
  ).match()
  row = df.filter(pl.col('Matched Sequence') == 'NALVETRFC')
  assert row.height == 1
  assert row['Indel Positions'].item() == 'd: A[6]'


def test_two_consecutive_deletions_found(tmp_path):
  # YYADGY loses the contiguous YA -> YDGY: one 2-residue deletion event.
  proteome_path = tmp_path / 'proteome.fasta'
  proteome_path.write_text('>P\nMKYDGYWW\n')
  df = Matcher(
    query=['YYADGY'], proteome_file=str(proteome_path), max_indels=2,
    preprocessed_files_path=str(tmp_path), output_format='dataframe'
  ).match()
  row = df.filter(pl.col('Matched Sequence') == 'YDGY')
  assert row.height == 1
  assert row['Indels'].item() == 2
  assert row['Indel Positions'].item() == 'd: YA[2]'


def test_two_non_consecutive_deletions_found(tmp_path):
  # ABCDEF loses C (3) and E (5) -> ABDF: two independent deletion events, so the
  # annotation reports them separately rather than as one chunk.
  proteome_path = tmp_path / 'proteome.fasta'
  proteome_path.write_text('>P\nMKABDFWW\n')
  df = Matcher(
    query=['ABCDEF'], proteome_file=str(proteome_path), max_indels=2,
    preprocessed_files_path=str(tmp_path), output_format='dataframe'
  ).match()
  row = df.filter(pl.col('Matched Sequence') == 'ABDF')
  assert row.height == 1
  assert row['Indels'].item() == 2
  assert row['Indel Positions'].item() == 'd: C[3], d: E[5]'


def test_two_insertions_found(tmp_path):
  proteome_path = tmp_path / 'proteome.fasta'
  proteome_path.write_text('>P\nMKABXCDYEFWW\n')
  df = Matcher(
    query=['ABCDEF'], proteome_file=str(proteome_path), max_indels=2,
    preprocessed_files_path=str(tmp_path), output_format='dataframe'
  ).match()
  row = df.filter(pl.col('Matched Sequence') == 'ABXCDYEF')
  assert row.height == 1
  assert row['Indels'].item() == 2
  assert row['Indel Positions'].item() == 'i: X[3], i: Y[5]'


def test_mixed_insertion_and_deletion_not_found(tmp_path):
  # max_indels=2 is HOMOGENEOUS: two insertions or two deletions, never one of each.
  # ABXCDF would need an inserted X plus a deleted E — that pair is a substitution,
  # which mismatch search covers, so the indel engine must not report it.
  proteome_path = tmp_path / 'proteome.fasta'
  proteome_path.write_text('>P\nMKABXCDFWW\n')
  df = Matcher(
    query=['ABCDEF'], proteome_file=str(proteome_path), max_indels=2,
    preprocessed_files_path=str(tmp_path), output_format='dataframe'
  ).match()
  matched = [m for m in df['Matched Sequence'].to_list() if m]
  assert 'ABXCDF' not in matched


def test_indel2_positions_annotation_unit():
  # Two edits at one site collapse to a chunk; otherwise they report separately. In a
  # repeat the position is a range — but only where the removed residues stay the same.
  assert format_indel_positions('YYADGY', 'YDGY') == 'd: YA[2]'        # consecutive chunk
  assert format_indel_positions('AAAAAA', 'AAAA') == 'd: AA[2,4]'      # homopolymer chunk
  assert format_indel_positions('AAAA', 'AAAAAA') == 'i: AA[2,4]'      # insertion chunk
  assert format_indel_positions('ABCDEF', 'ABDF') == 'd: C[3], d: E[5]'
  assert format_indel_positions('AAABAA', 'AABA') == 'd: A[2,3], d: A[5]'
  # Periodic repeat: the removed pair changes along the range (BA/AB/BA), so the range
  # must NOT be collapsed or it would misreport which residues went missing.
  assert format_indel_positions('ABABAB', 'ABAB') == 'd: BA[2], d: AB[3], d: BA[4]'


def test_indel_positions_insertion_end_to_end(tmp_path):
  # Query NALVEATRFC vs a protein with an extra X after E -> matched NALVEXATRFC,
  # annotated as an insertion of X before query position 6.
  proteome_path = tmp_path / 'proteome.fasta'
  proteome_path.write_text('>P\nMKVNALVEXATRFCGHI\n')
  df = Matcher(
    query=['NALVEATRFC'],
    proteome_file=str(proteome_path),
    max_indels=1,
    preprocessed_files_path=str(tmp_path),
    output_format='dataframe'
  ).match()
  row = df.filter(pl.col('Matched Sequence') == 'NALVEXATRFC')
  assert row.height == 1
  assert row['Indel Positions'].item() == 'i: X[6]'


def test_indel_honors_explicit_k_up_to_optimal_and_reuses_table(tmp_path, capsys):
  # Issue #28: a user who already built a k-mer table (e.g. for a mismatch run) and
  # passes that same k to indel search should have it reused, not silently swapped
  # for the optimized k and rebuilt. Query len 10, max_indels=1 -> optimized k = 5;
  # k=3 <= 5 so it is honored, and the pre-built 3-mer index is used as-is.
  proteome_path = tmp_path / 'proteome.fasta'
  proteome_path.write_text('>ProtDel\nMKVNALVETRFCGHI\n')
  # the user "already has" a 3-mer table from an earlier run
  Matcher(query=['NALVEATRFC'], proteome_file=str(proteome_path), max_mismatches=0,
          k=3, preprocessed_files_path=str(tmp_path), output_format='dataframe').match()
  assert (tmp_path / 'proteome_3mers.pepidx').exists()
  capsys.readouterr()  # discard the pre-build output
  df = Matcher(query=['NALVEATRFC'], proteome_file=str(proteome_path), max_indels=1,
               k=3, preprocessed_files_path=str(tmp_path), output_format='dataframe').match()
  out = capsys.readouterr().out
  assert '(k=3,' in out                       # honored the user's k, not the optimal 5
  assert 'No preprocessed file' not in out    # reused the existing table, no rebuild
  assert 'NALVETRFC' in df['Matched Sequence'].to_list()  # recall preserved at k=3


def test_indel_clamps_k_above_optimal_and_warns(tmp_path, capsys):
  # An explicit k above the pigeonhole optimal drops below max_indels+1 disjoint
  # seeds and would forfeit complete recall, so indel search clamps it to the
  # optimal and warns. Query len 10, max_indels=1 -> optimized k = 5; k=7 (still
  # <= min_len, so it passes the constructor guard) is clamped down to 5.
  proteome_path = tmp_path / 'proteome.fasta'
  proteome_path.write_text('>ProtDel\nMKVNALVETRFCGHI\n')
  df = Matcher(query=['NALVEATRFC'], proteome_file=str(proteome_path), max_indels=1,
               k=7, preprocessed_files_path=str(tmp_path), output_format='dataframe').match()
  out = capsys.readouterr().out
  assert 'Requested k=7 exceeds k=5' in out    # warned about the clamp
  assert '(k=5,' in out                        # searched at the clamped optimal
  assert 'NALVETRFC' in df['Matched Sequence'].to_list()   # recall still complete
  assert (tmp_path / 'proteome_5mers.pepidx').exists()     # built the optimal index
  assert not (tmp_path / 'proteome_7mers.pepidx').exists() # never the too-large one


def test_indel_unspecified_k_derives_optimal(tmp_path, capsys):
  # With no k passed, indel search derives the pigeonhole optimal: max(2, 10//2) = 5.
  proteome_path = tmp_path / 'proteome.fasta'
  proteome_path.write_text('>ProtDel\nMKVNALVETRFCGHI\n')
  Matcher(query=['NALVEATRFC'], proteome_file=str(proteome_path), max_indels=1,
          preprocessed_files_path=str(tmp_path), output_format='dataframe').match()
  out = capsys.readouterr().out
  assert '(k=5,' in out
  assert (tmp_path / 'proteome_5mers.pepidx').exists()


def test_indel_k_exceeding_shortest_query_raises():
  # The constructor rejects an explicit k larger than the shortest query peptide
  # (shorter peptides would be silently skipped). The guard is global, so it also
  # governs indel search -- an explicit k reaching indel_search is always <= min_len.
  with pytest.raises(ValueError, match='cannot exceed the shortest query peptide'):
    Matcher(query=['NALVEATRFC'], proteome_file='unused.fasta', max_indels=1, k=11)
