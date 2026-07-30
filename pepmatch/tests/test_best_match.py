import pytest
import polars as pl
from pathlib import Path
from pepmatch import Matcher

@pytest.fixture
def proteome_path() -> Path:
  return Path(__file__).parent / 'data' / 'proteome.fasta'

@pytest.fixture
def query_path() -> Path:
  return Path(__file__).parent / 'data' / 'best_match_query.fasta'

def test_best_match(proteome_path, query_path):
  df = Matcher(
    query=query_path,
    proteome_file=proteome_path,
    best_match=True,
    output_format='dataframe'
  ).match()

  assert df.height == 5
  assert df.filter(pl.col("Matched Sequence").is_not_null()).height >= 3
  assert df.filter(pl.col("Mismatches").is_not_null()).select("Mismatches").to_series().min() >= 0


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

def min_hamming(protein: str, query: str) -> int:
  """Brute-force oracle: fewest substitutions between `query` and any equal-length
  window of `protein`. Independent of the k-mer engine under test."""
  n, m = len(protein), len(query)
  if m > n:
    return m
  return min(
    sum(a != b for a, b in zip(protein[i:i + m], query))
    for i in range(n - m + 1)
  )


def header(acc, name, gene=None, pe=1, sv=1, db='sp'):
  gn = f' GN={gene}' if gene else ''
  return (f'>{db}|{acc}|{acc}_HUMAN {name} OS=Homo sapiens OX=9606'
          f'{gn} PE={pe} SV={sv}')


def run_best_match(tmp_path, fasta_text, queries, **kwargs):
  proteome_path = tmp_path / 'proteome.fasta'
  proteome_path.write_text(fasta_text)
  return Matcher(
    query=queries,
    proteome_file=str(proteome_path),
    best_match=True,
    preprocessed_files_path=str(tmp_path),
    output_format='dataframe',
    **kwargs
  ).match()


# ---------------------------------------------------------------------------
# A. the escalation ladder (best_match_search)
# ---------------------------------------------------------------------------

# One protein, no W anywhere, holding an exact site for EXACT, a 1-off site for
# NEAR and a 4-off site for FAR. Single record so the mismatch pass has no
# neighbouring record to bleed into.
LADDER_PROTEIN = 'ADEFGHIKLNALVEATRFCADEFGHIKLMNPQAAAAYYY'
LADDER_FASTA = header('P11111', 'Alpha protein', gene='AAA') + '\n' + LADDER_PROTEIN + '\n'


def test_best_match_one_row_per_query_across_the_whole_ladder(tmp_path):
  # The ladder matches peptides at different rungs (exact at k=10/mm=0, the 1-off
  # at k=5/mm=1, the 4-off at k=2/mm=4) and never matches the all-W peptide. Whatever
  # rung a peptide leaves on, best_match owes exactly one row per query -- no
  # duplicates from a later rung, no peptide silently dropped.
  queries = ['NALVEATRFC', 'NALVEATRFD', 'KLMNPQRSTV', 'WWWWWWWWWW']
  df = run_best_match(tmp_path, LADDER_FASTA, queries)

  assert df.height == len(queries)
  assert set(df['Query Sequence'].to_list()) == set(queries)
  assert df['Query ID'].n_unique() == len(queries)

  # each reported mismatch count is the brute-force optimum for that peptide
  for row in df.filter(pl.col('Matched Sequence').is_not_null()).iter_rows(named=True):
    assert row['Mismatches'] == min_hamming(LADDER_PROTEIN, row['Query Sequence'])


def test_best_match_unmatchable_peptide_kept_as_null_miss_row(tmp_path):
  # WWWWWWWWWW shares no residue with the proteome, so even the final k=2 rung
  # (budget climbing to len-1) cannot reach it. It must survive as an explicit
  # miss row rather than vanishing from the results.
  df = run_best_match(tmp_path, LADDER_FASTA, ['NALVEATRFC', 'WWWWWWWWWW'])

  miss = df.filter(pl.col('Query Sequence') == 'WWWWWWWWWW')
  assert miss.height == 1
  row = miss.row(0, named=True)
  for column in ('Matched Sequence', 'Protein ID', 'Mismatches', 'Index start', 'Index end'):
    assert row[column] is None, column
  assert row['Mutated Positions'] == '[]'
  # and the matchable peptide beside it is unaffected
  assert df.filter(pl.col('Query Sequence') == 'NALVEATRFC').row(0, named=True)['Mismatches'] == 0


def test_best_match_prefers_the_exact_occurrence_over_a_mismatched_one(tmp_path):
  # NALVEATRFC sits exactly in P22222 and 1-off in P11111. Both are reachable, so
  # the surviving row must be the exact one -- and it must name the protein that
  # actually holds it, not merely report Mismatches == 0.
  fasta = (
    header('P11111', 'Alpha protein', gene='AAA') + '\nDDDDDNALVEATRFDDDDDD\n' +
    header('P22222', 'Beta protein', gene='BBB') + '\nEEEEENALVEATRFCEEEEE\n'
  )
  df = run_best_match(tmp_path, fasta, ['NALVEATRFC'])

  assert df.height == 1
  row = df.row(0, named=True)
  assert row['Mismatches'] == 0
  assert row['Matched Sequence'] == 'NALVEATRFC'
  assert row['Protein ID'] == 'P22222.1'
  assert row['Mutated Positions'] == '[]'


def test_best_match_reports_the_minimum_achievable_mismatch_count(tmp_path):
  # A peptide with no exact site: its two candidate windows are 1 off and 3 off.
  # best_match must escalate only as far as it has to and report the 1, which the
  # brute-force Hamming scan independently confirms is optimal.
  protein = 'GGGGGNALVEATRFCGGGGGNALVEQQQFCGGGGG'
  fasta = header('P11111', 'Alpha protein', gene='AAA') + '\n' + protein + '\n'
  query = 'NALVEATRFD'
  df = run_best_match(tmp_path, fasta, [query])

  assert min_hamming(protein, query) == 1          # oracle sanity: 1 really is best
  assert df.height == 1
  row = df.row(0, named=True)
  assert row['Mismatches'] == 1
  assert row['Matched Sequence'] == 'NALVEATRFC'
  assert row['Mutated Positions'] == '[10]'
  assert protein[row['Index start'] - 1:row['Index end']] == row['Matched Sequence']


def test_best_match_with_explicit_k_still_returns_one_row_per_query(tmp_path):
  # An explicit k bypasses the ladder for a single search + filter. That shortcut
  # must honour the same contract: one row per query, hits and misses alike.
  queries = ['NALVEATRFC', 'NALVEATRFD', 'WWWWWWWWWW']
  df = run_best_match(tmp_path, LADDER_FASTA, queries, k=5, max_mismatches=1)

  assert df.height == len(queries)
  assert set(df['Query Sequence'].to_list()) == set(queries)
  assert df.filter(pl.col('Query Sequence') == 'WWWWWWWWWW').row(0, named=True)['Matched Sequence'] is None


def test_best_match_with_max_indels_raises():
  with pytest.raises(ValueError, match='not yet supported together'):
    Matcher(
      query=['NALVEATRFC'], proteome_file='unused.fasta',
      max_indels=1, best_match=True
    )


def test_best_match_with_counts_only_raises():
  with pytest.raises(ValueError, match='counts_only is not supported with best_match'):
    Matcher(
      query=['NALVEATRFC'], proteome_file='unused.fasta',
      best_match=True, counts_only=True
    )


# ---------------------------------------------------------------------------
# B. the tie-break chain (_best_match_filter)
#
# Each proteome below plants the SAME peptide in several proteins and varies one
# metadata dimension. Where an earlier rung of the chain would otherwise decide,
# a decoy protein carrying the gene first drains Gene Priority to 0 for every
# contender. The favoured protein is also given the WORSE value on a LATER rung
# (and, where possible, placed second in the file) so that dropping the rung
# under test would flip the winner instead of leaving it in place.
# ---------------------------------------------------------------------------

PEPTIDE = 'NALVEATRFC'
GENE_DECOY = header('P00000', 'Decoy protein', gene='AAA') + '\nQQQQQQQQQQQQQQQ\n'


def tie_break_winner(tmp_path, fasta):
  df = run_best_match(tmp_path, fasta, [PEPTIDE])
  assert df.height == 1
  return df.row(0, named=True)['Protein ID']


def test_tie_break_fewest_mismatches_wins(tmp_path):
  # Both proteins are reachable in a single k=5/mm=1 search (explicit k, so the
  # filter -- not the ladder -- is what discards the loser). The 1-off protein is
  # metadata-superior on every later rung, so only the mismatch rung can save the
  # exact hit in the (metadata-inferior, unreviewed, PE=5) protein.
  fasta = (
    header('P11111', 'Alpha protein', gene='AAA', pe=1) + '\nDDDDDNALVEATRFDDDDDD\n' +
    header('Q22222', 'Beta protein Fragment', gene='BBB', pe=5, db='tr') + '\nEEEEENALVEATRFCEEEEE\n'
  )
  proteome_path = tmp_path / 'proteome.fasta'
  proteome_path.write_text(fasta)
  df = Matcher(
    query=[PEPTIDE], proteome_file=str(proteome_path), best_match=True,
    k=5, max_mismatches=1, preprocessed_files_path=str(tmp_path),
    output_format='dataframe'
  ).match()

  assert df.height == 1
  row = df.row(0, named=True)
  assert row['Mismatches'] == 0
  assert row['Protein ID'] == 'Q22222.1'


def test_tie_break_gene_priority_beats_existence_level(tmp_path):
  # Same gene, so the first sp| canonical entry takes Gene Priority 1 and the
  # second gets 0. The priority holder is deliberately the PE=5 one: without the
  # gene-priority rung the PE rung would hand the win to P22222.
  fasta = (
    header('P11111', 'Alpha protein', gene='AAA', pe=5) + '\nWWWWWNALVEATRFCWWWWW\n' +
    header('P22222', 'Beta protein', gene='AAA', pe=1) + '\nYYYYYNALVEATRFCYYYYY\n'
  )
  assert tie_break_winner(tmp_path, fasta) == 'P11111.1'


def test_tie_break_canonical_protein_id_beats_existence_level(tmp_path):
  # A hyphenated (isoform) accession loses to the canonical one. The decoy drains
  # Gene Priority to 0 for both -- otherwise the canonical entry would win one rung
  # earlier and this rung would never be consulted. The canonical entry is the PE=5
  # one, so the later PE rung would pick the isoform if canonicality did not decide.
  fasta = (
    GENE_DECOY +
    header('P11111', 'Alpha protein', gene='AAA', pe=5) + '\nWWWWWNALVEATRFCWWWWW\n' +
    header('P11111-2', 'Alpha isoform', gene='AAA', pe=1) + '\nYYYYYNALVEATRFCYYYYY\n'
  )
  assert tie_break_winner(tmp_path, fasta) == 'P11111.1'


def test_tie_break_swissprot_reviewed_beats_existence_level(tmp_path):
  # sp| beats tr|. Again the decoy flattens Gene Priority, both accessions are
  # canonical, and the reviewed entry carries the worse PE so the PE rung alone
  # would prefer the TrEMBL one.
  fasta = (
    GENE_DECOY +
    header('P11111', 'Alpha protein', gene='AAA', pe=5) + '\nWWWWWNALVEATRFCWWWWW\n' +
    header('Q22222', 'Beta protein', gene='AAA', pe=1, db='tr') + '\nYYYYYNALVEATRFCYYYYY\n'
  )
  assert tie_break_winner(tmp_path, fasta) == 'P11111.1'


def test_tie_break_lower_existence_level_wins(tmp_path):
  # Distinct genes keep Gene Priority at 1 for both; both are canonical, reviewed
  # and non-fragment. PE=1 vs PE=5 is the only difference, and the winner is placed
  # second so a bare keep-first would pick the loser.
  fasta = (
    header('P11111', 'Alpha protein', gene='AAA', pe=5) + '\nWWWWWNALVEATRFCWWWWW\n' +
    header('P22222', 'Beta protein', gene='BBB', pe=1) + '\nYYYYYNALVEATRFCYYYYY\n'
  )
  assert tie_break_winner(tmp_path, fasta) == 'P22222.1'


def test_tie_break_non_fragment_protein_name_wins(tmp_path):
  # Everything ties through the PE rung; the only difference is the word Fragment
  # in the protein name. The fragment is listed first, so keep-first would choose
  # it and only the fragment rung can prevent that.
  fasta = (
    header('P11111', 'Alpha protein (Fragment)', gene='AAA', pe=1) + '\nWWWWWNALVEATRFCWWWWW\n' +
    header('P22222', 'Beta protein', gene='BBB', pe=1) + '\nYYYYYNALVEATRFCYYYYY\n'
  )
  assert tie_break_winner(tmp_path, fasta) == 'P22222.1'


def test_tie_break_keeps_one_row_when_every_dimension_ties(tmp_path):
  # Two indistinguishable proteins (and a repeat of the peptide inside one of them,
  # which the engine reports as its own row): the chain runs out of discriminators,
  # so the final unique() must still collapse everything to a single row.
  fasta = (
    header('P11111', 'Alpha protein', gene='AAA', pe=1) + '\nWWWWWNALVEATRFCWWWWWNALVEATRFCWWWWW\n' +
    header('P22222', 'Beta protein', gene='BBB', pe=1) + '\nYYYYYNALVEATRFCYYYYY\n'
  )
  df = run_best_match(tmp_path, fasta, [PEPTIDE])
  assert df.height == 1
  assert df.row(0, named=True)['Mismatches'] == 0
