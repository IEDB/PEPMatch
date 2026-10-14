import pytest
import polars.testing as plt

from pathlib import Path

import pepmatch.matcher as matcher_module
from pepmatch import Matcher, Preprocessor
from pepmatch._rs import rs_metadata

DATA = Path(__file__).parent / 'data'

# Three proteins. The second has no gene name (no GN=), so its Gene field is empty and must
# come back null; the third is never hit by the searches below, so its metadata is never read.
PROTEOME_FASTA = (
  '>sp|P11111|AAAA_HUMAN Alpha protein OS=Homo sapiens OX=9606 GN=AAA PE=1 SV=1\n'
  'MKVLATGGWQMKRRDEYH\n'
  '>sp|P22222|BBBB_HUMAN Beta protein OS=Homo sapiens OX=9606 PE=2 SV=3\n'
  'RRHFDSAMFFLLNNMKVLAT\n'
  '>tr|A0A333|A0A333_HUMAN Gamma protein OS=Homo sapiens OX=9606 GN=CCC PE=4 SV=1\n'
  'PQRSTWYVDEGHIKLMNPQ\n'
)


@pytest.fixture
def proteome(tmp_path) -> Path:
  """The three-protein FASTA, preprocessed at k=5 next to it."""
  fasta = tmp_path / 'three.fasta'
  fasta.write_text(PROTEOME_FASTA)
  Preprocessor(str(fasta), preprocessed_files_path=str(tmp_path)).preprocess(5)
  return fasta


def index_of(proteome):
  return str(proteome.with_name('three_5mers.pepidx'))


def rows(columns):
  return list(zip(*columns))


def search(proteome, query):
  return Matcher(
    query=query,
    proteome_file=str(proteome),
    k=5,
    preprocessed_files_path=str(proteome.parent),
    output_format='dataframe',
  ).match()


def test_requested_rows_equal_the_full_tables_rows(proteome):
  full = rows(rs_metadata(index_of(proteome)))
  assert [row[0] for row in full] == [1, 2, 3]
  assert rows(rs_metadata(index_of(proteome), [3, 1])) == [full[0], full[2]]   # ascending


def test_repeated_protein_numbers_collapse(proteome):
  # A repeated row would double every hit joined to it.
  assert [row[0] for row in rows(rs_metadata(index_of(proteome), [2, 2, 2]))] == [2]


def test_no_protein_numbers_return_empty_columns(proteome):
  assert rs_metadata(index_of(proteome), []) == tuple([] for _ in range(10))


@pytest.mark.parametrize('protein', [0, 4])
def test_a_protein_number_outside_the_index_raises(proteome, protein):
  with pytest.raises(ValueError, match=f'protein number {protein} '):
    rs_metadata(index_of(proteome), [protein])


def test_hits_carry_their_own_proteins_metadata(proteome):
  df = search(proteome, ['MKVLAT', 'FFLLNN'])
  got = {
    (row['Query Sequence'], row['Protein ID'], row['Gene'], row['Protein Existence Level'], row['SwissProt Reviewed'])
    for row in df.iter_rows(named=True)
  }
  assert got == {
    ('MKVLAT', 'P11111.1', 'AAA', 1, True),
    ('MKVLAT', 'P22222.3', None, 2, True),
    ('FFLLNN', 'P22222.3', None, 2, True),
  }


METADATA_COLUMNS = [
  'Protein ID', 'Protein Name', 'Species', 'Taxon ID', 'Gene',
  'Protein Existence Level', 'Gene Priority', 'SwissProt Reviewed',
]


def test_a_search_that_hits_nothing_keeps_the_metadata_column_types(proteome):
  # No protein is hit, so no metadata row is read. The metadata columns must still have the
  # types a search with hits gives them, not turn into untyped nulls.
  hit, miss = search(proteome, ['MKVLAT']), search(proteome, ['WWWWWW'])
  assert miss['Matched Sequence'].is_null().all()
  assert {c: miss.schema[c] for c in METADATA_COLUMNS} == {c: hit.schema[c] for c in METADATA_COLUMNS}


SEARCHES = {
  'exact': dict(query=DATA / 'exact_match_query.fasta'),
  'mismatch': dict(query=DATA / 'mismatching_query.fasta', max_mismatches=3),
  'best match': dict(query=DATA / 'best_match_query.fasta', best_match=True),
  'discontinuous with linear': dict(
    query=['L354, V420, G461, Q468, E486, K499, D501, M503, G509', 'T819, L822, A835, F840, S871, Y884, P886',
           'CNERFERLMKELSE', 'FAKEPEPTIDE'],
    max_mismatches=1),
  'one indel': dict(query=DATA / 'indel_query.fasta', max_indels=1),
  'two indels': dict(query=DATA / 'indel_query.fasta', max_indels=2),
  'nothing hit': dict(query=['WWWWWWWWW']),
}


@pytest.mark.parametrize('scenario', SEARCHES)
def test_reading_only_the_hit_proteins_gives_the_frame_reading_every_protein_gives(monkeypatch, tmp_path, scenario):
  # Capture the raw columns each search hands _to_dataframe, then build the frame from those
  # same columns twice: reading only the hit proteins' metadata, as searches now do, and
  # reading every protein's, as they did before. The input is identical, so the comparison
  # is exact even for searches whose hit order changes from run to run.
  captured = []
  to_dataframe = Matcher._to_dataframe

  def capture(self, cols, pepidx_path, is_indels=False):
    captured.append((self, cols, pepidx_path, is_indels))
    return to_dataframe(self, cols, pepidx_path, is_indels)

  monkeypatch.setattr(Matcher, '_to_dataframe', capture)
  Matcher(proteome_file=DATA / 'proteome.fasta', preprocessed_files_path=str(tmp_path), **SEARCHES[scenario]).match()
  monkeypatch.undo()
  assert captured

  hits = 0
  for matcher, cols, pepidx_path, is_indels in captured:
    hit_proteins_only = matcher._to_dataframe(cols, pepidx_path, is_indels)
    monkeypatch.setattr(matcher_module, 'rs_metadata', lambda path, protein_nums=None: rs_metadata(path))
    every_protein = matcher._to_dataframe(cols, pepidx_path, is_indels)
    monkeypatch.undo()
    plt.assert_frame_equal(hit_proteins_only, every_protein)
    hits += hit_proteins_only['Protein ID'].is_not_null().sum()
  assert (hits > 0) == (scenario != 'nothing hit')   # the comparison covered real hits


def test_a_protein_the_index_lacks_raises_with_a_rebuild_hint(proteome):
  # Protein numbers are the index's own FASTA ordinals, so a number this index lacks means the
  # hits came from an index built from a different version of the FASTA (best match without k
  # searches several indexes and reads metadata from the last). Fail loudly and say why.
  matcher = Matcher(query=['MKVLAT'], proteome_file=str(proteome), k=5, preprocessed_files_path=str(proteome.parent))
  with pytest.raises(ValueError, match='protein number 4 .* different versions'):
    matcher._metadata_table(index_of(proteome), [1, 4])
