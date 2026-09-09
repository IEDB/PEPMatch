import pytest

from pathlib import Path

from pepmatch import Matcher, Preprocessor

# Two references whose stems share a prefix: 'ref' and 'ref_extra'. Their sequences are
# identical, so an index built from either resolves the same hits at the same protein
# numbers -- only the metadata differs. That is exactly the case where reading metadata
# from the wrong index is silent: same rows, same coordinates, wrong protein identity.
REFERENCE_FASTA = (
  '>sp|P11111|AAAA_HUMAN Alpha protein OS=Homo sapiens OX=9606 GN=AAA PE=1 SV=1\n'
  'MKVLATGGWQMKRRDEYH\n'
  '>sp|P22222|BBBB_HUMAN Beta protein OS=Homo sapiens OX=9606 GN=BBB PE=1 SV=1\n'
  'RRHFDSAMFFLLNNMKVLAT\n'
)

DECOY_FASTA = (
  '>sp|Q99999|ZZZZ_MOUSE Decoy alpha OS=Mus musculus OX=10090 GN=ZZZ PE=1 SV=1\n'
  'MKVLATGGWQMKRRDEYH\n'
  '>sp|Q88888|YYYY_MOUSE Decoy beta OS=Mus musculus OX=10090 GN=YYY PE=1 SV=1\n'
  'RRHFDSAMFFLLNNMKVLAT\n'
)

REFERENCE_IDS = {'P11111.1', 'P22222.1'}
DECOY_IDS = {'Q99999.1', 'Q88888.1'}


@pytest.fixture
def references(tmp_path) -> Path:
  """Both FASTAs preprocessed into one directory, every k the tests could touch. The
  decoy indexes are written last so a directory listing is as likely to hand them back
  first as the reference's own -- which is precisely how a prefix-matching lookup picks
  the wrong reference."""
  reference_path = tmp_path / 'ref.fasta'
  reference_path.write_text(REFERENCE_FASTA)
  decoy_path = tmp_path / 'ref_extra.fasta'
  decoy_path.write_text(DECOY_FASTA)
  for fasta in (reference_path, decoy_path):
    for k in range(2, 7):
      Preprocessor(str(fasta), preprocessed_files_path=str(tmp_path)).preprocess(k)
  return reference_path


def hits(df):
  """Reported (query, protein id, matched sequence) tuples, miss rows dropped."""
  return {
    (row['Query Sequence'], row['Protein ID'], row['Matched Sequence'])
    for row in df.iter_rows(named=True)
    if row['Matched Sequence'] is not None
  }


def test_specified_k_reads_metadata_from_the_reference_it_searched(references, tmp_path):
  # k=6 is the index this search queries; the decoy's 2..6-mer indexes share the 'ref_'
  # prefix. Every hit must carry the reference's own accessions and species, and the hit
  # set must be complete -- nothing dropped by protecting the join.
  df = Matcher(
    query=['MKVLAT', 'GGWQMK', 'FFLLNN'],
    proteome_file=str(references),
    k=6,
    preprocessed_files_path=str(tmp_path),
    output_format='dataframe',
  ).match()

  assert hits(df) == {
    ('MKVLAT', 'P11111.1', 'MKVLAT'),
    ('MKVLAT', 'P22222.1', 'MKVLAT'),
    ('GGWQMK', 'P11111.1', 'GGWQMK'),
    ('FFLLNN', 'P22222.1', 'FFLLNN'),
  }
  assert set(df['Species']) == {'Homo sapiens'}
  assert not DECOY_IDS & set(df['Protein ID'])


def test_best_match_multi_k_walk_reads_metadata_from_the_reference(references, tmp_path):
  # best_match without an explicit k walks several k, so the metadata source has to be
  # one of the indexes that walk actually searched, not whatever the directory offers.
  df = Matcher(
    query=['MKVLAT', 'HFDSAM'],
    proteome_file=str(references),
    best_match=True,
    preprocessed_files_path=str(tmp_path),
    output_format='dataframe',
  ).match()

  assert set(df['Protein ID']) <= REFERENCE_IDS
  assert set(df['Species']) == {'Homo sapiens'}
  assert set(df['Query Sequence']) == {'MKVLAT', 'HFDSAM'}


def test_indel_search_reads_metadata_from_the_reference(references, tmp_path):
  # The indel path partitions the batch by required k and derives k itself; its metadata
  # must still come from an index it searched. MKVLT is MKVLAT with A deleted.
  df = Matcher(
    query=['MKVLT'],
    proteome_file=str(references),
    max_indels=1,
    preprocessed_files_path=str(tmp_path),
    output_format='dataframe',
  ).match()

  assert hits(df) == {
    ('MKVLT', 'P11111.1', 'MKVLAT'),
    ('MKVLT', 'P22222.1', 'MKVLAT'),
  }
  assert set(df['Species']) == {'Homo sapiens'}
