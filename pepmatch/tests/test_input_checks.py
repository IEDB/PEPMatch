"""Inputs PEPMatch used to mishandle silently or crash on: damaged or stale indexes, FASTA
files with empty headers, stray lines or lowercase residues, invalid parameters and queries,
discontinuous epitopes written without spaces, and custom proteome names."""
import hashlib
import json
import os
import sys

import polars as pl
import pytest
from pepmatch import Matcher, Preprocessor, preprocessor
from pepmatch.shell import run_matcher

PROTEOME = (
  '>sp|P11111|A_HUMAN Alpha OS=Homo sapiens OX=9606 GN=AAA PE=1 SV=1\nMKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQ\n'
  '>sp|P22222|B_HUMAN Beta OS=Homo sapiens OX=9606 GN=BBB PE=1 SV=1\nWWDDWWDDWWDDHHCCHHCCHHCCMMNNMMNN\n'
)


def proteome(tmp_path, text=PROTEOME, name='proteome.fasta'):
  path = tmp_path / name
  path.write_text(text)
  return path


def search(path, tmp_path, query, **kwargs):
  return Matcher(query=query, proteome_file=str(path), preprocessed_files_path=str(tmp_path), **kwargs).match()


def hits(df):
  rows = df.filter(pl.col('Matched Sequence').is_not_null())
  return sorted(zip(rows['Query Sequence'], rows['Protein ID'], rows['Index start']))


# --- damaged and stale indexes ----------------------------------------------------------

def damage(index, how):
  data = index.read_bytes()
  if how == 'truncated':
    index.write_bytes(data[:len(data) // 2])
  elif how == 'empty':
    index.write_bytes(b'')
  elif how == 'foreign':
    index.write_bytes(b'not an index ' * 50)
  elif how == 'one byte short':
    index.write_bytes(data[:-1])


@pytest.mark.parametrize('how', ['truncated', 'empty', 'foreign', 'one byte short'])
def test_a_damaged_index_raises_a_clear_error(tmp_path, how):
  path = proteome(tmp_path)
  Preprocessor(str(path), preprocessed_files_path=str(tmp_path)).preprocess(9)
  damage(tmp_path / 'proteome_9mers.pepidx', how)
  with pytest.raises(ValueError, match='Delete it and PEPMatch will rebuild it'):
    search(path, tmp_path, ['MKTAYIAKQ'])


def test_an_index_holding_another_k_raises(tmp_path):
  path = proteome(tmp_path)
  Preprocessor(str(path), preprocessed_files_path=str(tmp_path)).preprocess(3)
  os.replace(tmp_path / 'proteome_3mers.pepidx', tmp_path / 'proteome_9mers.pepidx')
  with pytest.raises(ValueError, match='holds 3-mers, not 9-mers'):
    search(path, tmp_path, ['MKTAYIAKQ'])


def test_a_missing_proteome_raises_file_not_found(tmp_path):
  with pytest.raises(FileNotFoundError, match='Proteome file not found'):
    search(tmp_path / 'nope.fasta', tmp_path, ['MKTAYIAKQ'])


def test_an_existing_index_still_works_without_its_fasta(tmp_path):
  path = proteome(tmp_path)
  Preprocessor(str(path), preprocessed_files_path=str(tmp_path)).preprocess(9)
  path.unlink()
  assert hits(search(path, tmp_path, ['MKTAYIAKQ'])) == [('MKTAYIAKQ', 'P11111.1', 1)]


def build(path, tmp_path, k=9):
  Preprocessor(str(path), preprocessed_files_path=str(tmp_path)).preprocess(k)
  return tmp_path / f'{path.stem}_{k}mers.pepidx'


def later(path):
  """Move the file's modification time forward, as a copy or a fresh checkout does."""
  stamp = os.stat(path).st_mtime_ns + 10**9
  os.utime(path, ns=(stamp, stamp))


GAMMA = '>sp|P33333|C_HUMAN Gamma OS=Homo sapiens OX=9606 GN=CCC PE=1 SV=1\nGGGGMKTAYIAKQGGGG\n'


def test_building_an_index_records_its_fasta(tmp_path):
  path = proteome(tmp_path)
  index = build(path, tmp_path)
  record = json.loads((tmp_path / 'proteome_9mers.pepidx.source').read_text())
  assert record['fasta_size'] == os.path.getsize(path)
  assert record['fasta_sha256'] == hashlib.sha256(path.read_bytes()).hexdigest()
  assert index.is_file()


def test_a_fasta_replaced_by_new_content_rebuilds_its_index(tmp_path, capsys):
  path = proteome(tmp_path)
  build(path, tmp_path)
  path.write_text(GAMMA)
  capsys.readouterr()
  assert hits(search(path, tmp_path, ['MKTAYIAKQ'])) == [('MKTAYIAKQ', 'P33333.1', 5)]
  assert 'changed since' in capsys.readouterr().out


def test_a_fasta_with_the_same_content_and_a_new_time_is_not_rebuilt(tmp_path, capsys):
  # A copy or a fresh checkout changes the time but not the content: no rebuild, and the
  # record takes the new time so the next search does not checksum again.
  path = proteome(tmp_path)
  index = build(path, tmp_path)
  built = os.stat(index).st_mtime_ns
  later(path)
  capsys.readouterr()
  assert hits(search(path, tmp_path, ['MKTAYIAKQ'])) == [('MKTAYIAKQ', 'P11111.1', 1)]
  assert 'building' not in capsys.readouterr().out
  assert os.stat(index).st_mtime_ns == built
  record = json.loads((tmp_path / 'proteome_9mers.pepidx.source').read_text())
  assert record['fasta_mtime_ns'] == os.stat(path).st_mtime_ns


def test_an_untouched_fasta_is_not_checksummed(tmp_path, monkeypatch):
  path = proteome(tmp_path)
  build(path, tmp_path)
  monkeypatch.setattr(preprocessor, 'fasta_sha256', lambda fasta: pytest.fail('checksummed an unchanged FASTA'))
  search(path, tmp_path, ['MKTAYIAKQ'])


def test_an_index_without_a_record_is_kept_when_its_counts_match(tmp_path, capsys):
  # Indexes built before records existed have none.
  path = proteome(tmp_path)
  index = build(path, tmp_path)
  (tmp_path / 'proteome_9mers.pepidx.source').unlink()
  later(path)
  built = os.stat(index).st_mtime_ns
  capsys.readouterr()
  search(path, tmp_path, ['MKTAYIAKQ'])
  assert 'building' not in capsys.readouterr().out and os.stat(index).st_mtime_ns == built
  assert (tmp_path / 'proteome_9mers.pepidx.source').is_file()


def test_an_index_without_a_record_is_rebuilt_when_its_counts_differ(tmp_path):
  path = proteome(tmp_path)
  build(path, tmp_path)
  (tmp_path / 'proteome_9mers.pepidx.source').unlink()
  path.write_text(GAMMA)
  assert hits(search(path, tmp_path, ['MKTAYIAKQ'])) == [('MKTAYIAKQ', 'P33333.1', 5)]


# --- FASTA reading -------------------------------------------------------------------------

def test_a_record_with_an_empty_header_is_its_own_protein(tmp_path):
  # It used to be glued onto the next record, shifting that record's coordinates.
  path = proteome(tmp_path, '>\nAAAAWWWW\n>p2 second\nMKTAYIAKQ\n')
  df = search(path, tmp_path, ['WWWWMKTA', 'MKTAYIAKQ'])
  assert hits(df) == [('MKTAYIAKQ', 'p2.1', 1)]


def test_sequence_before_the_first_header_is_left_out_with_a_warning(tmp_path, capsys):
  path = proteome(tmp_path, 'ORPHANSEQ\n\n>p1 first\nMKTAYIAKQ\n')
  df = search(path, tmp_path, ['ORPHANSEQMKTAYIAKQ', 'MKTAYIAKQ'])
  assert hits(df) == [('MKTAYIAKQ', 'p1.1', 1)]
  assert '1 sequence line(s) before its first ">" header' in capsys.readouterr().out


def test_a_lowercase_proteome_matches_uppercase_queries(tmp_path):
  path = proteome(tmp_path, '>p1 lower\nmktayiakqrqisfvk\n')
  df = search(path, tmp_path, ['MKTAYIAKQ', 'rqisfvk'])
  assert hits(df) == [('MKTAYIAKQ', 'p1.1', 1), ('RQISFVK', 'p1.1', 10)]


# --- parameters and queries ----------------------------------------------------------------

@pytest.mark.parametrize('kwargs', [dict(max_mismatches=-1), dict(max_indels=-1)])
def test_a_negative_edit_budget_raises(tmp_path, kwargs):
  with pytest.raises(ValueError, match='must be >= 0'):
    search(proteome(tmp_path), tmp_path, ['MKTAYIAKQ'], **kwargs)


@pytest.mark.parametrize('bad', ['MKTAY IAKQ', 'MKTAY\tIAKQ', ''])
@pytest.mark.parametrize('mode', [dict(), dict(max_mismatches=1), dict(max_indels=1), dict(best_match=True)])
def test_an_unsearchable_query_is_a_miss_and_does_not_stop_the_run(tmp_path, capsys, bad, mode):
  df = search(proteome(tmp_path), tmp_path, ['MKTAYIAKQ', bad], **mode)
  assert hits(df) == [('MKTAYIAKQ', 'P11111.1', 1)]
  assert df.filter(pl.col('Query Sequence') == bad)['Matched Sequence'].to_list() == [None]
  out = capsys.readouterr().out
  assert '1 query(s) cannot be searched' in out and 'is empty or contains whitespace' in out


def test_an_unsearchable_query_does_not_lower_k_for_the_rest(tmp_path, capsys):
  # An empty query used to make the shortest length 0 and drag the batch down to k=2.
  search(proteome(tmp_path), tmp_path, ['MKTAYIAKQ', ''])
  assert 'k=9' in capsys.readouterr().out


def test_a_run_of_only_unsearchable_queries_returns_misses(tmp_path):
  df = search(proteome(tmp_path), tmp_path, ['', 'A B'])
  assert df.height == 2 and df['Matched Sequence'].is_null().all()


@pytest.mark.parametrize('queries', [['M1, T3'], ['M1, T3', '']])
def test_counts_only_without_searchable_peptides_returns_an_empty_table(tmp_path, queries):
  df = search(proteome(tmp_path), tmp_path, queries, counts_only=True)
  assert df.height == 0 and df.columns == ['Query ID', 'Query Sequence', 'Mismatches', 'Count']


def test_a_discontinuous_epitope_may_omit_the_spaces(tmp_path):
  path = proteome(tmp_path)
  spaced = search(path, tmp_path, ['M1, T3, Y5'])
  compact = search(path, tmp_path, ['M1,T3,Y5'])
  assert hits(compact) == hits(spaced) == [('M1, T3, Y5', 'P11111.1', 1)]


@pytest.mark.parametrize('bad', ['M-1, T3', 'M-1,T3'])
def test_a_negative_discontinuous_position_is_a_miss_with_a_warning(tmp_path, capsys, bad):
  df = search(proteome(tmp_path), tmp_path, ['M1, T3, Y5', bad])
  assert hits(df) == [('M1, T3, Y5', 'P11111.1', 1)]
  assert df.filter(pl.col('Query Sequence') == bad)['Matched Sequence'].to_list() == [None]
  assert 'negative position' in capsys.readouterr().out


# --- custom proteome names -----------------------------------------------------------------

def test_matcher_finds_an_index_built_under_a_custom_name(tmp_path, capsys):
  path = proteome(tmp_path)
  Preprocessor(str(path), proteome_name='custom', preprocessed_files_path=str(tmp_path)).preprocess(9)
  capsys.readouterr()
  df = search(path, tmp_path, ['MKTAYIAKQ'], proteome_name='custom')
  assert hits(df) == [('MKTAYIAKQ', 'P11111.1', 1)]
  assert 'building' not in capsys.readouterr().out
  assert sorted(p.name for p in tmp_path.glob('*.pepidx')) == ['custom_9mers.pepidx']


def test_the_cli_accepts_a_custom_proteome_name(tmp_path, monkeypatch):
  path = proteome(tmp_path)
  Preprocessor(str(path), proteome_name='custom', preprocessed_files_path=str(tmp_path)).preprocess(9)
  queries = tmp_path / 'q.txt'
  queries.write_text('MKTAYIAKQ\n')
  out = tmp_path / 'out'
  monkeypatch.setattr(sys, 'argv', [
    'pepmatch-match', '-q', str(queries), '-p', str(path), '-P', str(tmp_path),
    '-n', 'custom', '-f', 'csv', '-o', str(out),
  ])
  run_matcher()
  assert pl.read_csv(f'{out}.csv')['Protein ID'].to_list() == ['P11111.1']
  assert sorted(p.name for p in tmp_path.glob('*.pepidx')) == ['custom_9mers.pepidx']
