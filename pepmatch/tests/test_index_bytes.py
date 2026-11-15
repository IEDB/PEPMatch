"""Byte-level guard for the index builder.

Searching reads a .pepidx by fixed offsets, so the builder's contract is the file itself. The
checksums below were recorded from the original builder (hash map per batch, chunk files,
k-way merge); any builder must reproduce them exactly. Rebuilding an index must therefore
never be needed after a builder change, and a search cannot tell which builder wrote its index.
"""
import hashlib
from pathlib import Path

import pytest
from pepmatch import Preprocessor

DATA = Path(__file__).parent / 'data'
KS = (2, 3, 5, 9, 16)

# UniProt and plain headers, an isoform, a homopolymer and a periodic repeat, proteins shorter
# than k, and a sequence wrapped over several lines.
EDGE_FASTA = (
  '>sp|P11111|AAAA_HUMAN Alpha protein OS=Homo sapiens OX=9606 GN=AAA PE=1 SV=2\n'
  'MKVLATGGWQMKRRDEYHAAAAAAAAAAAAAAAAAAAAMKVLATGGWQ\n'
  '>sp|P11111-2|AAAA_HUMAN Isoform 2 of Alpha protein OS=Homo sapiens OX=9606 GN=AAA\n'
  'MKVLATGGWQMKRR\n'
  '>tr|A0A000|A0A000_HUMAN Beta protein (Fragment) OS=Homo sapiens OX=9606 GN=BBB PE=4 SV=1\n'
  'EKEKEKEKEKEKEKEKEKEK\n'
  '>plain_header no uniprot fields\n'
  'ACDEFGHIKLMNPQRSTVWY\n'
  '>short\n'
  'AC\n'
  '>single\n'
  'M\n'
  '>wrapped\n'
  'MKTAYIAKQR\n'
  'QISFVKSHFS\n'
  'RQLEERLGLI\n'
)


def many_proteins_fasta():
  """60,000 short proteins: more than one 50,000-protein batch of the original builder, so
  the positions of a k-mer came from different batches and had to merge in protein order.
  A fixed linear congruential generator keeps the file identical on every Python."""
  residues = 'ACDEFGHIKLMNPQRSTVWY'
  state = 20261113
  lines = []
  for number in range(60_000):
    state = (state * 1103515245 + 12345) % 2**31
    length = 3 + (state >> 8) % 12
    sequence = []
    for _ in range(length):
      state = (state * 1103515245 + 12345) % 2**31
      sequence.append(residues[(state >> 16) % 20])
    lines.append(f'>p{number + 1}\n{"".join(sequence)}\n')
  return ''.join(lines)


EXPECTED = {
  ('proteome', 2): '318225e173586005b29ef05ab08ec9f79052192a90fa3cedabf61ec5730f9dd7',
  ('proteome', 3): 'f74f67104b7e1f6a7ebf43ba3b2743f6f012d415f6305cd446492bc4a24bd37b',
  ('proteome', 5): 'e059a945d6e9de6b95ad53fe64c8131f68565e0321b6a3f1248b036f23dd7473',
  ('proteome', 9): 'f6a9cb1b85a28847854eba5d05c6eced9fd60bd2cef25719d2aa60401711dec5',
  ('proteome', 16): 'fd3f66df7714efd7f316fc13992ae352604dc7645818c362b37c548b4e3ddc47',
  ('indel_golden_proteome', 2): 'a58988ffe7c674110dfa6b27e389a207d63c75570a750bcb6612a91e4d32e06c',
  ('indel_golden_proteome', 3): '2ae825b16b4856040b601cd20ce47c1c920ce552e1fc8b0dab0363debdb94037',
  ('indel_golden_proteome', 5): 'b429dc2ed53fca25e998c0440a27ee70eb4ac8ca42d098fb015e79ad9e0dce84',
  ('indel_golden_proteome', 9): 'af332df041ab85e8ab76c15ea34d7109bcdb0ac0ed442929d003613a0d9ecd98',
  ('indel_golden_proteome', 16): '2516abd6217f219e427290b227e9fa3715c8e74de14aeda7f64816f51f1c0d5b',
  ('edge', 2): '82baeb30421663af6ad2b31e284b840d6b40b29c589fccb1220527eaba906087',
  ('edge', 3): 'be86930b8aae0e7c0d234c83ae5b2947591b845ce6901ca637c250271d783c2e',
  ('edge', 5): '079ac7f1697fc39708385f28fbd8a774c60dd611fc7902db026eebc89854bfcf',
  ('edge', 9): '1b518128b7b4ecd63ec09c4387c4bcc6f529192a52b45d0d4b8ae343dbe9b44b',
  ('edge', 16): 'fc6f4a11fe5be44e0c873739ecc64ad4148634c3c8fc431da640a5bbd6c043cb',
  ('many', 2): '6bfc8e00afb35165383892ce6c494718576b5cb8665beed08c3d049d32ff0e5e',
  ('many', 3): '8e39a6c9b96a22115d7517acf2f0f0f959888c06036a2b76d63510ec98f9c1bb',
  ('many', 5): '93d4ebe0628d4a84882643ca48784696204a20cb64a60771ca88babf5220377d',
  ('many', 9): 'c59cbf7106513cbce0edc37a55f42170fe5d72a6442959e3679d6ad3f7ab793f',
  ('many', 16): '32b161b52406698903c13ca0694e456a9c288af2b4c7cdf417962e1361f51297',
}


def source_fasta(name, tmp_path):
  if name == 'edge':
    path = tmp_path / 'edge.fasta'
    path.write_text(EDGE_FASTA)
    return path
  if name == 'many':
    path = tmp_path / 'many.fasta'
    path.write_text(many_proteins_fasta())
    return path
  return DATA / f'{name}.fasta'


def index_sha256(name, k, tmp_path):
  fasta = source_fasta(name, tmp_path)
  Preprocessor(str(fasta), preprocessed_files_path=str(tmp_path)).preprocess(k)
  index = tmp_path / f'{fasta.stem}_{k}mers.pepidx'
  return hashlib.sha256(index.read_bytes()).hexdigest()


@pytest.mark.parametrize('name', ['proteome', 'indel_golden_proteome', 'edge', 'many'])
@pytest.mark.parametrize('k', KS)
def test_index_is_byte_identical_to_the_original_builder(name, k, tmp_path):
  assert index_sha256(name, k, tmp_path) == EXPECTED[(name, k)]


def test_build_leaves_no_chunk_merge_or_temporary_file(tmp_path):
  fasta = source_fasta('edge', tmp_path)
  Preprocessor(str(fasta), preprocessed_files_path=str(tmp_path)).preprocess(3)
  assert (tmp_path / 'edge_3mers.pepidx').is_file()
  assert not [p.name for p in tmp_path.iterdir() if '.chunk_' in p.name or '.merged' in p.name or '.tmp.' in p.name]
