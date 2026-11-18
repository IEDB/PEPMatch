import hashlib
import json
import os
import struct
from pathlib import Path
from ._rs import rs_preprocess

# preprocess.rs packs each k-mer into a u128, one byte per residue, so no index can be
# built for k above 16.
MAX_INDEX_K = 16

# The .pepidx header: magic, k, protein count, sequence length, unique k-mers, positions,
# metadata length. Every section's size follows from it (see preprocess.rs write_pepidx).
PEPIDX_MAGIC = b'PEPIDX05'
PEPIDX_HEADER = struct.Struct('<8sBIQQQQ')


def check_index(path, k):
  """Raise a ValueError if `path` is not a complete k-mer index for this k. The search reads
  the file by fixed offsets, so a truncated, foreign or mislabeled file would otherwise crash
  inside the engine or return wrong hits."""
  with open(path, 'rb') as f:
    head = f.read(PEPIDX_HEADER.size)
  if len(head) < PEPIDX_HEADER.size or head[:8] != PEPIDX_MAGIC:
    problem = 'is not a PEPMatch index from this version'
  else:
    _, index_k, proteins, residues, kmers, positions, metadata = PEPIDX_HEADER.unpack(head)
    expected = (PEPIDX_HEADER.size + residues + 16 * proteins + metadata
                + kmers * (8 + index_k + 4) + 8 * positions)
    size = os.path.getsize(path)
    if index_k != k:
      problem = f'holds {index_k}-mers, not {k}-mers'
    elif size != expected:
      problem = f'is {size} bytes but its header describes {expected}, so it is incomplete or damaged'
    else:
      return
  raise ValueError(f'{path} {problem}. Delete it and PEPMatch will rebuild it.')


def build_index(fasta, k, path):
  """Build the k-mer index for `fasta` at `path`. Sequence lines before the first header
  belong to no record, so the engine leaves them out; say so rather than drop them silently."""
  if not os.path.isfile(fasta):
    raise FileNotFoundError(f'Proteome file not found: {fasta}')
  orphan = 0
  with open(fasta) as f:
    for line in f:
      if line.startswith('>'):
        break
      orphan += bool(line.strip())
  if orphan:
    print(f'Warning: {fasta} has {orphan} sequence line(s) before its first ">" header; '
          f'they belong to no protein and are left out of the index.')
  rs_preprocess(str(fasta), k, str(path))
  record_source(path, fasta)


# Next to each index, a record of the FASTA it was built from: size, modification time and
# sha256. An index is named only by the FASTA's name and k, so without it a FASTA replaced by
# a new version under the same name would be searched through its old index.
SOURCE_SUFFIX = '.source'


def fasta_sha256(fasta):
  digest = hashlib.sha256()
  with open(fasta, 'rb') as f:
    for block in iter(lambda: f.read(1 << 20), b''):
      digest.update(block)
  return digest.hexdigest()


def record_source(path, fasta, sha256=None):
  """Write the source record for the index at `path`. Best effort: an index in a read-only
  directory still works, its FASTA is just checksummed again on the next search."""
  stat = os.stat(fasta)
  record = {'fasta_size': stat.st_size, 'fasta_mtime_ns': stat.st_mtime_ns,
            'fasta_sha256': sha256 or fasta_sha256(fasta)}
  try:
    with open(str(path) + SOURCE_SUFFIX, 'w') as f:
      json.dump(record, f)
  except OSError:
    pass


def _read_source(path):
  try:
    with open(str(path) + SOURCE_SUFFIX) as f:
      return json.load(f)
  except (OSError, ValueError):
    return None


def _fasta_counts(fasta):
  """(proteins, residues) as the engine reads the FASTA: every '>' line starts a record and
  lines before the first one belong to none."""
  proteins = residues = 0
  with open(fasta, 'rb') as f:
    for line in f:
      if line.startswith(b'>'):
        proteins += 1
      elif proteins:
        residues += len(line.strip())
  return proteins, residues


def _index_counts(path):
  with open(path, 'rb') as f:
    _, _, proteins, residues, _, _, _ = PEPIDX_HEADER.unpack(f.read(PEPIDX_HEADER.size))
  return proteins, residues


def built_from(path, fasta):
  """Whether the index at `path` was built from this FASTA's current content.

  The size and modification time are compared with the record first, which costs nothing
  and is the answer almost every time. Only when they differ (the file was replaced, copied
  or checked out again) is the FASTA checksummed: the same content means only the time
  moved, so the record is updated and nothing is rebuilt. An index from before records
  existed is accepted if its protein count and length match the FASTA's."""
  stat = os.stat(fasta)
  record = _read_source(path)
  if record and (record.get('fasta_size'), record.get('fasta_mtime_ns')) == (stat.st_size, stat.st_mtime_ns):
    return True
  if record:
    sha256 = fasta_sha256(fasta)
    if sha256 != record.get('fasta_sha256'):
      return False
  elif _fasta_counts(fasta) != _index_counts(path):
    return False
  else:
    sha256 = None
  record_source(path, fasta, sha256)
  return True


class Preprocessor:
  """Preprocesses a proteome FASTA file into a .pepidx binary index
  optimized for fast peptide searching. The index stores k-mer positions,
  protein sequences, and metadata in a single memory-mapped file."""

  def __init__(
    self,
    proteome,
    proteome_name='',
    preprocessed_files_path='.',
  ):
    if not os.path.isfile(str(proteome)):
      raise FileNotFoundError(f'Proteome file not found: {proteome}')

    os.makedirs(preprocessed_files_path, exist_ok=True)
    self.preprocessed_files_path = preprocessed_files_path
    self.proteome_file = str(proteome)

    if not proteome_name:
      self.proteome_name = Path(proteome).stem
    else:
      self.proteome_name = proteome_name

  def preprocess(self, k: int) -> None:
    """Build the k-mer index and write it to <proteome_name>_<k>mers.pepidx."""
    if k < 2:
      raise ValueError('k must be >= 2.')
    if k > MAX_INDEX_K:
      raise ValueError(f'k must be <= {MAX_INDEX_K}, the longest k-mer the index can store.')
    output_path = os.path.join(
      self.preprocessed_files_path, f'{self.proteome_name}_{k}mers.pepidx'
    )
    print(f"Preprocessing {self.proteome_name} with k={k}...")
    build_index(self.proteome_file, k, output_path)
    print(f"Saved to {output_path}")
