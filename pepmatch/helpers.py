import re
import polars as pl

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


class TqdmDummy:
  """A dummy class that mimics tqdm for when it's not installed."""
  def __init__(self, *args, **kwargs): pass
  def update(self, n=1): pass
  def __enter__(self): return self
  def __exit__(self, exc_type, exc_val, exc_tb): pass


def parse_fasta(file: str) -> list:
  """Return a parsed Biopython SeqRecord object from a FASTA file.
  Args:
    file: path to FASTA file.
  """
  return list(SeqIO.parse(file, 'fasta'))


def split_sequence(sequence: str, k: int) -> list:
  """
  Splits a peptide into equal sized k-mers on a rolling basis.
  Ex: k = 4, NSLFLTDLY --> ['NSLF', 'SLFL', 'LFLT', 'FLTD', 'LTDL', 'TDLY']

  Args:
    sequence: peptide sequence.
    k: k-mer length.
  """
  return [sequence[i:i+k] for i in range(len(sequence) - k + 1)]


def extract_metadata(record: SeqRecord, header_id: bool) -> list:
  """Extract metadata from a FASTA header.
  Args: 
    record: protein SeqRecord from proteome FASTA file.
  """
  regexes = {
    'protein_id': re.compile(r"\|([^|]*)\|"),     # between | and |
    'protein_name': re.compile(r"\s(.+?)\sOS"),   # between space and space before OS
    'species': re.compile(r"OS=(.+?)\sOX"),       # between OS= and space before OX
    'taxon_id': re.compile(r"OX=(.+?)(\s|$)"),         # between OX= and space
    'gene': re.compile(r"GN=(.+?)(\s|$)"),             # between GN= and space
    'pe_level': re.compile(r"PE=(.+?)(\s|$)"),         # between PE= and space
    'sequence_version': re.compile(r"SV=(.+?)(\s|$)"), # between SV= and space
    'gene_priority': re.compile(r"GP=(.+?)(\s|$)"),    # between GP= and space
    'swissprot': re.compile(r"^(tr|sp)\|"),        # between > and |
  }
  metadata = []
  for key in regexes: # loop through compiled regexes to extract metadata
    match = regexes[key].search(str(record.description))
    if match:
      if key == 'swissprot':
        metadata.append('1') if match.group(1) == 'sp' else metadata.append('0')
      elif key == 'protein_id':
        metadata.append(match.group(1)) if not header_id else metadata.append(str(record.id))
      else:
        metadata.append(match.group(1))
    else:
      if key == 'swissprot':
        metadata.append('0')
      elif key == 'protein_id':
        metadata.append(str(record.id)) # get record.id from FASTA header instead
      elif key == 'sequence_version':
        metadata.append('1')
      elif key in ['pe_level', 'gene_priority']:
        metadata.append('0') # zeros for integer columns
      else:
        metadata.append('')  # empty strings for string columns

  return metadata


def output_matches(df: pl.DataFrame, output_format: str, output_name: str) -> None:
  df_to_write = df.clone()  # for files that can't do nested data, we convert mutated positions column to string
  if "Mutated Positions" in df_to_write.columns and output_format in ['csv', 'tsv', 'xlsx']:
    df_to_write = df_to_write.with_columns(
      pl.when(pl.col("Mutated Positions").list.len() > 0)
        .then(pl.format("[{}]", pl.col("Mutated Positions").list.eval(pl.element().cast(pl.Utf8)).list.join(", ")))
        .otherwise(pl.lit("[]"))
        .alias("Mutated Positions")
    )

  # appends '.' + filetype if the name does not already contain it
  path = output_name.__str__()
  if not path.lower().endswith(f".{output_format}"):
    path += f".{output_format}"

  if output_format == 'csv':
    df_to_write.write_csv(path)
  elif output_format == 'tsv':
    df_to_write.write_csv(path, separator='\t')
  elif output_format == 'xlsx':
    df_to_write.write_excel(path)
  elif output_format == 'json':
    df_to_write.write_json(path)

