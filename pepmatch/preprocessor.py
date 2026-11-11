import os
from pathlib import Path
from ._rs import rs_preprocess

# preprocess.rs packs each k-mer into a u128, one byte per residue, so no index can be
# built for k above 16.
MAX_INDEX_K = 16

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
    rs_preprocess(self.proteome_file, k, output_path)
    print(f"Saved to {output_path}")
