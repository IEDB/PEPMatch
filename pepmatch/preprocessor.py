import os
from ._rs import rs_preprocess

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
      self.proteome_name = str(proteome).split('/')[-1].split('.')[0]
    else:
      self.proteome_name = proteome_name

  def preprocess(self, k: int) -> None:
    if k < 2:
      raise ValueError('k must be >= 2.')
    output_path = os.path.join(
      self.preprocessed_files_path, f'{self.proteome_name}_{k}mers.pepidx'
    )
    print(f"Preprocessing {self.proteome_name} with k={k}...")
    rs_preprocess(self.proteome_file, k, output_path)
    print(f"Saved to {output_path}")
