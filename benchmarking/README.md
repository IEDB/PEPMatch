# Benchmarking

Standardized framework for comparing peptide search tools. Results below were generated on an AMD Ryzen 5 5500U (6 cores / 12 threads).

## Results

Results now include the new Rust engine implementation vs the Python-only implementation (old) for PEPMatch.

### MHC Class I Dataset
1,000 9-mer peptides searched against the human proteome (~200,000 proteins) for exact matches.

| Method | Proteome Preprocessing (s) | Query Preprocessing (s) | Searching (s) | Total (s) | Recall (%) |
|--------|---------------------------|------------------------|---------------|-----------|------------|
| **PEPMatch** | **32.2** | **N/A** | **0.049** | **32.2** | **100** |
| PEPMatch (old) | 39.6 | N/A | 0.08 | 39.7 | 100 |
| NmerMatch | 53.7 | 0.006 | 12.3 | 66.0 | 100 |
| BLAST | 1.27 | N/A | 11.3 | 12.6 | 98.3 |
| DIAMOND | 0.25 | N/A | 5.01 | 5.26 | 1.5 |
| MMseqs2 | 2.65 | N/A | 0.50 | 3.15 | 0.0 |

### Neoepitopes Dataset
1,735 15-mer peptides searched against the human proteome with up to 3 mismatches.

| Method | Proteome Preprocessing (s) | Query Preprocessing (s) | Searching (s) | Total (s) | Recall (%) |
|--------|---------------------------|------------------------|---------------|-----------|------------|
| **PEPMatch** | **5.3** | **N/A** | **0.459** | **5.7** | **100** |
| PEPMatch (old) | 13.3 | N/A | 18.4 | 31.7 | 100 |
| NmerMatch | 50.4 | 0.002 | 40.1 | 90.5 | 100 |
| BLAST | 1.28 | N/A | 119.2 | 120.5 | 58.1 |
| DIAMOND | 0.24 | N/A | 4.93 | 5.17 | 34.0 |
| MMseqs2 | 2.28 | N/A | 0.59 | 2.87 | 24.6 |

### SARS-CoV-2 Dataset
628 peptides (8-15 mers) searched against betacoronaviruses with up to 2 mismatches.

| Method | Proteome Preprocessing (s) | Query Preprocessing (s) | Searching (s) | Total (s) | Recall (%) |
|--------|---------------------------|------------------------|---------------|-----------|------------|
| **PEPMatch** | **2.2** | **N/A** | **6.4** | **8.6** | **100** |
| PEPMatch (old) | 34.1 | N/A | 32.6 | 66.7 | 100 |
| NmerMatch | 214 | 0.003 | 21.3 | 235.9 | 100 |
| BLAST | 0.51 | N/A | 115.9 | 116.4 | 73.3 |
| DIAMOND | 0.15 | N/A | 3.45 | 3.60 | 6.4 |
| MMseqs2 | 1.83 | N/A | 0.77 | 2.60 | 7.4 |

### Milk Peptides Dataset
111 15-mer peptides searched against the human proteome for best match.

| Method | Proteome Preprocessing (s) | Query Preprocessing (s) | Searching (s) | Total (s) | Recall (%) |
|--------|---------------------------|------------------------|---------------|-----------|------------|
| **PEPMatch** | **69.6** | **N/A** | **3.7** | **73.3** | **100** |
| PEPMatch (old) | 45.5 | N/A | 600.3 | 645.8 | 100 |
| NmerMatch | 203.7 | 0.005 | 1,168 | 1,372 | 100 |
| BLAST | 1.30 | N/A | 105.3 | 106.6 | 84.2 |
| DIAMOND | 0.25 | N/A | 5.24 | 5.49 | 75.3 |
| MMseqs2 | 2.25 | N/A | 0.48 | 2.73 | 75.6 |

## Adding Your Tool

To benchmark your tool against PEPMatch, create a Python wrapper with a `Benchmarker` class:
```python
class Benchmarker:
  def __init__(self, benchmark, query, proteome, lengths, max_mismatches, method_parameters):
    """
    Args:
      benchmark: dataset name (mhc_ligands, milk, coronavirus, neoepitopes)
      query: path to query FASTA file
      proteome: path to proteome FASTA file
      lengths: list of peptide lengths in the query
      max_mismatches: max substitutions allowed (0 = exact, -1 = best match)
      method_parameters: dict of tool-specific parameters from benchmarking_parameters.json
    """
    self.query = query
    self.proteome = proteome
    self.max_mismatches = max_mismatches
    # raise ValueError if your tool can't handle this dataset:
    # if max_mismatches > 0:
    #   raise ValueError('MyTool cannot do mismatching.')
    # if max_mismatches == -1:
    #   raise ValueError('MyTool does not have best match.')

  def __str__(self):
    return 'Your Tool Name'

  def preprocess_proteome(self):
    """Preprocess the proteome. Raise TypeError if not applicable."""
    raise TypeError('MyTool does not preprocess proteomes.')

  def preprocess_query(self):
    """Preprocess the query. Raise TypeError if not applicable."""
    raise TypeError('MyTool does not preprocess queries.')

  def search(self):
    """Run the search. Return a pandas DataFrame with columns:
      Query Sequence, Matched Sequence, Protein ID, Index start
    Index start is 1-based."""
    ...
```

Then add your tool to `benchmarking_parameters.json`:
```json
{
  "name": "your_module_name",
  "text_shifting": 0,
  "method_parameters": {}
}
```

Place your wrapper in the `methods/` directory and run:
```bash
python benchmarking.py -b mhc_ligands
python benchmarking.py -b neoepitopes
python benchmarking.py -b coronavirus
python benchmarking.py -b milk
```

### Flags

* `-b` (Required): Dataset to benchmark (`mhc_ligands`, `milk`, `coronavirus`, `neoepitopes`)
* `-m`: Include memory benchmarking
* `-t`: Include text-shifting algorithms (Horspool, Boyer-Moore, KMP, Z-algorithm)
