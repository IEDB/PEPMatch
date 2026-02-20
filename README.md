<p align="center">
  <img src="docs/logo.png" alt="PEPMatch Logo">
</p>

--------------------------------------------------------------------

[![Unit Tests](https://github.com/IEDB/PEPMatch/actions/workflows/tests.yml/badge.svg)](https://github.com/IEDB/PEPMatch/actions/workflows/tests.yml)

**Author:** Daniel Marrama

`PEPMatch` is a high-performance peptide search tool for finding short peptide sequences within a reference proteome. Powered by a Rust engine with Python bindings, it delivers sub-second search times across entire proteomes while maintaining a simple Python API.

### Key Features

* **Blazing Fast**: Rust-powered search engine with automatic multi-core parallelization via Rayon. Search thousands of peptides against the entire human proteome in seconds.
* **Unified Index Format**: Single `.pepidx` binary format stores sequences, metadata, and k-mer index in one memory-mapped file. Preprocess once, search repeatedly.
* **Versatile Searching**: Exact matches, mismatch-tolerant searches, best match mode, and discontinuous epitope support.
* **Simple API**: Two classes — `Preprocessor` and `Matcher` — handle everything.
* **Flexible I/O**: Accepts queries from FASTA files, text files, or Python lists. Outputs to CSV, TSV, XLSX, JSON, or Polars DataFrame.

### Requirements

* Python 3.10+
* [Polars](https://pola.rs/)
* [Biopython](https://biopython.org/)

### Installation
```bash
pip install pepmatch
```

### Quick Start
```python
from pepmatch import Preprocessor, Matcher

# Preprocess a proteome (one-time step)
Preprocessor('human.fasta').preprocess(k=5)

# Search for exact matches
df = Matcher(
  query=['YLLDLHSYL', 'GLCTLVAML', 'FAKEPEPTIDE'],
  proteome_file='human.fasta',
  max_mismatches=0,
  k=5
).match()

print(df)
```

### Preprocessing

Preprocessing builds a `.pepidx` index from your proteome FASTA file. This only needs to be done once per proteome and k-mer size. If a `.pepidx` file doesn't exist when you search, `Matcher` will create it automatically.
```python
from pepmatch import Preprocessor

Preprocessor('human.fasta').preprocess(k=5)
```

**CLI:**
```bash
pepmatch-preprocess -p human.fasta -k 5
```

#### Flags

* `-p`, `--proteome` (Required): Path to the proteome FASTA file.
* `-k`, `--kmer_size` (Required): The k-mer size for indexing.
* `-n`, `--proteome_name`: Custom name for the proteome.
* `-P`, `--preprocessed_files_path`: Directory to save preprocessed files.

### Matching

#### Exact Matching
```python
from pepmatch import Matcher

df = Matcher(
  query='peptides.fasta',
  proteome_file='human.fasta',
  max_mismatches=0,
  k=5
).match()
```

#### Mismatch Searching
```python
df = Matcher(
  query='neoepitopes.fasta',
  proteome_file='human.fasta',
  max_mismatches=3,
  k=3
).match()
```

#### Best Match

Automatically finds the optimal match for each peptide by trying different k-mer sizes and mismatch thresholds. No manual preprocessing required.
```python
df = Matcher(
  query='peptides.fasta',
  proteome_file='human.fasta',
  best_match=True
).match()
```

#### Discontinuous Epitope Searching

Search for non-contiguous residues defined by their positions.
```python
df = Matcher(
  query=[
    "R377, Q408, Q432, H433, F436",
    "S2760, V2763, E2773, D2805, T2819"
  ],
  proteome_file='sars-cov-2.fasta',
  max_mismatches=1
).match()
```

#### Mixed Queries

Linear peptides and discontinuous epitopes can be searched together.
```python
df = Matcher(
  query=[
    'YLLDLHSYL',
    'R377, Q408, Q432, H433, F436',
    'GLCTLVAML',
  ],
  proteome_file='sars-cov-2.fasta',
  max_mismatches=0
).match()
```

#### Query Input Formats

* **Python list**: `['YLLDLHSYL', 'GLCTLVAML']`
* **FASTA file**: `.fasta`, `.fas`, `.fa`, `.fna`, `.ffn`, `.faa`, `.mpfa`, `.frn`
* **Text file**: `.txt` with one peptide per line

**CLI:**
```bash
pepmatch-match -q peptides.fasta -p human.fasta -m 0 -k 5
```

#### Flags

* `-q`, `--query` (Required): Path to the query file.
* `-p`, `--proteome_file` (Required): Path to the proteome FASTA file.
* `-m`, `--max_mismatches`: Maximum mismatches allowed (default: 0).
* `-k`, `--kmer_size`: K-mer size (default: 5).
* `-P`, `--preprocessed_files_path`: Directory containing preprocessed files.
* `-b`, `--best_match`: Enable best match mode.
* `-f`, `--output_format`: Output format — `csv`, `tsv`, `xlsx`, `json` (default: `csv`).
* `-o`, `--output_name`: Output file name (without extension).
* `-v`, `--sequence_version`: Disable sequence versioning on protein IDs.

### Output Formats

* **`dataframe`** (default for API): Returns a Polars DataFrame.
* **`csv`** (default for CLI): CSV file.
* **`tsv`**: Tab-separated file.
* **`xlsx`**: Excel file.
* **`json`**: JSON file.

### Performance

Benchmarked searching ~2,000 peptides against the human proteome (~200,000 proteins):

| Mode | Time |
|------|------|
| Exact match (k=5) | ~0.06s |
| 1 mismatch (k=3) | ~1.5s |
| 2 mismatches (k=3) | ~1.9s |
| 3 mismatches (k=3) | ~3.7s |

### Citation

If you use PEPMatch in your research, please cite:

Marrama D, Chronister WD, Westernberg L, et al. PEPMatch: a tool to identify short peptide sequence matches in large sets of proteins. *BMC Bioinformatics*. 2023;24(1):485. Published 2023 Dec 18. doi:10.1186/s12859-023-05606-4
