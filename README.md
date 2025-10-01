<p align="center">
  <img src="docs/logo.png" alt="PEPMatch Logo">
</p>

--------------------------------------------------------------------

[![Unit Tests](https://github.com/IEDB/PEPMatch/actions/workflows/tests.yml/badge.svg)](https://github.com/IEDB/PEPMatch/actions/workflows/tests.yml)

**Author:** Daniel Marrama

`PEPMatch` is a high-performance Python tool designed to find short peptide sequences within a reference proteome or other large protein sets. It is optimized for speed and flexibility, supporting exact matches, searches with a defined number of residue substitutions (mismatches), and a "best match" mode to find the most likely hit.

As a competition to improve tool performance, we created a benchmarking framework with instructions [here](./benchmarking).

### Key Features

* **Versatile Searching**: Find exact matches, matches with a specified tolerance for mismatches, or the single best match for each query peptide.
* **Discontinuous Epitope Support**: Search for non-contiguous residues in the format `"R377, Q408, Q432, ..."`.
* **High Performance**: Utilizes an efficient k-mer indexing strategy for rapid searching. The backend is powered by a C-based Hamming distance calculation for optimized mismatch detection.
* **Optimized Preprocessing**: Employs a two-step process. Proteomes are preprocessed once into a format optimized for the search type (SQLite for exact matching, Pickle for mismatching), making subsequent searches extremely fast.
* **Parallel Processing**: Built-in support for multicore processing to handle large query sets efficiently.
* **Flexible I/O**: Accepts queries from FASTA files or Python lists and can output results to multiple formats, including CSV, TSV, XLSX, JSON, or directly as a Polars DataFrame.

### Requirements

* Python 3.7+
* [Polars](https://pola.rs/)
* [Biopython](https://biopython.org/)

### Installation

```bash
pip install pepmatch
```

### Core Engine

`PEPMatch` operates using a two-step workflow:

1.  **Preprocessing**: First, the target proteome is processed into an indexed format. This step only needs to be performed once per proteome and k-mer size. `PEPMatch` uses SQLite databases for the speed of indexed lookups in exact matching and serialized Python objects (pickle) for the flexibility needed in mismatch searching.
2.  **Matching**: The user's query peptides are then searched against the preprocessed proteome.

This design ensures that the time-intensive task of parsing and indexing the proteome is separated from the search itself, allowing for rapid and repeated querying.

### Command-Line Usage

The tool provides two CLI commands: `pepmatch-preprocess` and `pepmatch-match`.

#### 1. Preprocessing

The `pepmatch-preprocess` command builds the necessary database from your proteome FASTA file.

* For **exact matching** (0 mismatches), use the `sql` format.
* For **mismatch matching**, use the `pickle` format.

```bash
# Preprocess for an exact match search using 5-mers
pepmatch-preprocess -p human.fasta -k 5 -f sql

# Preprocess for a mismatch search using 3-mers
pepmatch-preprocess -p human.fasta -k 3 -f pickle
```

#### 2. Matching

The `pepmatch-match` command runs the search against a preprocessed proteome.

```bash
# Find exact matches (-m 0) using the preprocessed 5-mer database
pepmatch-match -q peptides.fasta -p human.fasta -m 0 -k 5

# Find matches with up to 3 mismatches (-m 3) using the 3-mer database
pepmatch-match -q neoepitopes.fasta -p human.fasta -m 3 -k 3
```

### Python API Usage

For more control and integration into other workflows, `PEPMatch` provides a simple Python API.

#### 1. Exact Matching

```python
from pepmatch import Preprocessor, Matcher

# Preprocess the proteome into a SQLite DB for exact matching
Preprocessor('proteomes/human.fasta').sql_proteome(k=5)

# Initialize the Matcher for an exact search (0 mismatches)
matcher = Matcher(
  query='queries/mhc-ligands-test.fasta',
  proteome_file='proteomes/human.fasta',
  max_mismatches=0,
  k=5
)

# Run the search and get results
results_df = matcher.match()
```

#### 2. Mismatching

```python
from pepmatch import Preprocessor, Matcher

# Preprocess the proteome into pickle files for mismatching
Preprocessor('proteomes/human.fasta').pickle_proteome(k=3)

# Initialize the Matcher to allow up to 3 mismatches
matcher = Matcher(
  query='queries/neoepitopes-test.fasta',
  proteome_file='proteomes/human.fasta',
  max_mismatches=3,
  k=3
)

results_df = matcher.match()
```

#### 3. Best Match

The `best_match` mode automatically finds the optimal match for each peptide, trying different k-mer sizes and mismatch thresholds. No manual preprocessing is required.

```python
from pepmatch import Matcher

matcher = Matcher(
  query='queries/milk-peptides-test.fasta',
  proteome_file='proteomes/human.fasta',
  best_match=True
)

results_df = matcher.match()
```

#### 4. Parallel Processing

Use the `ParallelMatcher` class to run searches on multiple CPU cores. The `n_jobs` parameter specifies the number of cores to use.

```python
from pepmatch import Preprocessor, ParallelMatcher

# Preprocessing is the same
Preprocessor('proteomes/betacoronaviruses.fasta').pickle_proteome(k=3)

# Use ParallelMatcher to search with 4 jobs
parallel_matcher = ParallelMatcher(
  query='queries/coronavirus-test.fasta',
  proteome_file='proteomes/betacoronaviruses.fasta',
  max_mismatches=3,
  k=3,
  n_jobs=4
)

results_df = parallel_matcher.match()
```

#### 5. Discontinuous Epitope Searching

`PEPMatch` can search for epitopes defined by non-contiguous residues and their positions. Simply provide a query list where each item is a string in the format `"A1, B10, C15"`.

```python
from pepmatch import Matcher

# A list of discontinuous epitopes to find
discontinuous_query = [
  "R377, Q408, Q432, H433, F436",
  "S2760, V2763, E2773, D2805, T2819"
]

matcher = Matcher(
  query=discontinuous_query,
  proteome_file='proteomes/sars-cov-2.fasta',
  max_mismatches=1  # Allow 1 mismatch among the specified residues
)

results_df = matcher.match()
```

### Output Formats

You can specify the output format using the `output_format` parameter in the `Matcher` or `ParallelMatcher`.

* **`dataframe` (default for API)**: Returns a Polars DataFrame.
* **`csv` (default for CLI)**: Saves results to a CSV file.
* **`tsv`**: Saves results to a TSV file.
* **`xlsx`**: Saves results to an Excel file.
* **`json`**: Saves results to a JSON file.

To receive a DataFrame from the API, you can either omit the `output_format` parameter or set it explicitly:

```python
# The match() method will return a Polars DataFrame
df = Matcher(
  'queries/neoepitopes-test.fasta',
  'proteomes/human.fasta',
  max_mismatches=3,
  k=3,
  output_format='dataframe' # Explicitly request a DataFrame
).match()

print(df.head())
```

### Citation

If you use PEPMatch in your research, please cite the following paper:

Marrama D, Chronister WD, Westernberg L, et al. PEPMatch: a tool to identify short peptide sequence matches in large sets of proteins. *BMC Bioinformatics*. 2023;24(1):485. Published 2023 Dec 18. doi:10.1186/s12859-023-05606-4