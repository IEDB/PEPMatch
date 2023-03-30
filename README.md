<p align="center">
  <img src="docs/logo.png">
</p>

--------------------------------------------------------------------

![Unit Tests](https://github.com/IEDB/PEPMatch/actions/workflows/unittests.yml/badge.svg)

#### Author: Daniel Marrama

Peptide search against a reference proteome, or sets of proteins, with residue subtitutions.

Two step process: preprocessing and matching.

Preprocessed data is stored in a SQLite or pickle format and only has to be performed once.

As a competition to improve tool performance, we created a benchmarking framework with instructions [here](./benchmarking).

### Requirements

- Python 3.7+
- [Pandas](https://pandas.pydata.org/)
- [NumPy](https://numpy.org/)
- [Biopython](https://biopython.org/)
- [Levenshtein](https://pypi.org/project/python-Levenshtein/)
- [redis](https://redis.com/)

### Installation

```
pip install pepmatch
```


### Inputs

#### Preprocessor

```proteome``` - Path to proteome file to search against.\
```k``` - k-mer size to break up proteome into.\
```preprocessed_format``` - SQLite ("sqlite") or "pickle".\
```preprocessed_files_path``` - (optional) Directory where you want preprocessed files to go. Default is current directory.\
```gene_priority_proteome``` - (optional) Subset of ```proteome``` with prioritized protein IDs.\

#### Matcher

```query``` - Query of peptides to search either in .fasta file or as a Python list.\
```proteome_file``` - Name of preprocessed proteome to search against.\
```max_mismatches``` - Maximum number of mismatches (substitutions) for query.\
```k``` - (optional) k-mer size of the preprocessed proteome. If no k is selected, then a best k will be calculated and the proteome will be preprocessed\
```preprocessed_files_path``` - (optional) Directory where preprocessed files are. Default is current directory.\
```best_match``` - (optional) Returns only one match per query peptide. It will output the best match.\
```output_format``` - (optional) Outputs results into a file (CSV, XLSX, JSON, HTML) or just as a dataframe.\
```output_name``` - (optional) Specify name of file for output. Leaving blank will generate a name.

Note: For now, due to performance, SQLite is used for exact matching and pickle is used for mismatching.

Note: PEPMatch can also search for discontinuous epitopes in the residue:index format. Example: 

"R377, Q408, Q432, H433, F436, V441, S442, S464, K467, K489, I491, S492, N497"

### Command Line Example

```bash
pepmatch-preprocess -p human.fasta -k 5 -f sql
pepmatch-match -q peptides.fasta -p human.fasta -m 0 -k 5
```

### Exact Matching Example

```python
from pepmatch import Preprocessor, Matcher

# proteome, k, preprocessed_format, target directory, gene_priority_proteome
Preprocessor('proteomes/human.fasta', '.' 'proteomes/human_gp.fasta').preprocess('sql', k=5)
# PREPROCESSING ONLY NEEDS TO BE DONE ONCE!

# query, proteome, max_mismatches, k, preprocessed files directory
Matcher('queries/mhc_ligands_test.fasta', 'proteomes/human.fasta', 0, 5, '.').match()
```

### Mismatching Example 

```python
from pepmatch import Preprocessor, Matcher

# proteome, k, preprocessed_format, target directory
Preprocessor('proteomes/human.fasta').preprocess('pickle', k=3)
# PREPROCESSING ONLY NEEDS TO BE DONE ONCE!

# query, proteome, max_mismatches, k, preprocessed files directory
Matcher('queries/neoepitopes_test.fasta', 'proteomes/human.fasta', 3, 3).match()
```

### Best Match Example

```python
from pepmatch import Preprocessor, Matcher
Matcher('queries/milk_peptides.fasta', 'proteomes/human.fasta', best_match=True).match()
```

The best match parameter without k or mismatch inputs will produce the best match for each peptide in the query, meaning the match with the least number of mismatches, the best protein existence level, and if the match exists in the gene priority proteome. No preprocessing beforehand is required, as the Matcher class will do this for you to find the best match.

### Outputs

As mentioned above, outputs can be specified with the ```output_format``` parameter in the ```Matcher``` class. The following formats are allowed: `dataframe`, `csv`, `xlsx`, `json`, and `html`.

If specifying `dataframe`, the ```match()``` method will return a pandas dataframe which can be stored as a variable:

```python
df = Matcher('queries/neoepitopes_test.fasta', 'human.fasta', 3, 3, output_format='dataframe').match()
```

### TODO
- Test other key-value stores (Redis, Memcached, LMDB, etc.)
- Remove dependency on Levenshtein (this is not maintained very well)