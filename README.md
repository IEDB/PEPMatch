# PEPMatch

### Author: Daniel Marrama

Peptide search against reference proteome(s) with specified mismatches.


## Installation
```
pip install git+https://gitlab.lji.org/dmarrama/pepmatch.git
```

## How It Works
PEPMatch achieves highly sensitive and accurate searching of peptides and epitopes within a proteome by first performing a preprocessing step. The proteome in question in preprocessed by spliting it up into all possible k-mers for a given k and mapping them to the location of the protein and string index position within that protein. One all preprocessing is done, it does not have to be done again and can be stored either as a SQL database or simply as a pickle file. Preprocessing can be done as follows:

```
from pepmatch import Preprocessor
preprocessor = Preprocessor('proteome.fa', 5, 'sql') # pass path to file, k, and format
preprocessor.preprocess()
```

There are two ways that PEPMatch does matching: finding exact matches and finding matches with mismatches (substitutions). After a proteome has been preprocessed either in SQL or pickle format we can now submit a query to search against by specifying the # of mismatches:

```
from pepmatch import Matcher
matcher = Matcher('query.fa', 'proteome.fa', 3) # pass paths to query/proteome, and # of mismatches
results = matcher.match()
```

PEPMatch will output a file of your matches. You can specify the output format when initialzing the Matcher object:

```
# supports 'csv', xlsx', 'json' and 'html' 
matcher = Matcher('query.fa', 'proteome.fa', 3, output = True, output_format = 'csv')
matcher.match()
```

Note: SQL is used for exact matching for now and pickle format is used for mismatching. This will change in future versions.


## Applications
- Exact matching against a proteome (e.g. eluted MHC Ligands found in a proteome for IEDB)
- Similarity peptides against other proteomes (e.g. bovine milk peptides against human proteome)
- Peptide comparison of similar species (e.g. SARS-CoV-2 against other coronaviruses)
- Mutated peptide search (e.g. neoepitopes in cancer)


## Benchmarking
In order to improve on the tool, a system to compare current tools/algorithms (current list below) has been created. For each application above, certain benchmarks have been created with a specified query and proteome. To compare your tool/algorithm to the current ones, the following instructions should be implemented:

1. Tools should be broken down into initialization, query preprocessing, proteome preprocessing, and searching. Preprocessing is not required for every tool, but if it is, certain arguments must be passed in order to perform this task. Preprocessing is defined as performing some functions on the query or proteome before searching which helps speed this process up. Of course, searching is the actual runtime your tool performs to search each peptide through the specified proteome. In order to standardize the benchmarking smoothly, a Python wrapper between the benchmarking script and your tool's executable must be created.

2. A Python script called "benchmarking.py" is found within this repository and can be used to compare your tool with others. Currently, there are 4 benchmarking frameworks based on the 4 applications listed above. The names of the benchmarking frameworks are "mhc_ligands", "milk", "coronavirus" and "neoepitopes". In the benchmarking_parameters.json file, you can see the lengths of the query peptides, the mismatches required for the benchmarking, the query file, the proteome file, and the text file with the expected results from the search. Within this file is the list of algorithms which can be used for each application. Put your wrapper within the repository and add the module name to the algorithms list in the JSON file.

    You can call the script like this:

    ```
    ./benchmarking.py -d mhc_ligands -s -t
    ```
    Passing the "-d" argument with the framework you want to test your tool against will automatically pull the necessary parameters for that framework. Passing the "-s" argument in the command line will skip the memory benchmark which usually takes a long time. Note: the memory benchmark is limited since the recording method will only account for the process ID for the benchmarking script. In the future, the benchmarking may be put in a container in order to improve on this.

    Lastly, passing the "-t" argument will include the text-shifting algorithms, which are algorithms used in various datasets for searching text files. These algorithms are often slow, and do not do mismatching, so they are just an option.

    The following are the list of inputs for each dataset: 
    
    - The MHC ligand benchmarking will run 1,000 9-mers through the human proteome for exact matches (0 mismatches).
    - The milk benchmarking will run 111 15-mers through the human proteome for the best match. The parameter for mismatches is -1 in this case to indicate best match. If your tool does not have a best match optiion, raising an exception is sufficient and is described below.
    - The coronavirus benchmarking will run 628 peptides of varying lengths (8-15) through a large FASTA file of other coronaviruses for up to and including 2 mismatches.
    - The neoepitopes will run 1,735 15-mers through the human proteome for up to and including 3 mismatches.

    Simply pass "-d mhc_ligands", "-d milk", "-d coronavirus", or "-d neoepitopes" in the command line for the dataset you want to run.

3. Your wrapper needs to have a class called "Benchmarker" which will take the first two arguments for initilization: lengths and mismatches. The class can inherit from your actual tool or just be created as an object which will call your code's executable. The Benchmarker object will should then have the three methods: preprocess_proteome, preprocess_query, and search. Along with the three methods, you need an \_\_init\_\_ method taking in at least the lengths of the peptides and the corresponding mismatches. You will need a \_\_str\_\_ method that just returns the name of your tool as a string.

    The Python class with its methods should be written like this with the following arguments: 
    
    ```
    class Benchmarker(object):
        def __init__(self, lengths, mismatches):
            self.lengths = lengths
            self.mismatches = mismatches

        def __str__(self):
            return 'Name of your tool'

        def preprocess_proteome(self, proteome):
            ...
        def preprocess_query(self, query):
            ...
        def search(self, query, proteome):
            ...
    ```
    
    What you do with the arguments after they are passed is up to you. Having the other necessary arguments specific to your tool in the "benchmarking_parameters.json" file and importing them into your wrapper would be helpful so others can run your tool without having to alter your wrapper or code of your actual tool. The three methods are for the three separate calls that the script will perform in order to time each process.
    
4. If your tool does not do any mismatching, raise a ValueError within your \_\_init\_\_ method:

    Example:
    
    ```
    class Benchmarker(object):
        def __init__(self, lengths, mismatches, *args, **kwargs):
            self.lengths = lengths
            self.mismatches = mismatches
            
            if self.max_mismatches != 0:
                raise ValueError(self.__str__() + ' cannot do any mismatching.\n')
    ```
    
    And if your tool does not do query or proteome preprocessing, raise a TypeError within the method:
    
    ```
    class Benchmarker(object):
    .
    .
    .
        def preprocess_proteome(self, proteome):
            raise TypeError(self.__str__() + ' does not preprocess proteomes.\n')
    ```

5. Lastly the results from your search method should be returned in a standard format, which will then be put through another function for accuracy. The results should be in a Python list with the format as follows: comma-separated values with the query peptide first, followed by the matched peptide within the proteome, followed by the protein ID it is found in, followed by the number of mismatches, and then lastly, the index position where the peptide is found. Note: the index is 0-based, not 1-based, as is standard in Python.

    Example:
    
    ```
    YLLDLHSYL,YLLDLHSYL,sp|O60337|MARH6_HUMAN,0,561
    ```
    This example is also the first result from the eluted MHC ligands benchmarking. The benchmarking script will then use these results and compare them to expected values. The output will be a percentage of correct results.
    
6. Please take a look at the "benchmarking_example.py" script or the PEPMatch code to see how a wrapper should be look. If anything was unclear in this writeup, please email me at dmarrama@lji.org with any questions or bugs and I will clarify / make corrections.

Good luck!


## Current Algorithms
- Exact Matching
    - k-mer Mapping
    - BLAST
    - Text shifting algorithms
        - Horspool
        - Boyer-Moore
        - Knuth-Morris-Pratt
        - Z-algorithm
- Mismatching
    - Hamming distance (brute force)
    - k-mer Mapping
    - BLAST