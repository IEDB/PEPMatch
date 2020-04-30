# Nmer Match

## Dependencies
* Perl 5.28 or higher
* Perl modules in [requirements-notest.txt](requirements-notext.txt) (installed with cpanm -n)
* Perl modules in [requirements.txt](requirements.txt) (installed with cpanm)

## Installation
Installation is as easy as installing the required perl modules.  I recommend using [perlbrew](https://perlbrew.pl) to manage your perl environment and cpanm to install packages.  To install the requirements in your perl installation:

```bash
$ cat requirements.txt | xargs -I {} cpanm {}
```

Alternatively, you can use the Docker container, described below, that has all requierements installed.

## Usage
Identifying matches in a set of peptide sequences and a set of proteins (this can also work for sets of nucleotide sequences) is divided into two phases, __catalog building__ and __searching__.

### Catalog building
First, you should assemble the dataset that you want to search against.  We'll call this the protein _database_, which should be in fasta format ([example database](t/test_db.fasta)).  This can range from a small list of peptides, to a large proteome.

Building the catalog is the time-consuming step of the process, but would only need to be done once for each database and length for which you want to compare.  To build a catalog for the test example above, the command would be:

```perl
perl run_nmer_match.pl -a build -l 15 -s t/test_db.fasta -c catalogs/test
``` 

This will proceed to break up all protein sequences in the database file (-s) into all possible 15mers (-l) and store it in the directory 'catalogs/test' (-c).  Inside that directory, there will be a JSON file called 'catalog_info.json' with details about the catalog as well as serialized peptide data for later loading in the search step. **NOTE**: The serialization is currently Perl-specific, but could easily be updated to use a different backend (e.g., database, JSON, etc.).

### Searching for matches
Now that you have the database built, you are ready to search for matches.  Currenty, the tool supports a _query_ peptide file that is simply a _list_ (not FASTA) of peptides ([example query list](t/test_mm.lst)).  To search this list for peptides in the database file that match with up to 3 mismatches:

```bash
perl run_nmer_match.pl -a search -c catalogs/test -q t/test_mm.lst -o output.tsv
```

This will search the query file (-q) against the catalog built in the previous step (-c) and output the results to 'output.tsv' (-o).


### Important notes
* At this time, the tool is limited to searching for 1 length at a time.


## Docker container

To run the script from the docker container, you'll need to give it access to the directory where you have your input files with the -v option, e.g.:

```bash
docker run -v $PWD:/scratch -w /scratch IMAGE_ID -a build -l 15 -s test_data/human.fasta -c catalogs/humb1

```

Simply replace the 'perl run\_nmer\_match.pl' from the above commands with:

```bash
docker run -v $PWD:/scratch -w /scratch IMAGE_ID ...
```
Replace IMAGE_ID with the docker image ID (e.g., gitlab.lji.org:4567/jgbaum/nmer-match), or if you've built and tagged the image locally, you would use that ID.

### To validate tests are passing
To validate tests are passing in the container:

```bash
docker run --entrypoint "prove" -it IMAGE_ID /app/nmer-match/t
```

### Container issues

* The Docker container results in ~25% performance hit in terms of time, in limited testing.

## Test datasets

There are two small test datasets:

### Set 1

A semi-systematic test of catching mismatches at different positions in the peptides and query sequences.

 * [t/test\_db.fasta](t/test_db.fasta) - A set of 2 sequences against which to query.
 * [t/test\_mm.fasta](t/test_mm.fasta) - A set of 13 peptides with mismatches to the test_db at varying positions.
 * [t/test\_mm.lst](t/test_mm.lst) - Same content as the fasta, but in list format.

### Set 2

A more systematic test of catching mismatches at different positions in the peptides and query sequences.

 * [t/test\_db.fasta](t/test_db.fasta) - A set of 3 sequences against which to query.
 * [t/test\_mm.fasta](t/test_mm.fasta) - A set of 34 peptides with mismatches to the test_db at varying positions.
 * [t/test\_mm.lst](t/test_mm.lst) - Same content as the fasta, but in list format.

The headers of the query fasta file are in the following format:

```
{dbseq}_p{position}_mm{mismatches}_r{residues}
```

 * dbseq = the sequence identifier of the expected match in the database
 * position = the position in the database sequence of the expected match
 * mm = the number of mismatches to the database sequences
 * residues = the residues in the peptide with mismatches to the database sequence

So, for example, a peptide that should match sequence A at position 0 with 1 mismatch at position 4 of the peptide would have the following identifier:

```
A_p0_mm1_r4
```

For peptides with 2 or more mismatches, each residue that has a mismatch is separated by a '.'.  So, for example, a peptide that should match sequence A at position 0 with 2 mismatches at positions 4 and 6 of the peptide would have the following identifier:

```
A_p0_mm2_r4.6
```

Similarly, if a peptide matches more than 1 database protein the identifiers and positions are separated by '.'.  For example, a peptide that should match sequence A at position 0 and sequence B at position 10 with 1 mismatch at position 4 of the peptide would have the identifier:

```
A.B_p0.10_mm1_r4
```

Currently, there are no test peptides with a different number of mismatches between two sequences nor is there a way to specify this in the identifier.