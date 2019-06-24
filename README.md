# Nmer Match

## Docker container

To run the script from the docker container, you'll need to give it access to the directory where you have your input files with the -v option, e.g.:
```bash
docker run -v $PWD:/scratch -w /scratch IMAGE_ID -a build -l 15 -s test_data/human.fasta -c catalogs/humb1
```

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