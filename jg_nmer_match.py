#!/usr/bin/env python3

import tempfile
import collections
import pprint
import subprocess
import csv
from Bio import SeqIO

pp = pprint.PrettyPrinter(indent=4)

# define a function to allow nested defaultdicts
def rec_dd():
    return collections.defaultdict(rec_dd)

class Benchmarker(object):
    # if the required perl modules are installed in path that is not
    # in your PERL5LIB environment variable, that should be passed
    # as the perl_include_path here. E.g., if libraries are installed
    # to 'mylibs', the argument to be passed would be mylibs/lib/perl5
    def __init__(self, lengths, mismatches, nmer_script_path,
                 perl_exe='perl', perl_include_path=None,
                 catalog_master_dir=None,
                 output_master_dir=None):
        """Initialize the tool to run the benchmarks


        Positional arguments:
        mismatches -- the number of mismatches for which to run the benchmark
        nmer_script_path -- path to 'run_nmer_match.pl'

        Keyword arguments:
        perl_exe -- the perl executable to use (default: perl)
        perl_include_path -- additional paths to add to the perl include path
        catalog_master_dir -- directory in which to store the database catalogs
        output_master_dir -- directory in which to store the output files
        """
        self.mismatches = mismatches

        self.nmer_script_path = nmer_script_path

        # define the paths
        self.perl_exe = perl_exe

        # set the catalog & output directories
        if (catalog_master_dir):
            self.catalog_master_dir = catalog_master_dir
        else:
            self.catalog_master_dir = tempfile.mkdtemp(prefix='catalogs.')

        if (output_master_dir):
            self.output_master_dir = output_master_dir
        else:
            self.output_master_dir = tempfile.mkdtemp(prefix='outputs.')

        # set the extra perl include paths if given
        self.perl_include_path = perl_include_path

        # definte the system call to the script
        self.call_script = self.perl_exe
        if (self.perl_include_path):
            self.call_script += ' -I' + self.perl_include_path
        self.call_script += ' ' + self.nmer_script_path

        # set up the catlogs to be a nested dict keyed by proteome name
        # and length, pointing to the location of the catalog
        self.catalogs = collections.defaultdict(dict)

        # set up the results nested dictionary to track the output by:
        # proteome_file => nmer_length => number_of_mismatches
        self.result_sets = rec_dd()

        # initialize the dict that points to the query peptide lists
        # indexed first by query file, then by length
        self.query_fa2lists = collections.defaultdict(dict)

        # initilaize the different peptide lengths for which we need to search
        self.lengths = list()


    def __str__(self):
        return 'nmer_match - J. Greenbaum'

    def preprocess_proteome(self, proteome):

        # if lengths is empty, throw an exception
        if (len(self.lengths)) == 0:
            raise Exception('Query must be preprocessed before proteome, so that peptide lengths are known')

        # given a proteome file, process it for each length and store
        # it's path in the self.proteomes dict
        for l in self.lengths:

            # if this catalog already exists, skip creating it
            if proteome in self.catalogs.keys():
                if l in self.catalogs[proteome].keys():
                    print("Catalog for this proteome & length already exist...skipping")
                    continue

            catalog_name = self.catalog_master_dir + '/' + str(l)

            cmd = self.call_script + ' ' + self.nmer_script_path + ' -a build -l ' + str(l) + ' -c ' + catalog_name + ' -s ' + proteome
            print('Executing: ' + cmd)
            subprocess.run(cmd, shell=True, check=True)

            self.catalogs[proteome][l] = catalog_name

        return 0

    def preprocess_query(self, query):
        # given a query fasta file, preprocess it so that it is a peptide list
        # suitable for input to the tool

        fasta_sequences = SeqIO.parse(open(query),'fasta')

        # write the different length sequences to a different file

        # indexed by length, this will point to the filehandle
        query_lists = dict()
        for f in fasta_sequences:
            seq_len = len(f.seq)
            if seq_len not in query_lists.keys():
                # create the new filehandle
                query_temp_file = tempfile.NamedTemporaryFile(prefix=str(seq_len) + '_', mode='w', delete=False)
                query_lists[seq_len] = query_temp_file
            fh = query_lists[seq_len]
            fh.write(str(f.seq) + "\n")

        # close the files
        num_ql = len(query_lists)
        print("Created " + str(num_ql) + " input query lists:")
        for l in query_lists:
            print(str(l) + "\t" + query_lists[l].name)
            query_lists[l].close()

        self.query_fa2lists[query] = query_lists
        self.lengths = list(query_lists.keys())

        return 0

    def search(self, query, proteome):
        # for each length & corresponding number of mismatches, run the search
        query_files_by_length = self.query_fa2lists[query]
        m = self.mismatches

        # work through the lengths
        for l in self.lengths:

            catalog_name = self.catalogs[proteome][l]
            query_list_file = query_files_by_length[l].name

            out_temp_file = tempfile.NamedTemporaryFile(delete=False)
            output_file = out_temp_file.name

            print("Running search for length " + str(l) + " with " + str(m) + " mistmatches")
            print("Output file: " + output_file)

            self.result_sets[proteome][query][l][m] = output_file
            cmd = self.call_script + ' -a search -c ' + catalog_name + ' -q ' + query_list_file + ' -m ' + str(m) + ' -o ' + output_file
            print('Executing: ' + cmd)
            subprocess.run(cmd, shell=True, check=True)

        return self.get_results(query, proteome)


    def get_results(self, query, proteome):
        # retrieve the result for a given proteome, query, and max mm
        # and return a list with TSV rows:
        # query_peptide matching_protein match_position matching_peptide num_mm
        results = list()
        m = self.mismatches

        # iterate through results for all length
        for l in self.lengths:
            output_file = self.result_sets[proteome][query][l][m]
            with open(output_file) as tsvin:
                tsvin = csv.reader(tsvin, delimiter='\t')
                # skip the header row
                next(tsvin)
                for row in tsvin:
                    # pull out all matching proteins & positions
                    all_matching_proteins_positions = row[4].split('; ')
                    for mpp in all_matching_proteins_positions:
                        prot, pos = mpp.split(':')
                        #results.append("\t".join([row[0], prot, pos, row[1], row[2]]))
                        results.append("\t".join([row[0], prot, pos, row[2]]))

        return results




def main():
    mismatches = 1
    nmer_match_path = '/Users/jgbaum/projects/nmer-match/bin/run_nmer_match.pl'
    proteome_file = 'proteomes/tp.fa'
    query_file = 'queries/test_pep.fa'

    tool = benchmark_jg_nmer_match(mismatches, nmer_match_path)
    tool.preprocess_proteome(proteome_file)
    tool.preprocess_query(query_file)
    results = tool.search(query_file, proteome_file)

if __name__ == "__main__":
    main()