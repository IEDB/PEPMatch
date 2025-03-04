#!/usr/bin/env python3

import tempfile
import collections
import pprint
import subprocess
import csv
import os
import pandas as pd
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
  def __init__(
      self, benchmark, query, proteome, lengths, max_mismatches, method_parameters
    ):
    """Initialize the tool to run the benchmarks

    Positional arguments:
    mismatches -- the number of mismatches for which to run the benchmark
    nmer_script_path -- path to 'run_nmer_match.pl'

    method_parameters can include:
    perl_exe -- the perl executable to use (default: perl)
    perl_include_path -- additional paths to add to the perl include path
    catalog_master_dir -- directory in which to store the database catalogs
    output_master_dir -- directory in which to store the output files
    """
    self.benchmark = benchmark
    self.query = query
    self.proteome = proteome

    self.mismatches = max_mismatches

    # this must be defined in method_parameters; the rest of the
    # parameters are optional and have defaults
    self.nmer_script_path = method_parameters['nmer_script_path']

    # if it's a relative directory, we assume it is relative to this script
    if not os.path.isabs(self.nmer_script_path):
      print("Relative path to nmer script provided: " + self.nmer_script_path)
      self.nmer_script_path = os.path.join(
        os.path.dirname(os.path.abspath(__file__)), self.nmer_script_path
      )
      print("Absolute path: " + self.nmer_script_path)

    if 'perl_exe' in method_parameters and method_parameters['perl_exe']:
      self.perl_exe = method_parameters['perl_exe']
    else:
      self.perl_exe = 'perl'

    # set the catalog & output directories
    if 'catalog_master_dir' in method_parameters and \
        method_parameters['catalog_master_dir']:
    
      self.catalog_master_dir = method_parameters['catalog_master_dir']
    
    else:
      self.catalog_master_dir = tempfile.mkdtemp(prefix='catalogs.')

    if 'output_master_dir' in method_parameters and \
        method_parameters['output_master_dir']:
    
      self.output_master_dir = method_parameters['output_master_dir']
    
    else:
      self.output_master_dir = tempfile.mkdtemp(prefix='outputs.')

    # set the extra perl include paths if given
    if 'perl_include_path' in method_parameters and \
        method_parameters['perl_include_path']:
    
      self.perl_include_path = method_parameters['perl_include_path']
    
    else:
      self.perl_include_path = None

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
    return 'NmerMatch (J Greenbaum)'


  def preprocess_proteome(self):
    proteome = self.proteome

    # if lengths is empty, throw an exception
    if (len(self.lengths)) == 0:
      raise Exception(
        'Query must be preprocessed before proteome, so that peptide lengths are known'
      )

    # given a proteome file, process it for each length and store
    # it's path in the self.proteomes dict
    for length in self.lengths:

      # if this catalog already exists, skip creating it
      if proteome in self.catalogs.keys():
        if length in self.catalogs[proteome].keys():
          print("Catalog for this proteome & length already exist...skipping")
          continue

      catalog_name = f'{self.catalog_master_dir}/{str(length)}'
      cmd = f'{self.call_script} {self.nmer_script_path} -a build -l {length} ' \
            f'-c {catalog_name} -s {proteome}'
      
      print('Executing: ' + cmd)
      subprocess.run(cmd, shell=True, check=True)

      self.catalogs[proteome][length] = catalog_name

    return 0


  def preprocess_query(self):
    query = self.query

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
        query_temp_file = tempfile.NamedTemporaryFile(
          prefix=str(seq_len) + '_', mode='w', delete=False
        )
        query_lists[seq_len] = query_temp_file
      fh = query_lists[seq_len]
      fh.write(str(f.seq) + "\n")

    # close the files
    num_ql = len(query_lists)
    print("Created " + str(num_ql) + " input query lists:")
    for q_list in query_lists:
      print(str(q_list) + "\t" + query_lists[q_list].name)
      query_lists[q_list].close()

    self.query_fa2lists[query] = query_lists
    self.lengths = list(query_lists.keys())

    return 0


  def search(self):

    query = self.query
    proteome = self.proteome

    # for each length & corresponding number of mismatches, run the search
    query_files_by_length = self.query_fa2lists[query]
    m = self.mismatches

    # work through the lengths
    for length in self.lengths:
      catalog_name = self.catalogs[proteome][length]
      query_list_file = query_files_by_length[length].name

      out_temp_file = tempfile.NamedTemporaryFile(delete=False)
      output_file = out_temp_file.name

      print(f"Running search for length {length} with {m} mistmatches")
      print("Output file: " + output_file)

      self.result_sets[proteome][query][length][m] = output_file
      cmd = f'{self.call_script} -a search -c {catalog_name} -q {query_list_file} ' \
            f' -m {m} -o {output_file}'

      # if we're searching for the best match, regardless of mismatches we change the
      # action to 'search-deep' and ignore the mismatch parameter
      if (self.mismatches == -1):
        cmd = f'{self.call_script} -a search-deep -c {catalog_name} '\
              f'-q {query_list_file} -o {output_file}'

      print('Executing: ' + cmd)
      subprocess.run(cmd, shell=True, check=True)

    return self.get_results()


  def get_results(self):

    query = self.query
    proteome = self.proteome

    # retrieve the result for a given proteome, query, and max mm
    # and return a list with TSV rows:
    # query_peptide matching_protein match_position matching_peptide num_mm
    results = list()
    m = self.mismatches

    # iterate through results for all length
    for length in self.lengths:
      output_file = self.result_sets[proteome][query][length][m]
      with open(output_file) as tsvin:
        tsvin = csv.reader(tsvin, delimiter='\t')
        # skip the header row
        next(tsvin)
        for row in tsvin:
          # pull out all matching proteins & positions
          all_matching_proteins_positions = row[4].split('; ')
          for mpp in all_matching_proteins_positions:
            prot, pos = mpp.split(':')
            try:
              results.append(
                ",".join([row[0], row[1], prot.split('|')[1], pos])
              )
            except IndexError:
              results.append(
                ",".join([row[0], row[1], prot, pos])
              )
    
    columns = ['Query Sequence', 'Matched Sequence', 'Protein ID', 'Index start']
    results_df = pd.DataFrame([s.split(',') for s in results], columns = columns)
    results_df.to_csv('/home/dan/Desktop/nmermatch_results.csv', sep='\t', index=False)
    results_df['Index start'] = results_df['Index start'].astype(int) + 1
    return results_df


def main():
  mismatches = 1
  method_parameters = {
    'nmer_script_path': 'NmerMatch/bin/run_nmer_match.pl',
    'perl_include_path': '/Users/jgbaum/perl_envs/nmer_match/lib/perl5'
  }
  # lengths are currently ignored, so we set it to empty
  lengths=[]
  proteome_file = 'proteomes/tp.fa'
  query_file = 'queries/test_pep.fa'

  tool = Benchmarker(query_file, proteome_file, lengths, mismatches, method_parameters)
  tool.preprocess_query()
  tool.preprocess_proteome()
  tool.search()

if __name__ == "__main__":
    main()