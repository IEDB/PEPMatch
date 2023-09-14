#!/usr/bin/env python3

import os
import time
import sys
import pandas as pd
import tracemalloc
import importlib
import json
import glob

from pathlib import Path


# add pepmatch and other methods to sys path for importing
for path in [str(Path(__file__).parent), str(Path(__file__)) + '/methods']:
  sys.path.insert(0, path)


def run_benchmark(benchmark: str, benchmark_memory: bool, include_text_shifting: bool):
  """Run the benchmarking for the specified benchmark. Return the benchmarking
  speeds and memory usage for each method.

  Args:
    benchmark (str): name of the benchmark to run (mhc_ligands, milk, coronavirus,
      neoepitopes)
    benchmark_memory (bool): whether or not to benchmark memory usage.
    include_text_shifting (bool): whether or not to include text shifting methods.
  """  

  print('Benchmarking %s dataset...\n' % benchmark)

  with open('benchmarking_parameters.json', 'r') as file:
    benchmarking_parameters = json.load(file)

  inputs = benchmarking_parameters['datasets'][benchmark]
  methods = benchmarking_parameters['methods']

  if not include_text_shifting: # skip text shifting (horspool, boyer_moore, etc.)
    methods = [x for x in methods if not x['text_shifting']]

  columns = [
    'Name', 'Preprocessing Proteome (s)', 'Preprocessing Query (s)',
    'Searching Time (s)', 'Total Time (s)', 'Memory Usage (MB)', 'Accuracy (%)'
  ]

  benchmark_df = pd.DataFrame(columns = columns)
  for method in methods:
    print('Initializing method...: ' + method['name'] + '\n')
    try:
      if (method['name'] == 'PEPMatch'):
        get_benchmark_object = getattr(
          importlib.import_module('pepmatch.benchmarker'), 'Benchmarker'
        )
      else:
        get_benchmark_object = getattr(
          importlib.import_module(method['name']), 'Benchmarker'
        )
      
      benchmark_tool = get_benchmark_object(
        benchmark,
        Path(__file__).parent / inputs['query'],  
        Path(__file__).parent / inputs['proteome'],
        inputs['lengths'], inputs['mismatches'], 
        method['method_parameters']
      )
    
    except ValueError as error:
      print(error)
      continue

    total_time = 0
    print('Preprocessing query...\n')
    try:
      preprocess_query_start = time.time()
      benchmark_tool.preprocess_query()
      preprocess_query_end = time.time()
      preprocess_query_time = preprocess_query_end - preprocess_query_start
      total_time += preprocess_query_time
    except TypeError as error:
      print(error)
      preprocess_query_time = 'N/A'

    print('Preprocessing proteome...\n')
    try:
      preprocess_proteome_start = time.time()
      benchmark_tool.preprocess_proteome()
      preprocess_proteome_end = time.time()
      preprocess_proteome_time = preprocess_proteome_end - preprocess_proteome_start
      total_time += preprocess_proteome_time
    except TypeError as error:
      print(error)
      preprocess_proteome_time = 'N/A'

    print('Searching...\n')
    search_time_start = time.time()
    results = benchmark_tool.search()
    search_time_end = time.time()
    search_time = search_time_end - search_time_start
    total_time += search_time

    memory_use = 'N/A'
    if benchmark_memory:
      print('Checking memory usage...\n')
      tracemalloc.start()
      benchmark_tool.search()
      _, peak = tracemalloc.get_traced_memory()
      memory_use = peak / (10**6)
      tracemalloc.stop()

    print('Calculating accuracy...\n')
    accuracy_result = accuracy(results, inputs['expected'])

    benchmark_stats = [
      str(benchmark_tool), 
      preprocess_proteome_time, 
      preprocess_query_time, 
      search_time, total_time, memory_use, accuracy_result
    ]

    new_df = pd.DataFrame([benchmark_stats], columns = columns)
    benchmark_df = pd.concat([benchmark_df, new_df], ignore_index = True)

    print('Done benchmarking', str(benchmark_tool), '\n\n')
    print(benchmark_df)

  return benchmark_df


def accuracy(results, expected_file): 
  """Function that calculates the accuracy of your tool from the query
  that is being used.

  Args:
    results: list of results from the search.
    expected_file: file with expected matches.
  """
  expected = []

  # open file with expected matches for checking accuracy
  with open(expected_file, 'r') as f:
    lines = f.readlines()
    for line in lines:
      expected.append(line.replace('\n', ''))


  # return intersection of real and expected matches 
  # divided by number of expected times 100 for percentage
  
  return len(set(results).intersection(set(expected))) / (len(expected)) * 100


def main():
  import argparse
  parser = argparse.ArgumentParser()

  parser.add_argument(
    '-b', '--benchmark',
    choices=['mhc_ligands', 'milk', 'coronavirus', 'neoepitopes'],
    required=True
  )
  parser.add_argument(
    '-m', '--benchmark_memory', action='store_true', default=False
  )
  parser.add_argument(
    '-t', '--include_text_shifting', action='store_true', default=False
  )

  args = parser.parse_args()
  benchmark = args.benchmark
  benchmark_memory = args.benchmark_memory
  include_text_shifting = args.include_text_shifting
  
  master_df = run_benchmark(
    benchmark, benchmark_memory, include_text_shifting
  )
  master_df['Searching Time (s)'] = pd.to_numeric(master_df['Searching Time (s)'])
  master_df.round(3).to_csv('%s_benchmarking.tsv' % benchmark, sep='\t', index=False)

  # remove files after benchmarking
  for file in glob.glob('*.db') + glob.glob('*.pkl'):
    os.remove(file)

if __name__ == '__main__':
    main()