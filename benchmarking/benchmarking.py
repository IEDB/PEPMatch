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


# add methods to sys path for importing
methods_dir = str(Path(__file__).parent / 'methods')
if methods_dir not in sys.path:
  sys.path.insert(0, methods_dir)


def run_benchmark(
  benchmark: str, benchmark_memory: bool, include_text_shifting: bool
) -> pd.DataFrame:
  """Run the benchmarking for the specified benchmark. Return the benchmarking
  speeds and memory usage for each method.

  Args:
    benchmark (str): name of the benchmark to run (mhc_ligands, milk, coronavirus,
      neoepitopes)
    benchmark_memory (bool): whether or not to benchmark memory usage.
    include_text_shifting (bool): whether or not to include text shifting methods."""  

  print('Benchmarking %s dataset...\n' % benchmark)

  with open('benchmarking_parameters.json', 'r') as file:
    benchmarking_parameters = json.load(file)

  inputs = benchmarking_parameters['datasets'][benchmark]
  methods = benchmarking_parameters['methods']

  if not include_text_shifting: # skip text shifting (horspool, boyer_moore, etc.)
    methods = [x for x in methods if not x['text_shifting']]

  columns = [
    'Name', 'Preprocessing Proteome (s)', 'Preprocessing Query (s)',
    'Searching Time (s)', 'Total Time (s)', 'Memory Usage (MB)', 'Recall (%)'
  ]

  benchmark_df = pd.DataFrame(columns = columns)
  for method in methods:
    print('Initializing method...: ' + method['name'] + '\n')
    try:
      if (method['name'] == 'PEPMatch'):
        benchmark_object = getattr(
          importlib.import_module('pepmatch.benchmarker'), 'Benchmarker'
        )
      else:
        benchmark_object = getattr(
          importlib.import_module(method['name']), 'Benchmarker'
        )

      benchmark_tool = benchmark_object(
        benchmark=benchmark,
        query=Path(__file__).parent / inputs['query'],
        proteome=Path(__file__).parent / inputs['proteome'],
        lengths=inputs['lengths'], 
        max_mismatches=inputs['mismatches'],
        method_parameters=method['method_parameters']
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
    results_df = benchmark_tool.search()
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

    print('Calculating recall...\n')
    expected_df = pd.read_csv(inputs['expected'], sep='\t')
    recall_result = recall(results_df, expected_df)

    benchmark_stats = [
      str(benchmark_tool),
      preprocess_proteome_time,
      preprocess_query_time,
      search_time,
      total_time,
      memory_use,
      recall_result
    ]

    new_df = pd.DataFrame([benchmark_stats], columns = columns)
    benchmark_df = pd.concat([benchmark_df, new_df], ignore_index = True)

    print('Done benchmarking', str(benchmark_tool), '\n\n')
    print(benchmark_df)

  return benchmark_df


def recall(results_df: pd.DataFrame, expected_df: pd.DataFrame) -> float:
  """Function that calculates the recall of your tool from the query
  that is being used.

  Args:
    results: pandas dataframe with results from the benchmarking.
    expected_df: pandas dataframe with expected matches for the benchmarking."""

  columns = ['Query Sequence', 'Matched Sequence', 'Protein ID', 'Index start']
  results = results_df[columns].drop_duplicates(subset=columns)
  expected = expected_df[columns].drop_duplicates(subset=columns)

  results['Index start'] = results['Index start'].fillna(0).astype(int)
  
  matched_rows = pd.merge(results, expected, how='inner', on=columns)
  matched_rows = matched_rows.drop_duplicates(subset=columns)

  # calculate the recall
  total_expected = len(expected)
  total_matched = len(matched_rows)

  recall = (total_matched / total_expected) * 100

  return min(recall, 100)


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