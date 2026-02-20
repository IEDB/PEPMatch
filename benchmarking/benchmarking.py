import argparse
import importlib
import json
import sys
import time
import tracemalloc
from pathlib import Path

import pandas as pd

METHODS_DIR = str(Path(__file__).parent / 'methods')
if METHODS_DIR not in sys.path:
  sys.path.insert(0, METHODS_DIR)

RESULT_COLUMNS = ['Query Sequence', 'Matched Sequence', 'Protein ID', 'Index start']


def load_config():
  with open(Path(__file__).parent / 'benchmarking_parameters.json') as f:
    return json.load(f)


def load_method(name, benchmark, dataset, method_params):
  try:
    if name == 'PEPMatch':
      module = importlib.import_module('pepmatch.benchmarker')
    else:
      module = importlib.import_module(name)

    return module.Benchmarker(
      benchmark=benchmark,
      query=Path(__file__).parent / dataset['query'],
      proteome=Path(__file__).parent / dataset['proteome'],
      lengths=dataset['lengths'],
      max_mismatches=dataset['mismatches'],
      method_parameters=method_params,
    )
  except ValueError as e:
    print(f'  Skipping {name}: {e}')
    return None


def time_step(fn):
  try:
    start = time.perf_counter()
    result = fn()
    return time.perf_counter() - start, result
  except TypeError:
    return None, None


def recall(results_df, expected_df):
  results = results_df[RESULT_COLUMNS].drop_duplicates()
  expected = expected_df[RESULT_COLUMNS].drop_duplicates()

  results['Index start'] = results['Index start'].fillna(0).astype(int)
  expected['Index start'] = expected['Index start'].fillna(0).astype(int)

  matched = pd.merge(results, expected, how='inner', on=RESULT_COLUMNS).drop_duplicates()
  return min((len(matched) / len(expected)) * 100, 100)


def run_benchmark(benchmark, include_memory=False, include_text_shifting=False):
  config = load_config()
  dataset = config['datasets'][benchmark]
  methods = config['methods']

  if not include_text_shifting:
    methods = [m for m in methods if not m['text_shifting']]

  expected_df = pd.read_csv(
    Path(__file__).parent / dataset['expected'], sep='\t'
  )

  rows = []

  for method in methods:
    name = method['name']
    print(f'\n{"=" * 60}')
    print(f'  {name}')
    print(f'{"=" * 60}')

    tool = load_method(name, benchmark, dataset, method['method_parameters'])
    if tool is None:
      continue

    # preprocess proteome
    print('  Preprocessing proteome...')
    preprocess_proteome_time, _ = time_step(tool.preprocess_proteome)
    if preprocess_proteome_time is not None:
      print(f'  -> {preprocess_proteome_time:.3f}s')
    else:
      print(f'  -> N/A')

    # preprocess query
    print('  Preprocessing query...')
    preprocess_query_time, _ = time_step(tool.preprocess_query)
    if preprocess_query_time is not None:
      print(f'  -> {preprocess_query_time:.3f}s')
    else:
      print(f'  -> N/A')

    # search
    print('  Searching...')
    search_time, results_df = time_step(tool.search)
    print(f'  -> {search_time:.3f}s')

    # total
    total = sum(t for t in [preprocess_proteome_time, preprocess_query_time, search_time] if t is not None)

    # memory
    memory = None
    if include_memory:
      print('  Measuring memory...')
      tracemalloc.start()
      tool.search()
      _, peak = tracemalloc.get_traced_memory()
      tracemalloc.stop()
      memory = peak / 1e6
      print(f'  -> {memory:.1f} MB')

    # recall
    recall_pct = recall(results_df, expected_df)
    print(f'  Recall: {recall_pct:.1f}%')

    # cleanup
    if hasattr(tool, 'cleanup'):
      tool.cleanup()

    rows.append({
      'Method': str(tool),
      'Proteome Preprocessing (s)': f'{preprocess_proteome_time:.3f}' if preprocess_proteome_time is not None else 'N/A',
      'Query Preprocessing (s)': f'{preprocess_query_time:.3f}' if preprocess_query_time is not None else 'N/A',
      'Searching (s)': f'{search_time:.3f}',
      'Total (s)': f'{total:.3f}',
      'Memory (MB)': f'{memory:.1f}' if memory is not None else 'N/A',
      'Recall (%)': f'{recall_pct:.1f}',
    })

  results = pd.DataFrame(rows)

  print(f'\n\n{"=" * 80}')
  print(f'  {benchmark.upper()} RESULTS')
  print(f'{"=" * 80}\n')
  print(results.to_string(index=False))

  output_path = f'{benchmark}_benchmarking.tsv'
  results.to_csv(output_path, sep='\t', index=False)
  print(f'\nSaved to {output_path}')

  return results


def main():
  parser = argparse.ArgumentParser(description='PEPMatch Benchmarking Framework')
  parser.add_argument(
    '-b', '--benchmark',
    choices=['mhc_ligands', 'milk', 'coronavirus', 'neoepitopes'],
    required=True,
  )
  parser.add_argument('-m', '--memory', action='store_true', default=False)
  parser.add_argument('-t', '--text_shifting', action='store_true', default=False)
  args = parser.parse_args()

  run_benchmark(args.benchmark, args.memory, args.text_shifting)


if __name__ == '__main__':
  main()
