#!/usr/bin/env python3

import os
import argparse
import time
import pandas as pd
import tracemalloc
import importlib
import json

##########################################################################

# global variables for generalization
benchmark_columns = ['Name', 'Preprocessing Proteome (s)', 'Preprocessing Query (s)',
                     'Searching Time (s)', 'Total Time (s)', 'Memory Usage (MB)', 'Accuracy (%)']

directory = os.path.dirname(os.path.abspath(__file__))

valid_datasets = ['mhc_ligands', 'milk', 'coronavirus', 'neoepitopes']

def parse_arguments():
    # add arguments
    parser = argparse.ArgumentParser()

    # select the benchmarking you want to use (mhc_ligands, milk, coronavirus, neoepitopes)
    parser.add_argument('-d', '--dataset', nargs=1, type=str, default = ['mhc_ligands'])

    # skip_memory_benchmark - a Boolean for whether to skip the memory benchmark
    parser.add_argument('-s', '--skip_mem', action='store_true')

    # include text-shifting algorithms for benchmarking
    parser.add_argument('-t', '--text_shifting', action='store_true')

    arguments = parser.parse_args()

    dataset = arguments.dataset[0]
    skip_mem = arguments.skip_mem
    text_shifting = arguments.text_shifting

    if dataset not in valid_datasets:
        raise ValueError('Invalid dataset. Please pass "mhc_ligands", "milk", "coronavirus", or "neoepitopes".')

    benchmark_options = [dataset, skip_mem, text_shifting]

    return benchmark_options


def accuracy(results, expected_file): 
    '''
    Function that calculates the accuracy of your tool from the query
    that is being used.
    '''
    expected = []

    # open file with expected matches for checking accuracy
    with open(expected_file, 'r') as f:
        lines = f.readlines()
        for line in lines:
            expected.append(line.replace('\n', ''))


    # return intersection of real and expected matches divided by number of expected 
    # times 100 for percentage

    return len(set(results).intersection(set(expected))) / (len(expected)) * 100


def benchmark_algorithms(benchmark_options):
    dataset = benchmark_options[0]
    skip_mem = benchmark_options[1]
    include_text_shifting = benchmark_options[2]

    with open('benchmarking_parameters.json', 'r') as file:
        benchmarking_parameters = json.load(file)

    inputs = benchmarking_parameters['datasets'][dataset]
    algorithms = benchmarking_parameters['algorithms']

    if not include_text_shifting:
        algorithms = [x for x in algorithms if not x['text_shifting']]

    benchmark_df = pd.DataFrame(columns = benchmark_columns)
    for algorithm in algorithms:
        print('Initializing algorithm...\n')
        try:
            get_benchmark_object = getattr(importlib.import_module(algorithm['name']), 'Benchmarker')
            benchmark_tool = get_benchmark_object(
                os.path.join(directory, inputs['query']),  
                os.path.join(directory, inputs['proteome']), 
                inputs['lengths'], inputs['mismatches'], 
                algorithm['algorithm_parameters'])
        except ValueError as error:
            print(error)
            continue

        print('Benchmarking with', str(benchmark_tool))

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
        if not skip_mem:
            print('Checking memory usage...\n')
            tracemalloc.start()
            benchmark_tool.search()
            current, peak = tracemalloc.get_traced_memory()
            memory_use = peak / (10**6)
            tracemalloc.stop()

        print('Calculating accuracy...\n')
        accuracy_result = accuracy(results, inputs['expected'])
        
        benchmarks = [str(benchmark_tool), preprocess_proteome_time, preprocess_query_time, search_time, total_time, memory_use, accuracy_result]

        new_df = pd.DataFrame([benchmarks], columns = benchmark_columns)
        benchmark_df = pd.concat([benchmark_df, new_df], ignore_index = True)

        print('Done benchmarking', str(benchmark_tool), '\n\n')

    return benchmark_df



def main():
    benchmark_options = parse_arguments()
    
    master_df = benchmark_algorithms(benchmark_options)
    master_df['Searching Time (s)'] = pd.to_numeric(master_df['Searching Time (s)'])

    print(master_df.round(3))
    master_df.round(3).to_excel('benchmarking.xlsx')

if __name__ == '__main__':
    main()