import argparse
from .preprocessor import Preprocessor
from .matcher import Matcher

def run_preprocessor():
  parser = argparse.ArgumentParser()
  parser.add_argument('-p', '--proteome', required=True)
  parser.add_argument('-n', '--proteome_name', default='')
  parser.add_argument('-k', '--kmer_size', type=int, required=True)
  parser.add_argument('-P', '--preprocessed_files_path', default='.')
  args = parser.parse_args()

  Preprocessor(
    proteome=args.proteome,
    proteome_name=args.proteome_name,
    preprocessed_files_path=args.preprocessed_files_path,
  ).preprocess(k=args.kmer_size)

def run_matcher():
  parser = argparse.ArgumentParser()
  parser.add_argument('-q', '--query', required=True)
  parser.add_argument('-p', '--proteome_file', required=True)
  parser.add_argument('-m', '--max_mismatches', type=int, default=0)
  parser.add_argument('-k', '--kmer_size', type=int, default=5)
  parser.add_argument('-P', '--preprocessed_files_path', default='.')
  parser.add_argument('-b', '--best_match', action='store_true', default=False)
  parser.add_argument('-f', '--output_format', default='csv')
  parser.add_argument('-o', '--output_name', default='')
  parser.add_argument('-v', '--sequence_version', action='store_false', default=True)
  args = parser.parse_args()

  Matcher(
    query=args.query,
    proteome_file=args.proteome_file,
    max_mismatches=args.max_mismatches,
    k=args.kmer_size,
    preprocessed_files_path=args.preprocessed_files_path,
    best_match=args.best_match,
    output_format=args.output_format,
    output_name=args.output_name,
    sequence_version=args.sequence_version,
  ).match()
