from preprocessor import Preprocessor
from matcher import Matcher

from scipy.stats import fisher_exact

import pandas as pd
import math
import os
import glob

class ConservationAnalysis(object):
  '''
  Class that takes data with two columns: peptides and binary values for any given feature
  and produces output from Fisher's exact test for a conservation analysis. Output can be
  either the p-value, odds ratio, 2x2 table, or all of the above.

  data = 
  proteome = path to .fasta file of reference proteome

  For getting homology levels, the user can input EITHER:
  
  max_mismatches = # of mismatches to use a a threshold for conservation levels
  
  or 
  
  homology_threshold = # between 0 and 1 representing up to and including that homology
                       threshold. Recommended to not go much below 0.5 (50%).
  '''
  def __init__(self, data, proteome, homology_threshold):
    if data.split('.')[1] == 'csv':
      df = pd.read_csv(data)
    elif data.split('.')[1] == 'tsv':
      df = pd.read_csv(data, sep='\t')
    elif type(data) == pd.core.frame.DataFrame:
      df = data

    assert df.shape[1] == 2, 'Data received is not two columns.'

    for i in list(df.iloc[:, 0]):
      assert type(i) == str, '1st column contains non-string data.'
    
    for i in list(df.iloc[:, 1]):
      assert i in [0, 1], 'Not all values in 2nd column are binary.'

    if homology_threshold > 1 or 0 > homology_threshold:
      raise ValueError('Homology threshold outside of [0,1] boundary.')
    else:
      # specify the maximum number of mismatches from homology threshold value and
      # shortest sequence length from list of peptides
      min_length = len(min(list(df.iloc[:,0]), key=len))
      self.max_mismatches = math.ceil(min_length - (min_length * homology_threshold))

    self.df = df
    self.peptides = list(df.iloc[:,0])
    self.proteome = proteome
    self.split = math.floor(min_length / (self.max_mismatches + 1))

  def preprocess(self):
    print('Preprocessing proteome...')
    Preprocessor(self.proteome, self.split, 'pickle').preprocess()
    print('Finished preprocessing.')

  def remove_preprocessed_data(self):
    for file in glob.glob('./*.pickle'):
      os.remove(file)
    print('Removed preprocessed files.')
  
  def search(self):
    return Matcher(self.peptides, self.proteome, self.split, output_format='').match()

  def odds_ratio(self):
    return fisher_exact(table)[0]

  def p_value(self, table):
    return fisher_exact(table)[1]
  
  def create_2x2_table(self):
    pass

  def run(self):
    self.preprocess()
    df = self.search()
    self.remove_preprocessed_data()

ConservationAnalysis('test.csv', './proteomes/9606_uniprot_small.fa', 0.8).run()
