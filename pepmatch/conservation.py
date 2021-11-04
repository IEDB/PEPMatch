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
      self.df = pd.read_csv(data)
    elif data.split('.')[1] == 'tsv':
      self.df = pd.read_csv(data, sep='\t')
    elif type(data) == pd.core.frame.DataFrame:
      self.df = data

    assert df.shape()[1] == 2, 'Data received is not two columns.'
    assert df.iloc[:, 1].dtype == bool, '2nd column not binary values.'
    for i in list(df.iloc[:, 0]):
          assert type(i) == str, '1st column contains non-string data.'
    
    if homology_threshold > 1 or 0 > homology_threshold:
      raise ValueError('Homology threshold outside of [0,1] boundary.')
    else:
      # specify the maximum number of mismatches from homology threshold value and
      # shortest sequence length from list of peptides
      min_length = min(list(df.iloc[:,0]), key=len)
      self.max_mismatches = math.ceil(min_length - (min_length / homology_threshold))

    self.proteome = proteome

  def preprocess(self, proteome, split):
    print('Preprocessing proteome...')
    Preprocessor(self.proteome, self.split, 'pickle').preprocess()
    print('Finished preprocessing.')

  def remove_preprocessed_data(self, proteome):
    for file in glob.glob(os.path.dirname(self.proteome) + '/*.pickle'):
      os.remove(file)

  def odds_ratio(self):
    return fisher_exact(table)[0]

  def p_value(self, table):
    return fisher_exact(table)[1]
  
  def create_2x2_table(self):
    pass
