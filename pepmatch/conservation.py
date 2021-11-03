from .preprocessor import Preprocessor
from .matcher import Matcher

from scipy.stats import fisher_exact

import pandas as pd
import os
import glob

class ConservationAnalysis(object):
  '''
  Class that takes data with two columns: peptides and binary values for any given feature
  and produces output from Fisher's exact test for a conservation analysis. Output can be
  either the p-value, odds ratio, 2x2 table, or all of the above.
  '''
  def __init__(self, data, proteome, max_mismatches=-1, homology_threshold=-1):
    if data.split('.')[1] == 'csv':
      df = pd.read_csv(data)
    elif data.split('.')[1] == 'tsv':
      df = pd.read_csv(data, sep='\t')
    elif type(data) == pd.core.frame.DataFrame:
      df = data

    assert df.shape()[1] == 2, 'Data received is not two columns.'
    assert df.iloc[:, 1].dtype == bool, '2nd column not binary values.'
    for i in list(df.iloc[:, 0]):
          assert type(i) == str, '1st column contains non-string data.'
    
    if max_mismatches == -1 and homology_threshold == -1:
      raise ValueError('Neither mismatch or homology threshold specified.')

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
