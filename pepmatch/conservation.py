from .preprocessor import Preprocessor
from .matcher import Matcher

from scipy.stats import fisher_exact

import pandas as pd
import os
import glob

class ConservationAnalysis(object):
  '''
  Class that takes data as a 2 column pandas dataframe
  '''
  def __init__(self, data, proteome, max_mismatches=-1, homology_threshold=-1):
    df = pd.read_csv(data)
    assert df.shape()[1] == 2, 'Data received is not two columns.'

    if max_mismatches == -1 and homology_threshold == -1:
      raise ValueError('Neither mismatch or homology threshold specified.')

  def preprocess(self, proteome, split):
    print('Preprocessing proteome...')
    Preprocessor(proteome, split, 'pickle').preprocess()
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
