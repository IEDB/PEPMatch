from pepmatch import Preprocessor
from pepmatch import Matcher

from scipy.stats import fisher_exact

import pandas as pd
import numpy as np
import math
import os
import glob


class ConservationAnalysis(object):
  '''
  Class that takes data with two columns: peptides and binary values for any given feature
  and produces output from Fisher's exact test for a conservation analysis. Output can be
  either the p-value, odds ratio, 2x2 table, or all of the above.
  data = two column dataset (from .csv, .tsv or pandas DataFrame) with peptides in one column
         and binary values in the other column separating the peptides into two groups.
  proteome = path to .fasta file of reference proteome
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
    else:
      raise ValueError('Data passed not supported. Please pass .csv file, .tsv file, or a pandas dataframe of peptides and binary values.')

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
      self.homology_thresholds = [(min_length - i) / min_length for i in range(self.max_mismatches)]

    self.df = df
    self.query_name = os.path.basename(data).split('.')[0] + '_to_' + os.path.basename(proteome).split('.')[0]
    self.group_name = df.columns[-1]
    self.binary_map = dict(zip(df.iloc[:, 0], df.iloc[:, 1]))
    self.peptides = list(df.iloc[:,0])
    self.proteome = proteome
    self.max_homology_threshold = homology_threshold
    self.split = max(2, math.floor(min_length / (self.max_mismatches + 1)))

  def preprocess(self):
    '''Creates preprocessed files needed for PEPMatch search.'''
    Preprocessor(self.proteome, self.split, 'pickle').preprocess()

  def remove_preprocessed_data(self):
    '''Removes preprocessed files used by PEPMatch.'''
    for file in glob.glob(os.path.dirname(self.proteome) + '/*.pickle'):
      os.remove(file)
    print('Removed preprocessed files.')
  
  def search(self):
    '''Searches peptides in proteome with PEPMatch.'''
    return Matcher(self.peptides, self.proteome, self.split, output_format='dataframe').match()

  def threshold_map(self, df):
    '''
    Create dictionary with mismatches as keys and peptides that have matches with that 
    mismatch threshold found within to use to create binary dataframe.
    '''
    threshold_map = {}
    for i in range(self.max_mismatches + 1):
      threshold_map[i] = list(df[df['Mismatches'] <= i]['Peptide Sequence'])
    return threshold_map

  def create_binary_df(self, df):
    '''
    Takes the PEPMatch results and creates a dataframe where the columns represent
    each mismatch threshold as well as the original binary grouping. Binary values
    are placed in each column for each peptide if the peptide is found at that
    mismatch threshold.
    Example output:
          Peptides | 0 | 1 | 2 | Group |
          ------------------------------
          KLQNLNIFL  1   1   1     1
          LEDEERVVRL 0   0   1     0
          RLLDEWFTL  0   1   1     1
    '''
    # replace empty strings with NaN values for conditional below
    df['Mismatches'] = (df['Mismatches'].replace(r'^\s*$', np.NaN, regex=True))

    # 
    threshold_map = self.threshold_map(df)

    data = []
    for peptide in self.peptides:
      values = [peptide]
      for i in range(self.max_mismatches + 1):
        if peptide in threshold_map[i]:
          values.append(1)
        else:
          values.append(0)

      data.append(values)

    binary_df =  pd.DataFrame(data=data, 
                              columns = ['Peptides'] + [str(i) for i in range(self.max_mismatches+1)])
    binary_df[self.df.columns[1]] = binary_df['Peptides'].map(self.binary_map)

    return binary_df

  def create_2x2_table(self, binary_df, threshold):
    conserved_g1     = sum(binary_df[binary_df.iloc[:, -1] == 1][str(threshold)].isin([1]))
    conserved_g2     = sum(binary_df[binary_df.iloc[:, -1] == 0][str(threshold)].isin([1]))
    not_conserved_g1 = sum(binary_df[binary_df.iloc[:, -1] == 1][str(threshold)].isin([0]))
    not_conserved_g2 = sum(binary_df[binary_df.iloc[:, -1] == 0][str(threshold)].isin([0]))

    table = [ [conserved_g1, not_conserved_g1], 
              [conserved_g2, not_conserved_g2] ]

    return table

  def output_2x2_table(self, table, p_value, odds_ratio, homology_threshold, path):
    fout = open(path + '.csv', 'a')
    pd.DataFrame(table, columns = ['Conserved', 'Not Conserved'],
                        index   = [self.group_name, 'Not ' + self.group_name]).to_csv(fout)
    fout.write('Homology threshold: ' + str(round(homology_threshold * 100, 2)) + '%')
    fout.write(' p-value: ' + str(round(p_value, 4)))
    fout.write(' Odds ratio: ' + str(round(odds_ratio, 4)))
    fout.close()

  def p_value(self, table):
    return fisher_exact(table)[1]

  def odds_ratio(self, table):
    return fisher_exact(table)[0]

  def graph_p_values(self):
    pass

  def graph_odds_ratios(self):
    pass

  def run(self):
    print('Running conservation analysis with', len(self.peptides),
          'peptides at >=', str(self.max_homology_threshold * 100) + '%',
          'max homology threshold. Evaluating at', len(self.homology_thresholds), 'thresholds.\n')

    print('Preprocessing proteome...\n')
    self.preprocess()
    print('Finished preprocessing.\n')

    print('Searching peptides in proteome...\n')
    df = self.search()
    print('Finished searching peptides.\n')

    print('Creating binary data for analysis...\n')
    binary_df = self.create_binary_df(df)
    print('Done.\n')

    p_values = []
    odds_ratios = []

    for i in range(self.max_mismatches):
      print("Evaluating Fisher's exact test at", 
             round(self.homology_thresholds[i] * 100, 2), '% homology threshold...')
      
      # create 2x2 table and get Fisher's p-value and odds ratio
      table = self.create_2x2_table(binary_df, i)
      p_value = self.p_value(table)
      odds_ratio = self.odds_ratio(table)

      # write 2x2 table to .csv file
      self.output_2x2_table(table, p_value, odds_ratio, self.homology_thresholds[i],
                            self.query_name + '_' + str(round(self.homology_thresholds[i], 2)))

      # save p-values and odds ratios for charts
      p_values.append(p_value)
      odds_ratios.append(odds_ratio)
      
      print('p-value: ', round(p_value, 4))
      print('Odds ratio: ', round(odds_ratio, 4), '\n')

    ############################
    # TODO: 
    print('Creating charts...\n')
    self.graph_p_values()
    self.graph_odds_ratios()
    print('Done.\n')
    ############################
    
    print('Removing preprocessed data...\n')
    self.remove_preprocessed_data()
    print('Done with analysis.\n')


ConservationAnalysis('/home/dan/projects/pertussis/one_donor_reactive.csv', 
                     '/home/dan/projects/pertussis/D422.fasta', 0.5).run()
