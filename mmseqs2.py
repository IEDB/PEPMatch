#!/usr/bin/env python3

import os, glob

directory = os.path.dirname(os.path.abspath(__file__))

class MMseqs2(object):
    def __init__(self, query, proteome, max_mismatches, algorithm_parameters):
        # if max_mismatches == -1:
        #     raise ValueError(self.__str__() + ' does not have a best match feature.\n')

        self.query = query
        self.proteome = proteome
        self.proteome_name = proteome.replace('.fa', '')

        self.max_mismatches = max_mismatches

        bin_directory = algorithm_parameters['bin_directory']
        self.bin_file = os.path.join(bin_directory, 'mmseqs')

    def __str__(self):
        return 'MMseqs2'
    
    def preprocesss(self):
        os.system(self.bin_file + ' createdb ' + self.proteome + ' ' + self.proteome_name)
        os.system(self.bin_file + ' createindex ' + self.proteome_name + ' tmp')
    
    def mmseqs_search(self, query, proteome):

        os.system(self.bin_file + ' easy-search ' + self.query + ' ' + self.proteome_name + ' alnRes.m8' + 
                  ' tmp' + ' -s 7.0' + ' -e 10000' + ' --format-output "qseq,taln,theader,mismatch,tstart"')

        all_matches = []
        with open('alnRes.m8', 'r') as file:
            lines = file.readlines()

            for line in lines:
                match = []
                result = line.split('\t')
                for i in range(len(result)):
                    if i == 2:
                        match.append(result[i].split(' ')[0])
                    elif i == 4:
                        match.append(int(result[i].replace('\n', '')) - 1)
                    else:
                        match.append(result[i])

                all_matches.append(match)

        return all_matches


class Benchmarker(MMseqs2):
    def __init__(self, query, proteome, lengths, max_mismatches, algorithm_parameters):
        self.query = query
        self.proteome = proteome
        self.lengths = lengths
        self.max_mismatches = max_mismatches
        self.algorithm_parameters = algorithm_parameters
        
        super().__init__(query, proteome, max_mismatches, algorithm_parameters)

    def __str__(self):
        return 'MMseqs2'

    def preprocess_proteome(self):
        return self.preprocesss()

    def preprocess_query(self):
        raise TypeError(self.__str__() + ' does not preprocess queries.\n')

    def search(self):
        matches = self.mmseqs_search(self.query, self.proteome)

        all_matches = []
        for match in matches:
            match_string = ''
            for i in match:
                if i == match[-1]:
                    match_string += str(i)
                else:
                    match_string += str(i) + ','

            all_matches.append(match_string)

        # for extension in ['source', 'dbtype', 'index', 'idx', 'lookup', 'pot', 'pto']:
        #     os.remove(glob.glob(os.path.dirname(self.proteome) + '/*.' + extension)[0])

        os.remove('alnRes.m8')

        return all_matches