#!/usr/bin/env python3

import os, glob, shutil

directory = os.path.dirname(os.path.abspath(__file__))

class MMseqs2(object):
    def __init__(self, query, proteome, max_mismatches, method_parameters):
        if max_mismatches == -1:
            raise ValueError(self.__str__() + ' does not have a best match feature.\n')

        self.query = query
        self.proteome = proteome
        self.proteome_name = proteome.replace('.fa', '')

        self.max_mismatches = max_mismatches

        bin_directory = method_parameters['bin_directory']
        self.bin_file = os.path.join(bin_directory, 'mmseqs')

    def __str__(self):
        return 'MMseqs2'
    
    def preprocesss(self):
        os.system(self.bin_file + ' createdb ' + self.proteome + ' ' + self.proteome_name)
        os.system(self.bin_file + ' createindex ' + self.proteome_name + ' tmp')
    
    def mmseqs_search(self, query, proteome):

        os.system(self.bin_file + ' easy-search ' + self.query + ' ' + self.proteome_name + ' results.m8' + 
                  ' tmp' + ' -s 7.0' + ' -e 10000' + ' --format-output "qseq,taln,theader,mismatch,tstart"')

        all_matches = []
        with open('results.m8', 'r') as file:
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
    def __init__(self, query, proteome, lengths, max_mismatches, method_parameters):
        self.query = query
        self.proteome = proteome
        self.lengths = lengths
        self.max_mismatches = max_mismatches
        self.method_parameters = method_parameters
        
        super().__init__(query, proteome, max_mismatches, method_parameters)

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

        # remove all files generated by MMseqs2
        for extension in ['source', 'dbtype', 'index', 'idx', 'lookup']:
            for file in glob.glob(os.path.dirname(self.proteome) + '/*.' + extension):
                os.remove(file)
        os.remove(self.proteome.replace('.fasta', 'sta'))
        os.remove(self.proteome.replace('.fasta', 'sta_h'))
        os.remove('results.m8')

        # remove tmp folder
        shutil.rmtree('tmp')

        return all_matches