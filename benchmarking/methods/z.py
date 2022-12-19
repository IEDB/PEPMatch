#!/usr/bin/env python3

from Bio import SeqIO


def parse_fasta(file):
    return SeqIO.parse(file, 'fasta')

class Z(object):
    def __init__(self, query, proteome):
        self.query = query
        self.proteome = proteome

    def getZarr(self, string, z): 
        n = len(string) 

        l, r, k = 0, 0, 0
        for i in range(1, n): 

            if i > r: 
                l, r = i, i 
      
                while r < n and string[r - l] == string[r]: 
                    r += 1
                z[i] = r - l 
                r -= 1
            else:
                k = i - l 

                if z[k] < r - i + 1: 
                    z[i] = z[k] 
                else: 
                    l = i 
                    while r < n and string[r - l] == string[r]: 
                        r += 1
                    z[i] = r - l 
                    r -= 1

    def exact_search(self):
        query = list(parse_fasta(self.query))
        proteome = list(parse_fasta(self.proteome))

        all_matches = []
        for peptide in query:
            for protein in proteome:

                pattern = str(peptide.seq)
                text = str(protein.seq)

                concat = pattern + "$" + text 
                l = len(concat) 
              

                z = [0] * l 
                self.getZarr(concat, z) 
              
                for i in range(l): 

                    if z[i] == len(pattern): 
                        all_matches.append((str(peptide.seq), str(peptide.seq), str(protein.id), 0, i - len(pattern) - 1))

        return all_matches

class Benchmarker(Z):
    def __init__(self, query, proteome, lengths, max_mismatches, method_parameters):
        if max_mismatches > 0:
            raise ValueError(self.__str__() + ' cannot do any mismatching.\n')
        elif max_mismatches == -1:
            raise ValueError(self.__str__() + ' does not have a best match feature.\n')

        self.query = query
        self.proteome = proteome
        self.lengths = lengths
        self.max_mismatches = max_mismatches
        self.method_parameters = method_parameters
        
        super().__init__(query, proteome)

    def __str__(self):
        return 'Z algorithm'

    def preprocess_proteome(self):
        raise TypeError(self.__str__() + ' does not preprocess proteomes.\n')

    def preprocess_query(self):
        raise TypeError(self.__str__() + ' does not preprocess queries.\n')

    def search(self):
        matches = self.exact_search()

        all_matches = []
        for match in matches:
            match = list(match)
            match[2] = match[2].split('|')[1]
            all_matches.append(','.join([str(i) for i in match]))

        return all_matches