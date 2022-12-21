#!/usr/bin/env python3

from Bio import SeqIO
import pandas as pd


def parse_fasta(file):
    return SeqIO.parse(file, 'fasta')

class Horspool(object):
    def __init__(self, query, proteome):
        self.query = query
        self.proteome = proteome

    def exact_search(self):
        query = list(parse_fasta(self.query))
        proteome = list(parse_fasta(self.proteome))

        all_matches = []
        for peptide in query:
            for protein in proteome:

                p_len = len(str(peptide.seq))
                t_len = len(str(protein.seq))

                pattern = str(peptide.seq)
                text = str(protein.seq)

                if p_len > t_len:
                    continue

                shift = []
                for i in range(91):
                    shift.append(p_len)

                for j in range(p_len - 1):
                    shift[ord(pattern[j])] = p_len - j - 1

                shift = tuple(shift)
                k = p_len - 1
                
                while k < t_len:
                    j = p_len - 1
                    i = k
                    while j >= 0 and text[i] == pattern[j]:
                        j -= 1
                        i -= 1
                    if j == -1:
                        all_matches.append((str(peptide.seq), str(peptide.seq), str(protein.id), 0, i + 1))

                    k += shift[ord(text[k])]

        return all_matches


class Benchmarker(Horspool):
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
        return 'Horspool algorithm'

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
