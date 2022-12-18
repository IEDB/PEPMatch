#!/usr/bin/env python3

from Bio import SeqIO

def parse_fasta(file):
    return SeqIO.parse(file, 'fasta')

class HammingDistance(object):
    def __init__(self, query, proteome, max_mismatches):
        self.query = query
        self.proteome = proteome
        self.max_mismatches = max_mismatches

    def hamming_distance(self, string1, string2): 
        len1 = len(string1)
        len2 = len(string2)
        assert len1 == len2, 'Strings not the same length'

        mismatches = 0
        for i in range(len(string1)):
            if string1[i] != string2[i]:
                mismatches +=1
                if mismatches > self.max_mismatches:
                    return None

        return mismatches

    def mismatching_search(self):
        query = list(parse_fasta(self.query))
        proteome = list(parse_fasta(self.proteome))

        peptide_sequences = []
        matched_peptides = []
        protein_hits = []
        index_start = []
        mismatches = []

        for peptide in query:
            k = len(peptide.seq)
            for protein in proteome:
                for i in range(len(protein.seq) - k + 1):
                    distance = self.hamming_distance(str(peptide.seq), str(protein.seq[i:i+k]))
                    if distance != None:
                        peptide_sequences.append(str(peptide.seq))
                        matched_peptides.append(str(protein.seq[i:i+k]))
                        protein_hits.append(str(protein.description).split(' ')[0])
                        index_start.append(i)
                        mismatches.append(distance)

        all_matches = []
        for i in range(len(peptide_sequences)):
            all_matches.append((
                peptide_sequences[i], 
                matched_peptides[i], 
                protein_hits[i], 
                index_start[i],
                mismatches[i])
            )

        return all_matches

class Benchmarker(HammingDistance):
    def __init__(self, query, proteome, lengths, max_mismatches, method_parameters):
        self.query = query
        self.proteome = proteome
        self.lengths = lengths
        self.max_mismatches = max_mismatches
        self.method_parameters = method_parameters
        super().__init__(query, proteome, max_mismatches)

    def __str__(self):
        return 'Hamming distance search algorithm'

    def preprocess_proteome(self):
        raise TypeError(self.__str__() + ' does not preprocess proteomes.\n')

    def preprocess_query(self):
        raise TypeError(self.__str__() + ' does not preprocess queries.\n')

    def search(self, query, proteome):
        matches = self.mismatching_search()

        all_matches = []
        for match in matches:
            if match[1] == '':
                continue
            
            match_string = ''
            for i in range(0, len(match)):

                if i == len(match) - 1:
                    match_string += str(match[i])
                elif i == 2:
                    match_string += match[i].split('|')[1]
                else:
                    match_string += str(match[i]) + ','

            all_matches.append(match_string)
        
        return all_matches