#!/usr/bin/env python3

from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastpCommandline
import os, glob
import pandas as pd
import json

directory = os.path.dirname(os.path.abspath(__file__))

def parse_fasta(file):
    return SeqIO.parse(file, 'fasta')


class BLAST(object):
    def __init__(self, query, proteome, max_mismatches, algorithm_parameters):
        if max_mismatches == -1:
            raise ValueError(self.__str__() + ' does not have a best match feature.\n')

        self.query = query
        self.proteome = proteome

        self.max_mismatches = max_mismatches

        bin_directory = algorithm_parameters['bin_directory']
        self.makeblastdb_path = os.path.join(bin_directory, 'makeblastdb')
        self.blastp_path = os.path.join(bin_directory, 'blastp')

    def __str__(self):
        return 'BLAST'
    
    def preprocesss(self):
        os.system(self.makeblastdb_path + ' -in ' + self.proteome + ' -dbtype prot')
    
    def blast_search(self, query, proteome):
        peptides = parse_fasta(self.query)
        proteins = parse_fasta(self.proteome)

        peptide_dict = {}

        i = 0
        for peptide in peptides:
            peptide_dict[str(peptide.id)] = str(peptide.seq)
            i += 1

        protein_dict = {}
        i = 0
        for protein in proteins:
            protein_dict[str(protein.id)] = str(protein.seq)
            i += 1

        if self.max_mismatches == 0:
            blastx_cline = NcbiblastpCommandline(cmd=self.blastp_path,
                                                query = query, 
                                                db = proteome, 
                                                evalue=100, outfmt=10, out='output.csv')

            stdout, stderr = blastx_cline()

        else:
            blastx_cline = NcbiblastpCommandline(cmd=self.blastp_path, 
                                                query = query, 
                                                db = proteome, 
                                                evalue=10000, outfmt=10, out='output.csv')

            stdout, stderr = blastx_cline()

        df = pd.read_csv('output.csv', names = ['Peptide Sequence', 'Protein ID', 'Sequence Identity', 
                                                'Length', 'Mismatches', 'Gap Openings', 'Query start', 
                                                'Query end', 'Index start', 'Index end', 'e value', 'bit score'])

        
        df['Peptide Sequence'] = df['Peptide Sequence'].apply(str)
        df = df.replace({'Peptide Sequence': peptide_dict})

        df['Index start'] = df['Index start'].apply(lambda x: x - 1)


        all_matches = []
        for i, row in df.iterrows():
            all_matches.append((
                row['Peptide Sequence'], 
                protein_dict[row['Protein ID']][row['Index start']:int(row['Index start'])+len(row['Peptide Sequence'])], 
                row['Protein ID'],
                row['Mismatches'],
                row['Index start'],
                ))

        return all_matches


class Benchmarker(BLAST):
    def __init__(self, query, proteome, lengths, max_mismatches, algorithm_parameters):
        self.query = query
        self.proteome = proteome
        self.lengths = lengths
        self.max_mismatches = max_mismatches
        self.algorithm_parameters = algorithm_parameters
        
        super().__init__(query, proteome, max_mismatches, algorithm_parameters)

    def __str__(self):
        return 'BLAST'

    def preprocess_proteome(self):
        return self.preprocesss()

    def preprocess_query(self):
        raise TypeError(self.__str__() + ' does not preprocess queries.\n')

    def search(self):
        matches = self.blast_search(self.query, self.proteome)

        all_matches = []
        for match in matches:
            match_string = ''
            for i in match:
                if i == match[-1]:
                    match_string += str(i)
                else:
                    match_string += str(i) + ','

            all_matches.append(match_string)

        for extension in ['pdb', 'phr', 'pin', 'psq', 'ptf', 'pot', 'pto']:
            os.remove(glob.glob(os.path.dirname(self.proteome) + '/*.' + extension)[0])

        os.remove('output.csv')

        return all_matches
