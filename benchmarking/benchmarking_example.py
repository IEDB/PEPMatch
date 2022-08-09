#!/usr/bin/env python3


class Benchmarker(object):
    def __init__(self, lengths, mismatches):
        '''
        The only required arguments here is the lengths of
        the query and the number of mismatches.
        However, you can initialize whatever is necessary to
        interface with your tool for optimization.
        If your tool does not do mismatching, raise a TypeError.
        '''
        self.lengths = lengths
        self.mismatches = mismatches

        if self.max_mismatches != 0:
            raise ValueError(self.__str__() + ' cannot do any mismatching.\n')

    def __str__(self):
        '''
        Overwrite __str__ method to give the name of your tool.
        '''
        return 'Example tool'

    def preprocess_proteome(self, proteome):
        '''
        Call method from your tool that preprocesses the 
        proteome for benchmarking. Do not return anything.
        This is not applicable for every tool.
        '''
        raise TypeError(self.__str__() + ' does not preprocess queries.\n')

    def preprocess_query(self, query):
        '''
        Call method from your tool that preprocesses the 
        query for benchmarking. Do not return anything.
        This is not applicable for every tool.
        '''
        raise TypeError(self.__str__() + ' does not preprocess queries.\n')

    def search(self, query, proteome):
        '''
        Call method from your tool that does the actual
        searching from the query to the proteome. Return the
        list of results with the query peptide, matched peptide, 
        protein ID, # of mismatches, and index position.
        '''

        # Example:
        results = ['YLLDHLSYL,YLLDHLSYL,sp|O60337|MARH6_HUMAN,0,561', ]

        return results