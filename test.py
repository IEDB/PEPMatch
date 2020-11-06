#!/usr/bin/env python3

from Bio import SeqIO


def parse_fasta(file):
    return SeqIO.parse(file, 'fasta')

protein_dict = {}

for item in [protein for protein in list(parse_fasta('./proteomes/9606_uniprot_small.fa'))]:
    sequence = str(item.seq).replace('X', '')
    sequence = sequence.replace('U', '')

    if len(sequence) < 18:
    	continue

    protein_dict[str(item.description)] = sequence

with open('9606_uniprot_small_test.fa', 'w') as file:
    for key, value in protein_dict.items():
        file.write('>' + key + '\n')
        file.write(value + '\n')
