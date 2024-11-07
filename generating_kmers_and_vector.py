from itertools import product
from Bio import SeqIO
import os
import numpy as np
import pandas as pd
import pickle

dir_data = r"C:\Users\Saranya Sundar\OneDrive\Desktop\Python\IISER Pune Assignments\Bi4164\\"

os.chdir(dir_data)
fasta_sequences_neg = SeqIO.parse(open("negative.fa"), 'fasta')
fasta_sequences_pos = SeqIO.parse(open("positive.fa"), 'fasta')

def read_fasta(fasta_sequences_neg, fasta_sequences_pos):
    '''Returns a list of sequences and their corresponding genomic activities i.e., promoter or not
    
        Parameters
        ----------
        fasta_sequences_neg: class 'Bio.SeqIO.FastaIO.FastaIterator'
            The list of sequences that do not show promoter activity
        fasta_sequences_pos: class 'Bio.SeqIO.FastaIO.FastaIterator'
            The list of sequences that show promoter activity
        
        Returns
        -------
        list_of_sequences: list
            The list of combined negative and positive sequences
        list_of_ids: list
            The list of ids, which is whether the corresponding sequence has promoter activity or not'''
    
    list_of_sequences = []
    list_of_ids = []
    for fasta in fasta_sequences_neg:
        name, sequence = fasta.id, str(fasta.seq)
        list_of_sequences.append(sequence)
        list_of_ids.append('neg')
    for fasta in fasta_sequences_pos:
        name, sequence = fasta.id, str(fasta.seq)
        list_of_sequences.append(sequence)
        list_of_ids.append('pos')
    print('read_fasta has run!')
    return list_of_sequences, list_of_ids

def make_kmers(sequence, kmer_size):
    '''Returns the list of kmers for a particular sequence and size of kmer inputted by the user
    
        Parameters
        ----------
        sequence: string
            The sequence for which kmers are to be generated 
        kmer_size: int
            The size of the sliding window to generate kmers of that size
        
        Returns
        -------
        kmers: list
            The list of kmers for the sequence'''

    kmers = []
    number_of_kmers = len(sequence) + 1 - kmer_size
    for i in range (0, number_of_kmers):
        kmer = sequence[i: 1+kmer_size]
        kmers.append(kmer)
    return kmers

def possible_kmers(kmer_size):
    '''Returns all possible combinations of the four letter nucleotide for a particular window size
    
        Parameters
        ----------
        kmer_size: int
            The sliding window size or size of each kmer considered
        
        Returns
        --------
        combinations: list
            List of all possible kmers'''

    nucleotides = ['A', 'T', 'G', 'C']
    combinations = [''.join(combinations) for combinations in product(nucleotides, repeat=kmer_size)]
    print('possible_kmers has run!')
    return combinations

def define_vector(kmer_from_sequence, possible_kmers):
    '''Returns the normalised vector for the sequence
        
        Parameters
        ----------
        kmer_from_sequence: list
            The list of kmers for the particular sequence
        possible_kmers: list
            List of all possible kmers for sliding window size kmer_size
        
        Returns
        -------
        normalised_vector: numpy array
            An array of counts of each kmer occuring in the sequence and normalised to the number of kmers in the sequence to be as little biased to size of sequence as possible'''
    
    vector = np.zeros(len(possible_kmers))
    for i in range(0, len(kmer_from_sequence)):
        if kmer_from_sequence[i] in possible_kmers:
            index = possible_kmers.index(kmer_from_sequence[i])
            vector[index] += 1
    normalised_vector = vector/len(kmer_from_sequence)
    return normalised_vector

def get_vectors_for_all_sequences(sequences_list, kmer_size):
    normalised_vectors = []
    kmers_possible_for_kmer_size = possible_kmers(kmer_size)
    for i in range(0, len(sequences_list)):
        kmers_for_sequence = make_kmers(sequences_list[i], kmer_size)
        normalised_vectors.append(define_vector(kmers_for_sequence, kmers_possible_for_kmer_size))
    print('kmers have been made and normalised vectors have been acquired!')
    return normalised_vectors
        
kmer_size = 6
sequences, ids = read_fasta(fasta_sequences_neg, fasta_sequences_pos)
normalised_vectors_for_all_data = get_vectors_for_all_sequences(sequences, kmer_size)
print(normalised_vectors_for_all_data)