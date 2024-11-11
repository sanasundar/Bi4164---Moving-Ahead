#insert your fetch pickle file code
clf= #something idk

#NEW CODE BTW PLS UNCOMMENT
'''
from itertools import product
from Bio import SeqIO
import os
import numpy as np
import pandas as pd
import pickle
from collections import Counter
from textwrap import wrap
from sklearn.feature_extraction import DictVectorizer

dir_data=r"d:\\github\\Bi4164---Moving-Ahead\\data_set"


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

    output= []
    for i in range(kmer_size):
        output+=wrap(sequence[i:], width=kmer_size)
    output= [i for i in output if len(i)==kmer_size]
    dict_output= Counter(output)
    return output, dict_output

def create_matrix_labels_given_order( order):
    
    nucleotides= ["A", "T", "G", "C"]
    counter=0
    if order==0:
        return ["NULL"]
    if order==1:
        return nucleotides
    else:
        past_output= create_matrix_labels_given_order(order-1)
        new_output= []
        for entries in past_output:
            for nucleotide in nucleotides:
                new_output.append(entries + nucleotide)
        return new_output

def get_vectors_for_all_sequences(sequences_list, kmer_size):
    
    v = DictVectorizer()
    keys_null= create_matrix_labels_given_order(kmer_size)
    values_null= np.zeros(shape= (len(keys_null)))
    dict_null= dict(zip(keys_null, values_null))
    D = [dict_null]
    for i in range(0, len(sequences_list)):
        kmers_for_sequence = dict(make_kmers(sequences[i], kmer_size)[1])
        D.append(kmers_for_sequence)
    print('kmers have been made and vectors acquired!')
    X = v.fit_transform(D)
    X= X[1:]
    return X, v
        
kmer_size = 6
sequences, ids = read_fasta(fasta_sequences_neg, fasta_sequences_pos)
# normalised_vectors_for_all_data = get_vectors_for_all_sequences(sequences, kmer_size)
# print(normalised_vectors_for_all_data[0])
# print(sequences[0])
# print(dict(make_kmers(sequences[0], kmer_size)[1]))
X, v= get_vectors_for_all_sequences(sequences, kmer_size)
onehot_ids= [1 if i=="pos" else 0 for i in ids]
y=onehot_ids

list_of_sequences = []
dir_data=r"d:\\github\\Bi4164---Moving-Ahead\\data_set"

os.chdir(dir_data)
fasta_sequences_test = SeqIO.parse(open("test.fa"), 'fasta')

for fasta in fasta_sequences_test:
    name, sequence = fasta.id, str(fasta.seq)
    list_of_sequences.append(sequence)
    
X_testing, v_testing= get_vectors_for_all_sequences(list_of_sequences, kmer_size)
'''


model_name= "linearSVC" #change model name
dir_outputs=r"d:\\github\\Bi4164---Moving-Ahead\\model_outputs"
prob_testingFA= clf.predict(X_testing)
# neg_prob= prob_testingFA[:,0]
# pos_prob= prob_testingFA[:,1]
os.chdir(dir_outputs)
np.save(model_name+"_negative_0_positive_1.npy", neg_prob)