from itertools import product
from Bio import SeqIO
import os
import numpy as np
import pandas as pd
import pickle
from collections import Counter
from textwrap import wrap
from sklearn.feature_extraction import DictVectorizer
from sklearn.model_selection import KFold
from sklearn.model_selection import StratifiedKFold
from sklearn import datasets, metrics, model_selection, svm
import sklearn
import matplotlib.pyplot as plt
from sklearn.metrics import RocCurveDisplay
from sklearn.naive_bayes import BernoulliNB
from sklearn.metrics import precision_score

dir_data=r"d:\\github\\Bi4164---Moving-Ahead\\data_set"
plot= r"d:\\github\\Bi4164---Moving-Ahead\\plot"

os.chdir(dir_data)
fasta_sequences_neg = SeqIO.parse(open("negative.fa"), 'fasta')
fasta_sequences_pos = SeqIO.parse(open("positive.fa"), 'fasta')

def make_dir(new_dir):
    if not os.path.exists(new_dir):
        os.makedirs(new_dir)
    return "made dir " + str(new_dir)

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

def get_vectors_for_all_sequences(sequences_list, kmer_size):
    D = []
    v = DictVectorizer()
    for i in range(0, len(sequences_list)):
        kmers_for_sequence = dict(make_kmers(sequences[i], kmer_size)[1])
        D.append(kmers_for_sequence)
    print('kmers have been made and normalised vectors have been acquired!')
    X = v.fit_transform(D)
    return X, v

        
kmer_sizes = [12,14,16,18,20,22,24,26]
sequences, ids = read_fasta(fasta_sequences_neg, fasta_sequences_pos)

for kmer_size in kmer_sizes:

    
    X, v= get_vectors_for_all_sequences(sequences, kmer_size)
    onehot_ids= [1 if i=="pos" else 0 for i in ids]
    y=onehot_ids
    kf = StratifiedKFold(n_splits=3)
    kf.get_n_splits(X, y)
    # clf, clf_label= KNeighborsClassifier(n_jobs=4, n_neighbors=3), "knn"
    
    # clf, clf_label= LogisticRegression(), "linear"
    # clf, clf_label= RidgeClassifier(alpha= 1.0), "ridge"
    # clf, clf_label= RidgeClassifier(alpha= 0.5), "ridge_0.5"
    # X_test, Y_test, X_train, Y_train, Pred, Pred_train= [], [], [], [], [], []
    kmer_plot= plot+f"\\kmer_size_{kmer_size}"

    for i, (train_index, test_index) in enumerate(kf.split(X, y)):
        clf = BernoulliNB()
        x_train, x_test = X[train_index.astype(int)], X[test_index.astype(int)]
        y_train, y_test = np.array(y)[train_index.astype(int)], np.array(y)[test_index.astype(int)]
        clf.fit(x_train, y_train)
        
        y_pred = clf.predict(x_test)
        y_train_pred = clf.predict(x_train) 
        print(f"fold number: {i}, precision test: {precision_score(y_test, y_pred)}, precision train: {precision_score(y_train, y_train_pred)}")
        y_pred = clf.predict_proba(x_test)[:,1]
        fpr, tpr, thresholds = metrics.roc_curve(y_test, y_pred)
        print(fpr, tpr, thresholds)
        # roc_auc = metrics.auc(fpr, tpr)
        
        RocCurveDisplay.from_predictions(y_test, y_pred)
        # display.plot()
        make_dir(kmer_plot)
        os.chdir(kmer_plot)
        plt.savefig(f"ROC_fold_number_{i}.png")
