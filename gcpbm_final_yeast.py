# -*- coding: utf-8 -*-
### conda code v1

!pip install biopython
!pip install folium
!curl -O https://raw.githubusercontent.com/jperkel/example_notebook/master/NC_005816.gb

##################

import sys
import os
from itertools import product
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt1
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestRegressor, RandomForestClassifier
from sklearn.preprocessing import MinMaxScaler
from sklearn import metrics
from sklearn.linear_model import LinearRegression
from sklearn.linear_model import Ridge
from sklearn.svm import SVR
from Bio import pairwise2

from collections import defaultdict
import random

protein = 'cbf1'
concentration = '100' # '100' or '200'
experiment = 'pb' # 'chip' or 'pb' for in vivo or in vitro resp.

# yeast

raw_data = pd.read_csv(f'drive/MyDrive/ML/gcPBM/{protein}/{protein}_{concentration}.txt', sep='\t')
proc_data = pd.concat([raw_data[["Sequence", "Alexa488"]]])
proc_data = proc_data.dropna()
proc_data = proc_data.reset_index()
strings = list(proc_data["Sequence"])
strings = [seq[3:33] for seq in strings]
values = proc_data["Alexa488"]

min_val = min(values)
normalization= max(values)-min(values)
values = list((values-min_val)/normalization)

dictionary = dict(zip(values, strings))

plt.plot(sorted(values))
strings[0][12:18]



# # hooman

# # Read and process raw data from input files 

# raw_data = pd.read_csv(f'drive/MyDrive/ML/gcPBM/sequences.txt', sep='\t')
# raw_values = pd.read_csv(f'drive/MyDrive/ML/gcPBM/{protein}/{protein}_values.txt', sep='\t')

# proc_data = pd.concat([raw_data[["ID", "Name", "SEQUENCE"]],raw_values [["VALUE"]]],axis=1)

# proc_data = proc_data.dropna()
# proc_data = proc_data.reset_index()

# indices = []
# for i in range(len(proc_data)):
#   if proc_data["Name"][i][0] == 'U' or proc_data["Name"][i][0] == 'T':
#     indices.append(i)
# proc_data = proc_data.drop(indices)

# # Processed data

# strings = list(proc_data["SEQUENCE"])
# #strings = [seq[13:36] for seq in strings]
# values = proc_data["VALUE"]

# min_val = min(values)
# normalization= max(values)-min(values)
# values = list((values-min_val)/normalization)

# dictionary = dict(zip(values, strings))

# plt.plot(sorted(values))



freq_matrix = pd.read_csv(f'drive/MyDrive/ML/gcPBM/cbf1/cbf1_freq_matrix_6.txt', delim_whitespace=True)
freq_matrix = np.array(freq_matrix)
translate = {'A':0, 'C':1, 'G':2, 'T':3}
reverse = {0:'A', 1:'C', 2:'G', 3:'T'}

len_aln = len(strings[0])

# scores

l =  list(np.arange(0,12-5))
for k in np.arange(18,30-5):
  l.append(k)

discarded = []
for s in np.arange(len(strings)):
    string = strings[s]
    scores = []
    for i in l:
        score = 0
        for j in np.arange(0,len(freq_matrix)):
            a = translate[string[i+j]]
            score = score + freq_matrix[j][a]
        scores.append(score)
    if any(sc > 2.5 for sc in scores):
        discarded.append(s)

# discarded = set(discarded)
# for s in np.arange(len(strings)):
#     string = strings[s]
#     i = 12
#     score = 0
#     for j in np.arange(0,len(freq_matrix)):
#           a = translate[string[i+j]]
#           score = score + freq_matrix[j][a]
#     if score < 1.1:
#         discarded.add(s)

discarded = list(discarded)
print('We have discarded',len(discarded)/len(strings)*100,'% of the data')



# ## get freq mat

# max_aff = max(values)
# sort_aff = sorted(list(values))
# upper_bound = sort_aff[-50]

# with open(f'drive/MyDrive/ML/gcPBM/{protein}/{protein}_FASTA_seq.txt','w') as file:
#     for i in np.arange(0,len(strings)):
#         if values[i]>upper_bound:
#             file.write(f'>{i}\n{strings[i]}\n')

# ## read freq mat

# freq_matrix = pd.read_csv(f'drive/MyDrive/ML/gcPBM/{protein}/{protein}_freq_matrix_36.txt', delim_whitespace=True)
# freq_matrix = np.array(freq_matrix)

# translate = {'A':0, 'C':1, 'G':2, 'T':3}
# reverse = {0:'A', 1:'C', 2:'G', 3:'T'}

# len_aln = len(strings[0])

# ## scores

# scores = []
# for seq in strings:
#     score = 0
#     for j in np.arange(0,len_aln):
#         a = translate[seq[j]]
#         score = score + freq_matrix[j][a]
#     scores.append(score)

# # boundaries
# upp_score = sorted(scores)[::-1][50]

# # file
# to_write = []
# weights = []

# for seq, val, score in zip(strings, values, scores):
#     if score<upp_score:
#         weights.append(1)
#         to_write.append((seq,val, 1))
#     if score>=upp_score:
#         weights.append(10)
#         to_write.append((seq,val, 10))

# #### undersampling       ## comment if not using it


# # number of divisions we are using
# L = 25000

# # affinitiy boundaries
# sort_aff = sorted(values)
# # pick top points ~ proportion 1:500 before undersampling i.e. we pick 50 top points not to undersample
# upper_bound = sort_aff[-1000]

# # interval dictionary -- 
# interval_dict = defaultdict(list)
# for value in values:
#     interval_dict[int(value*L)].append(value)

# # creates a matrix of values: we pick one value for the intervals where i<L/2
# # and up to two for i>L/2

# training_ordered = []
# start = 1 # i.e. starts picking values larger than start/L from the normalized affinites

# for i in np.arange(start,L+1):
#     if len(interval_dict[i])>0:
#         if i <= L*upper_bound:
#             random_values = [random.choice(interval_dict[i]) for j in np.arange(0,20)]
#             random_values = set(random_values)
#             value = random.choice(interval_dict[i])
#             training_ordered.append([dictionary[value], value])
#         if i > L*upper_bound:
#             for value in interval_dict[i]:
#                 training_ordered.append([dictionary[value], value])

# # writes the file w/o randomization (not used but in case we need it)

# with open(f'drive/MyDrive/ML/gcPBM/{protein}/{protein}_training_ordered.txt','w') as file:
#     file.write('ID_REF\tVALUE\n')

# with open(f'drive/MyDrive/ML/gcPBM/{protein}/{protein}_training_ordered.txt','a') as file:
#     for vector in training_ordered:
#         file.write("%s\t" % vector[0])
#         file.write("%s\n" % vector[1])


# # writes the randomized file ready for training


# np.random.shuffle(training_ordered)

# with open(f'drive/MyDrive/ML/gcPBM/{protein}/{protein}_training.txt','w') as file:
#     file.write('ID_REF\tVALUE\n')
    
# with open(f'drive/MyDrive/ML/gcPBM/{protein}/{protein}_training.txt','a') as file:
#     for vector in training_ordered:
#         file.write("%s\t" % vector[0])
#         file.write("%s\n" % vector[1])

# #### CLASSES: NO WEIGHTING
l = []
features = []
class Dataset:
    def __init__(self, protein, fitxer, model_target='tetramers', randomize_fce=False,
                 chosen_features=['full_fce', 'avg'], score='Median_intensity',
                 selected_tetramers = [0, 1, 2, 3, 4, 5, 6, 7, 8], pbm = 10):

        self.protein = protein
        self.pbm35 = fitxer

        ##self.pbm35 = pd.read_csv(f'drive/My Drive/ML/{protein}/{fitxer}.txt', sep='\t')
        ##self.pbm35 = pd.read_csv(f'fitxer, sep='\t')

        self.pbm35 = self.pbm35[["ID_REF","VALUE"]]
        self.pbm35.columns = ["ID_REF", "VALUE"]
        self.pbm35 = self.pbm35.dropna()


        inv_seq_data = self.pbm35.copy()
        inv_seq_data["ID_REF"] = inv_seq_data["ID_REF"].apply(self.inv_seq)
        self.pbm35 = pd.concat([self.pbm35, inv_seq_data])
        self.pbm35.reset_index(inplace = True)

        self.feats = chosen_features
        self.score = score
        self.selected_tetramers = selected_tetramers
        self.randomize = randomize_fce
        self.model_target = model_target           # TODO include inverse sequences as well

        #Get list of all kmers for indices purposes
        self.sequences = list(self.pbm35["ID_REF"])
        self.p = len(self.sequences[0])

        #Get list of all tetramers per each octamer in order to get features
        self.featurize()

    def __str__(self):
        return f"{self.model_target}-based dataset for protein {self.protein} using features {self.feats}"

    def featurize(self):
        if 'avg' in self.feats:
            tetra_avg = {line.split()[0] : np.array([float(x) for x in line.split()[1:]])  for line in open(f'drive/MyDrive/ML/input/avg_tetramer.dat') if 'SHIFT' not in line}
            self.tetra_avg = np.array([np.concatenate([tetra_avg[otmer[i:i+4]] for i in self.selected_tetramers]) for otmer in self.sequences])
        if 'full_fce' in self.feats or 'diagonal_fce' in self.feats:
            tetra_fce = {line.split()[0] : np.array([float(x) for x in line.split()[1:]])  for line in open(f'drive/MyDrive/ML/input/fce_tetramer.dat') if 'SHIFT' not in line}
            if self.randomize:
                keys = list(tetra_fce.keys())
                permut_keys = np.random.permutation(keys)
                tetra_fce = {key: tetra_fce[val] for key, val in zip(keys, permut_keys)}
            self.tetra_fce = np.array([np.concatenate([tetra_fce[otmer[i:i+4]] for i in self.selected_tetramers]) for otmer in self.sequences])
        # We might want to scramble the matchings to do 'negative control', i.e., remove all physical information from the dataset
        # Let's also keep track of the reduced matrix
        if 'diagonal_fce' in self.feats:
            tetra_fce_reduced = {tt: tetra_fce[tt][list(range(0,36,7))] for tt in tetra_fce.keys()}
            self.tetra_fce_reduced = np.array([np.concatenate([tetra_fce_reduced[otmer[i:i+4]] for i in self.selected_tetramers]) for otmer in self.sequences])
        if 'onehot_1mer' in self.feats:
#             self.onehot_1mer = np.array([self.onehot_encoding(otmer) for otmer in self.sequences]).astype(np.int8)
            self.onehot_1mer = np.array([self.onehot_encoding(otmer) for otmer in self.sequences])
        if 'integer' in self.feats:
            # define encoding input values
            nucleotides = product('ACGT', repeat=4)
            char_to_int = dict((''.join(c), i) for i, c in enumerate(nucleotides))
            self.integer = np.array([[char_to_int[otmer[i:i+4]] for i in self.selected_tetramers] for otmer in self.sequences]).astype(int)
        if 'presence_tetramer' in self.feats:
            self.presence_tetra = np.array([self.presence(otmer, 4) for otmer in self.sequences]).astype(np.int8)
        #####
        if 'mgw' in self.feats:
            self.mgw = {line.split()[0] : [float(x) for x in line.split()[1:]]  for line in open(f'drive/MyDrive/ML/input/mgw_rohs.txt') if 'SHIFT' not in line}
            self.mgw = np.array([np.concatenate([self.mgw[otmer[i:i+4]] for i in self.selected_tetramers]) for otmer in self.sequences])

        if 'electrostatic' in self.feats:
            self.electrostatic = {line.split()[0] : [x for x in line.split()[1:]]  for line in open(f'drive/MyDrive/ML/input/electrostatic.txt') if 'SHIFT' not in line}
            self.electrostatic = np.array([np.concatenate([self.electrostatic[otmer[i:i+4]] for i in self.selected_tetramers]) for otmer in self.sequences])


    @staticmethod
    def onehot_encoding(sequence):
        # define encoding input values
        nucleotides = 'ACGT'
        # define mapping of chars to integers and viceversa
        char_to_int = dict((c, i) for i, c in enumerate(nucleotides))
        # integer encode input data
        integer_encoded = [char_to_int[char] for char in sequence]
        # one hot encode
        onehot_encoded = list()
        for value in integer_encoded:
            letter = [0 for _ in range(len(nucleotides))]
            letter[value] = 1
            onehot_encoded.extend(letter)
        return(np.array(onehot_encoded))

    @staticmethod
    def presence(sequence, k):
        kmers = np.zeros(4**k)
        positions = {''.join(x): n for n, x in enumerate(list(product('ACTG',repeat=k)))}
        for i in range(len(sequence)-k+1):
#            kmers[positions[sequence[i:i+k]]] = 1
            kmers[positions[sequence[i:i+k]]] += 1
        return kmers

    @property
    def scores(self):
        keyw = "VALUE" if self.score == 'Median_intensity' else "Z-score"
        if self.model_target == "octamers":
            vals = self.pbm35[keyw].values
            return MinMaxScaler().fit_transform(vals.reshape(-1,1))
        elif self.model_target == "tetramers":
            mean_scores = self.mean_score()
            vals = np.array([[mean_scores[i][otmer[i:i+4]] for i in self.selected_tetramers] for otmer in self.sequences])
            return MinMaxScaler().fit_transform(vals)

    @staticmethod
    def inv_seq(seq):
        complementary = {"A": "T", "T": "A", "C": "G", "G": "C"}
        return ''.join([complementary[x] for x in seq[::-1]])


    @property
    def features(self):
        avail = {'avg': 'tetra_avg', 'full_fce': 'tetra_fce', 'diagonal_fce': 'tetra_fce_reduced', 'onehot_1mer': 'onehot_1mer', \
                 'integer': 'integer', 'presence_tetramer': 'presence_tetra', 'mgw': 'mgw', 'electrostatic': 'electrostatic'}
        return np.hstack([self.__getattribute__(avail[feat]) for feat in self.feats])

    def mean_score(self):
        """
        generates a list of dictionaries, each containing the position-wise
        score per tetramer
        """
        keyw = "VALUE" if self.score == 'Median_intensity' else "Z-score"
        if self.model_target == "tetramers":
            expt_scores = self.pbm35[keyw].values
            from collections import defaultdict
            position_scores = [defaultdict(lambda: 0) for _ in range(self.p - 3)]
            position_counts = [defaultdict(lambda: 0) for _ in range(self.p - 3)]
            for seq, score in zip(self.sequences, expt_scores):
                for i in range(self.p - 3):
                    position_scores[i][seq[i:i+4]] += score
                    position_counts[i][seq[i:i+4]] += 1
        return [{ttmer: position_scores[i][ttmer]/position_counts[i][ttmer] for ttmer in position_scores[i].keys()} for i in range(self.p - 3)]

class Model:
    regressors = {'random_forest': RandomForestRegressor(n_estimators = 10, random_state =None), \
                  'linear': LinearRegression(), 'ridge': Ridge(alpha=0.0001), 'svr': SVR()}
    def __init__(self, dataset, regressor='random_forest', training_set_size=0.8):
        plt.clf()
        self.data = dataset
        self.X = self.data.features
        self.y = self.data.scores

        # self.weights = self.data.pbm35['WEIGHT']

        print("features have shape: ", self.X.shape)
        print("target values have shape: ", self.y.shape)
        self.regressor = Model.regressors[regressor]
        # splits w/o randomness for training and testing
        self.X_train, self.X_test, self.y_train, self.y_test = train_test_split(self.X, self.y, train_size=training_set_size, random_state=0, shuffle=False)
        self.regressor.fit(self.X_train, self.y_train)
        # self.X_train, self.X_test, self.y_train, self.y_test, self.w_train, self.w_test = train_test_split(self.X, self.y, self.weights, train_size=training_set_size, random_state=0, shuffle=False)
        # self.regressor.fit(self.X_train, self.y_train, self.w_train)  # with weight
        importances = self.regressor.feature_importances_
        std = np.std([tree.feature_importances_ for tree in self.regressor.estimators_],axis=0)
        indices = np.argsort(importances)[::-1]

        # Print the feature ranking
        print("Feature ranking:")

        for f in range(self.X.shape[1]):
            if indices[f]<256:
              feat = 'Presence'
            if indices[f] in np.arange(256,256+(len_aln-3)*6,1):
              if indices[f]%6 == 0:
                feat = 'AVG Shift'
              if indices[f]%6 == 1:
                feat = 'AVG Slide'
              if indices[f]%6 == 2:
                feat = 'AVG Rise'
              if indices[f]%6 == 3:
                feat = 'AVG Tilt'
              if indices[f]%6 == 4:
                feat = 'AVG Roll'
              if indices[f]%6 == 5:
                feat = 'AVG Twist'
            if indices[f] in np.arange(256+(len_aln-3)*6,256+(len_aln-3)*12,1):
              if indices[f]%6 == 0:
                feat = 'Diagonal Shift'
              if indices[f]%6 == 1:
                feat = 'Diagonal Slide'
              if indices[f]%6 == 2:
                feat = 'Diagonal Rise'
              if indices[f]%6 == 3:
                feat = 'Diagonal Tilt'
              if indices[f]%6 == 4:
                feat = 'Diagonal Roll'
              if indices[f]%6 == 5:
                feat = 'Diagonal Twist'
            if indices[f]>=256+(len_aln-3)*12:
              feat = 'Electro'
 #       for f in range(self.X.shape[1]):
            # print(feat)
            # print("%d. feature %d (%f)" % (f + 1, indices[f], importances[indices[f]]))
            # print(feat,indices[f],importances[indices[f]])
            l.append(importances[indices[f]])
            features.append([feat,importances[indices[f]]])

        # Plot the impurity-based feature importances of the forest
        # plt.figure()
        # plt.title("Feature importances")
        # plt.bar(range(self.X.shape[1]), importances[indices], color="r", yerr=std[indices], align="center")
        # plt.xticks(range(self.X.shape[1]), indices)
        # plt.xlim([-1, self.X.shape[1]])
        # plt.savefig(f'drive/MyDrive/ML/{protein}/feature_importance.png')

    def predict(self, testing_dataset=None):
        if testing_dataset is not None:
            self.X_test, self.y_test = testing_dataset.features, testing_dataset.scores
        self.y_pred = self.regressor.predict(self.X_test).reshape(-1,1)
        self.y_test = self.y_test.reshape(-1,1)
        # self.r2 = np.round(metrics.r2_score(self.y_test, self.y_pred, self.w_test), 3) # with weights
        self.r2 = np.round(metrics.r2_score(self.y_test, self.y_pred), 3)           # w/o weights
        self.mse = np.round(metrics.mean_absolute_error(self.y_test, self.y_pred), 6)

    def plot(self):
        plt.rc('xtick', labelsize=15)
        plt.rc('ytick', labelsize=15)
        plt.xlabel('Testing')
        plt.ylabel('Predicted')
        #plt.title(f'MSE = {self.mse}, r$^2$ = {self.r2}')
        plt.xlim = (-1,1)
        plt.ylim = (-1,1)
        plt.scatter(self.y_test, self.y_pred, label=f'r$^2$ = {self.r2}')

    def show(self):
        plt.legend()
        plt.show()

# #### CLASSES: WEIGHTING

# l = []
# features = []
# class Dataset:
#     def __init__(self, protein, fitxer, model_target='tetramers', randomize_fce=False,
#                  chosen_features=['full_fce', 'avg'], score='Median_intensity',
#                  selected_tetramers = [0, 1, 2, 3, 4, 5, 6, 7, 8], pbm = 10):

#         self.protein = protein
#         self.pbm35 = fitxer

#         ##self.pbm35 = pd.read_csv(f'drive/My Drive/ML/{protein}/{fitxer}.txt', sep='\t')
#         ##self.pbm35 = pd.read_csv(f'fitxer, sep='\t')

#         self.pbm35 = self.pbm35[["ID_REF","VALUE",'WEIGHT']]
#         self.pbm35.columns = ["ID_REF", "VALUE", 'WEIGHT']
#         self.pbm35 = self.pbm35.dropna()


#         inv_seq_data = self.pbm35.copy()
#         inv_seq_data["ID_REF"] = inv_seq_data["ID_REF"].apply(self.inv_seq)
#         self.pbm35 = pd.concat([self.pbm35, inv_seq_data])
#         self.pbm35.reset_index(inplace = True)

#         self.feats = chosen_features
#         self.score = score
#         self.selected_tetramers = selected_tetramers
#         self.randomize = randomize_fce
#         self.model_target = model_target           # TODO include inverse sequences as well

#         #Get list of all kmers for indices purposes
#         self.sequences = list(self.pbm35["ID_REF"])
#         self.p = len(self.sequences[0])

#         #Get list of all tetramers per each octamer in order to get features
#         self.featurize()

#     def __str__(self):
#         return f"{self.model_target}-based dataset for protein {self.protein} using features {self.feats}"

#     def featurize(self):
#         if 'avg' in self.feats:
#             tetra_avg = {line.split()[0] : np.array([float(x) for x in line.split()[1:]])  for line in open(f'drive/MyDrive/ML/input/avg_tetramer.dat') if 'SHIFT' not in line}
#             self.tetra_avg = np.array([np.concatenate([tetra_avg[otmer[i:i+4]] for i in self.selected_tetramers]) for otmer in self.sequences])
#         if 'full_fce' in self.feats or 'diagonal_fce' in self.feats:
#             tetra_fce = {line.split()[0] : np.array([float(x) for x in line.split()[1:]])  for line in open(f'drive/MyDrive/ML/input/fce_tetramer.dat') if 'SHIFT' not in line}
#             if self.randomize:
#                 keys = list(tetra_fce.keys())
#                 permut_keys = np.random.permutation(keys)
#                 tetra_fce = {key: tetra_fce[val] for key, val in zip(keys, permut_keys)}
#             self.tetra_fce = np.array([np.concatenate([tetra_fce[otmer[i:i+4]] for i in self.selected_tetramers]) for otmer in self.sequences])
#         # We might want to scramble the matchings to do 'negative control', i.e., remove all physical information from the dataset
#         # Let's also keep track of the reduced matrix
#         if 'diagonal_fce' in self.feats:
#             tetra_fce_reduced = {tt: tetra_fce[tt][list(range(0,36,7))] for tt in tetra_fce.keys()}
#             self.tetra_fce_reduced = np.array([np.concatenate([tetra_fce_reduced[otmer[i:i+4]] for i in self.selected_tetramers]) for otmer in self.sequences])
#         if 'onehot_1mer' in self.feats:
# #             self.onehot_1mer = np.array([self.onehot_encoding(otmer) for otmer in self.sequences]).astype(np.int8)
#             self.onehot_1mer = np.array([self.onehot_encoding(otmer) for otmer in self.sequences])
#         if 'integer' in self.feats:
#             # define encoding input values
#             nucleotides = product('ACGT', repeat=4)
#             char_to_int = dict((''.join(c), i) for i, c in enumerate(nucleotides))
#             self.integer = np.array([[char_to_int[otmer[i:i+4]] for i in self.selected_tetramers] for otmer in self.sequences]).astype(int)
#         if 'presence_tetramer' in self.feats:
#             self.presence_tetra = np.array([self.presence(otmer, 4) for otmer in self.sequences]).astype(np.int8)
#         #####
#         if 'mgw' in self.feats:
#             self.mgw = {line.split()[0] : [float(x) for x in line.split()[1:]]  for line in open(f'drive/MyDrive/ML/input/mgw_rohs.txt') if 'SHIFT' not in line}
#             self.mgw = np.array([np.concatenate([self.mgw[otmer[i:i+4]] for i in self.selected_tetramers]) for otmer in self.sequences])

#         if 'electrostatic' in self.feats:
#             self.electrostatic = {line.split()[0] : [x for x in line.split()[1:]]  for line in open(f'drive/MyDrive/ML/input/electrostatic.txt') if 'SHIFT' not in line}
#             self.electrostatic = np.array([np.concatenate([self.electrostatic[otmer[i:i+4]] for i in self.selected_tetramers]) for otmer in self.sequences])


#     @staticmethod
#     def onehot_encoding(sequence):
#         # define encoding input values
#         nucleotides = 'ACGT'
#         # define mapping of chars to integers and viceversa
#         char_to_int = dict((c, i) for i, c in enumerate(nucleotides))
#         # integer encode input data
#         integer_encoded = [char_to_int[char] for char in sequence]
#         # one hot encode
#         onehot_encoded = list()
#         for value in integer_encoded:
#             letter = [0 for _ in range(len(nucleotides))]
#             letter[value] = 1
#             onehot_encoded.extend(letter)
#         return(np.array(onehot_encoded))

#     @staticmethod
#     def presence(sequence, k):
#         kmers = np.zeros(4**k)
#         positions = {''.join(x): n for n, x in enumerate(list(product('ACTG',repeat=k)))}
#         for i in range(len(sequence)-k+1):
# #            kmers[positions[sequence[i:i+k]]] = 1
#             kmers[positions[sequence[i:i+k]]] += 1
#         return kmers

#     @property
#     def scores(self):
#         keyw = "VALUE" if self.score == 'Median_intensity' else "Z-score"
#         if self.model_target == "octamers":
#             vals = self.pbm35[keyw].values
#             return MinMaxScaler().fit_transform(vals.reshape(-1,1))
#         elif self.model_target == "tetramers":
#             mean_scores = self.mean_score()
#             vals = np.array([[mean_scores[i][otmer[i:i+4]] for i in self.selected_tetramers] for otmer in self.sequences])
#             return MinMaxScaler().fit_transform(vals)

#     @staticmethod
#     def inv_seq(seq):
#         complementary = {"A": "T", "T": "A", "C": "G", "G": "C"}
#         return ''.join([complementary[x] for x in seq[::-1]])


#     @property
#     def features(self):
#         avail = {'avg': 'tetra_avg', 'full_fce': 'tetra_fce', 'diagonal_fce': 'tetra_fce_reduced', 'onehot_1mer': 'onehot_1mer', \
#                  'integer': 'integer', 'presence_tetramer': 'presence_tetra', 'mgw': 'mgw', 'electrostatic': 'electrostatic'}
#         return np.hstack([self.__getattribute__(avail[feat]) for feat in self.feats])

#     def mean_score(self):
#         """
#         generates a list of dictionaries, each containing the position-wise
#         score per tetramer
#         """
#         keyw = "VALUE" if self.score == 'Median_intensity' else "Z-score"
#         if self.model_target == "tetramers":
#             expt_scores = self.pbm35[keyw].values
#             from collections import defaultdict
#             position_scores = [defaultdict(lambda: 0) for _ in range(self.p - 3)]
#             position_counts = [defaultdict(lambda: 0) for _ in range(self.p - 3)]
#             for seq, score in zip(self.sequences, expt_scores):
#                 for i in range(self.p - 3):
#                     position_scores[i][seq[i:i+4]] += score
#                     position_counts[i][seq[i:i+4]] += 1
#         return [{ttmer: position_scores[i][ttmer]/position_counts[i][ttmer] for ttmer in position_scores[i].keys()} for i in range(self.p - 3)]

# class Model:
#     regressors = {'random_forest': RandomForestRegressor(n_estimators = 10, random_state =None), \
#                   'linear': LinearRegression(), 'ridge': Ridge(alpha=0.0001), 'svr': SVR()}
#     def __init__(self, dataset, regressor='random_forest', training_set_size=0.8):
#         plt.clf()
#         self.data = dataset
#         self.X = self.data.features
#         self.y = self.data.scores

#         self.weights = self.data.pbm35['WEIGHT']

#         print("features have shape: ", self.X.shape)
#         print("target values have shape: ", self.y.shape)
#         self.regressor = Model.regressors[regressor]
#         # splits w/o randomness for training and testing
#         self.X_train, self.X_test, self.y_train, self.y_test, self.w_train, self.w_test = train_test_split(self.X, self.y, self.weights, train_size=training_set_size, random_state=0, shuffle=False)
#         self.regressor.fit(self.X_train, self.y_train, self.w_train)  # with weight
#         importances = self.regressor.feature_importances_
#         std = np.std([tree.feature_importances_ for tree in self.regressor.estimators_],axis=0)
#         indices = np.argsort(importances)[::-1]

#         # Print the feature ranking
#         print("Feature ranking:")

#         for f in range(self.X.shape[1]):
#             if indices[f]<256:
#               feat = 'Presence'
#             if indices[f] in np.arange(256,292,1):
#               if indices[f]%6 == 0:
#                 feat = 'AVG Shift'
#               if indices[f]%6 == 1:
#                 feat = 'AVG Slide'
#               if indices[f]%6 == 2:
#                 feat = 'AVG Rise'
#               if indices[f]%6 == 3:
#                 feat = 'AVG Tilt'
#               if indices[f]%6 == 4:
#                 feat = 'AVG Roll'
#               if indices[f]%6 == 5:
#                 feat = 'AVG Twist'
#             if indices[f] in np.arange(292,298,1):
#               if indices[f]%6 == 0:
#                 feat = 'Diagonal Shift'
#               if indices[f]%6 == 1:
#                 feat = 'Diagonal Slide'
#               if indices[f]%6 == 2:
#                 feat = 'Diagonal Rise'
#               if indices[f]%6 == 3:
#                 feat = 'Diagonal Tilt'
#               if indices[f]%6 == 4:
#                 feat = 'Diagonal Roll'
#               if indices[f]%6 == 5:
#                 feat = 'Diagonal Twist'
#             if indices[f]>=298:
#               feat = 'Electro'
#  #       for f in range(self.X.shape[1]):
#             # print(feat)
#             # print("%d. feature %d (%f)" % (f + 1, indices[f], importances[indices[f]]))
#             # print(feat,indices[f],importances[indices[f]])
#             l.append(importances[indices[f]])
#             features.append([feat,importances[indices[f]]])

#         # Plot the impurity-based feature importances of the forest
#         # plt.figure()
#         # plt.title("Feature importances")
#         # plt.bar(range(self.X.shape[1]), importances[indices], color="r", yerr=std[indices], align="center")
#         # plt.xticks(range(self.X.shape[1]), indices)
#         # plt.xlim([-1, self.X.shape[1]])
#         # plt.savefig(f'drive/MyDrive/ML/{protein}/feature_importance.png')

#     def predict(self, testing_dataset=None):
#         if testing_dataset is not None:
#             self.X_test, self.y_test = testing_dataset.features, testing_dataset.scores
#         self.y_pred = self.regressor.predict(self.X_test).reshape(-1,1)
#         self.y_test = self.y_test.reshape(-1,1)
#         self.r2 = np.round(metrics.r2_score(self.y_test, self.y_pred, self.w_test), 3) # with weights
#         #self.r2 = np.round(metrics.r2_score(self.y_test, self.y_pred), 3)           # w/o weights
#         self.mse = np.round(metrics.mean_absolute_error(self.y_test, self.y_pred), 6)

#     def plot(self):
#         plt.rc('xtick', labelsize=15)
#         plt.rc('ytick', labelsize=15)
#         plt.xlabel('Testing')
#         plt.ylabel('Predicted')
#         #plt.title(f'MSE = {self.mse}, r$^2$ = {self.r2}')
#         plt.xlim = (-1,1)
#         plt.ylim = (-1,1)
#         plt.scatter(self.y_test, self.y_pred, label=f'r$^2$ = {self.r2}')

#     def show(self):
#         plt.legend()
#         plt.show()

## options:

# no smart under, no weighting
df_train = pd.DataFrame({'ID_REF': strings, 'VALUE': values})
df_train = df_train.drop(discarded)
df_train = df_train.reset_index(drop=True)
# df_train = df_train.sample(frac=1).reset_index(drop=True)[0:15000]
len(df_train)
# # smart under, no weighting
# df_train = pd.read_csv(f'drive/MyDrive/ML/gcPBM/{protein}/{protein}_training.txt', sep='\t')

# # no smart under, weighting
# df_train = pd.DataFrame({'ID_REF': strings, 'VALUE': values, 'WEIGHT':weights})
# df_train = df_train.sample(frac=1).reset_index(drop=True)[0:15000]



the_features = defaultdict(list)
the_features = {0: ['electrostatic'], 1: ['presence_tetramer'], 2: ['avg'], 3: ['diagonal_fce'], 4: ['presence_tetramer', 'avg', 'diagonal_fce','electrostatic'], 5: ['avg', 'diagonal_fce']}
len_aln = len(strings[0])

chosen_features = the_features[4]
model_target = 'octamers'
regressor = 'random_forest'
training_set_size = 0.9
randomize_fce = False
score = 'Median_intensity'
selected_tetramers = list(np.arange(0,len_aln-3))

## time

import time
start_time = time.time()
dset_train = Dataset(protein, df_train, model_target, randomize_fce, chosen_features, score, selected_tetramers)
model = Model(dset_train, regressor, training_set_size)
print('the model has been trained')
print("It took %s seconds" % (time.time() - start_time))

fig = plt.gcf()
fig.set_size_inches(5, 5)
model.predict()
print(model.y_test.shape, model.y_pred.shape)
print('The correlation is ', model.r2)
model.plot()
plt.plot([0,1],color='red')

# with open(f'output/{protein}/electro_correlations.txt','a') as file:
#    file.write("0vs%s\t" % cycle)
#    file.write("%s\n" % model.r2)

df = pd.DataFrame(features)

index = [k for k in range(model.X.shape[1]) if df[0][k] == 'Presence']
p = sum([df[1][k] for k in index])

index = [k for k in range(model.X.shape[1]) if df[0][k] == 'Electro']
e = sum([df[1][k] for k in index])

shape = 1-e-p

y = np.array([p,e,shape])

fig = plt.gcf()
fig.set_size_inches(8, 8)
labels = ["Presence", "Electrostatic", "Shape"]
patches, texts = plt.pie(y, startangle=0)
plt.legend(patches, labels, loc="best")
# plt.pie(y)
plt.show()

plt.xlabel('Feature number')
plt.ylabel('Relative importance (%)')
plt.legend('Top 30 features')
plt.bar(range(len(l[0:10])),l[0:10],color='red',align="center",)

df[0:10]

## genomic

raw_peaks = pd.read_csv(f'drive/MyDrive/ML/gcPBM/{protein}/{protein}_{experiment}-exo.gff', sep='\t', header = None)
proc_peaks = raw_peaks[[0,3,4,5,8]]
proc_peaks = proc_peaks.rename(columns={0: "chr", 3: "start", 4:"end", 5:"score", 8:"distance"})
proc_peaks["distance"] = proc_peaks["distance"].str[12:].astype(int)
proc_peaks["start"] = proc_peaks["start"]-proc_peaks["distance"]
for i in range(len(proc_peaks)):
    if proc_peaks["start"][i]<0:
        proc_peaks["start"][i]=0
proc_peaks["end"] = proc_peaks["end"]+proc_peaks["distance"]
proc_peaks = proc_peaks[["chr","start","end","score"]]
min_val = min(proc_peaks["score"])
normalization= max(proc_peaks["score"])-min(proc_peaks["score"])
proc_peaks["score"] = list(((proc_peaks["score"]-min_val)/normalization))

with open(f'drive/MyDrive/ML/gcPBM/{protein}/{protein}_{experiment}-exo.bedgraph','w') as file:
    for i in range(len(proc_peaks)):
        file.write("%s\t" % proc_peaks.iloc[i][0])
        file.write("%s\t" % proc_peaks.iloc[i][1])
        file.write("%s\t" % proc_peaks.iloc[i][2])
        file.write("%s\n" % proc_peaks.iloc[i][3])



for k in range(1,17):
  chrom = 'chr'+str(k)
  print(chrom)
  chr = proc_peaks[proc_peaks['chr']==chrom] 
  chr = chr.sort_values("start") 
  chr = chr.reset_index(drop=True)
  print(len(chr)/(chr['end'][len(chr)-1]-chr['start'][0])*100000)

# intervals: 25-26 (1), 27-28 (2), 29-30 (3), 32-33 (4), 34-35 (5)

from Bio import SeqIO

chr = 'chr14'
chr_r = 'chrXIV'
number = 2
start_bp = 100000*(number-1)+1 
end_bp= 100000*(number)
#end_bp = 230218 # chr1
#end_bp = 813184 # chr2
#end_bp = 316620 # chr3
#end_bp = 1531933 # chr4
#end_bp = 576874 # chr5
#end_bp = 270161 # chr6
#end_bp = 1090940 # chr7
#end_bp = 562643 # chr8
#end_bp = 439888 # chr9
#end_bp = 745752 # chr10
#end_bp = 666816 # chr11
#end_bp = 1078177 # chr12
#end_bp = 924431 # chr13
#end_bp = 784333 # chr14
#end_bp = 1091291 # chr15
#end_bp =  948066 # chr16

print(start_bp,end_bp,end_bp-start_bp)

fasta = SeqIO.read(f"drive/MyDrive/ML/gcPBM/{chr}_{number}.fa", "fasta")
sequence = str(fasta.seq)

# start_end = fasta.id[12:]
# end_bp = int(start_end[-6:])
# start_bp =int(start_end[:-7])
# start_bp, end_bp

# start_end = fasta.id[12:]
# end_bp = start_end[-6:]
# start_bp =start_end[:-7]

experiment='chip'

## genomic

raw_peaks = pd.read_csv(f'drive/MyDrive/ML/gcPBM/{protein}/{protein}_{experiment}-exo.gff', sep='\t', header = None)
proc_peaks = raw_peaks[[0,3,4,5,8]]
proc_peaks = proc_peaks.rename(columns={0: "chr", 3: "start", 4:"end", 5:"score", 8:"distance"})
proc_peaks["distance"] = proc_peaks["distance"].str[12:].astype(int)
proc_peaks["start"] = proc_peaks["start"]-proc_peaks["distance"]
for i in range(len(proc_peaks)):
    if proc_peaks["start"][i]<0:
        proc_peaks["start"][i]=0
proc_peaks["end"] = proc_peaks["end"]+proc_peaks["distance"]
proc_peaks = proc_peaks[["chr","start","end","score"]]
min_val = min(proc_peaks["score"])
normalization= max(proc_peaks["score"])-min(proc_peaks["score"])
proc_peaks["score"] = list(((proc_peaks["score"]-min_val)/normalization))

with open(f'drive/MyDrive/ML/gcPBM/{protein}/{protein}_{experiment}-exo.bedgraph','w') as file:
    for i in range(len(proc_peaks)):
        file.write("%s\t" % proc_peaks.iloc[i][0])
        file.write("%s\t" % proc_peaks.iloc[i][1])
        file.write("%s\t" % proc_peaks.iloc[i][2])
        file.write("%s\n" % proc_peaks.iloc[i][3])


chr4 = proc_peaks[proc_peaks['chr']==chr] 
chr4 = chr4.sort_values("start") # [1:12], [13:32]
chr4 = chr4.reset_index(drop=True)
chr4

experiment='pb'

## genomic

raw_peaks = pd.read_csv(f'drive/MyDrive/ML/gcPBM/{protein}/{protein}_{experiment}-exo.gff', sep='\t', header = None)
proc_peaks = raw_peaks[[0,3,4,5,8]]
proc_peaks = proc_peaks.rename(columns={0: "chr", 3: "start", 4:"end", 5:"score", 8:"distance"})
proc_peaks["distance"] = proc_peaks["distance"].str[12:].astype(int)
proc_peaks["start"] = proc_peaks["start"]-proc_peaks["distance"]
for i in range(len(proc_peaks)):
    if proc_peaks["start"][i]<0:
        proc_peaks["start"][i]=0
proc_peaks["end"] = proc_peaks["end"]+proc_peaks["distance"]
proc_peaks = proc_peaks[["chr","start","end","score"]]
min_val = min(proc_peaks["score"])
normalization= max(proc_peaks["score"])-min(proc_peaks["score"])
proc_peaks["score"] = list(((proc_peaks["score"]-min_val)/normalization))

with open(f'drive/MyDrive/ML/gcPBM/{protein}/{protein}_{experiment}-exo.bedgraph','w') as file:
    for i in range(len(proc_peaks)):
        file.write("%s\t" % proc_peaks.iloc[i][0])
        file.write("%s\t" % proc_peaks.iloc[i][1])
        file.write("%s\t" % proc_peaks.iloc[i][2])
        file.write("%s\n" % proc_peaks.iloc[i][3])

chr4 = proc_peaks[proc_peaks['chr']==chr] 
chr4 = chr4.sort_values("start")[34:49] # [12:29], [30:42]
chr4 = chr4.reset_index(drop=True)
chr4

new_peaks = pd.read_csv(f'drive/MyDrive/ML/gcPBM/{protein}/{protein}_rep1_peaks.bed', sep='\t', header = None)
new_peaks = new_peaks[[0,1,2,4]]
new_peaks = new_peaks.rename(columns={0: "chr", 1: "start", 2:"end", 4:"score"})
min_val = min(new_peaks["score"])
normalization= max(new_peaks["score"])-min(new_peaks["score"])
new_peaks["score"] = list(((new_peaks["score"]-min_val)/normalization))

with open(f'drive/MyDrive/ML/gcPBM/{protein}/{protein}_rep1_peaks.bedgraph','w') as file:
    for i in range(len(new_peaks)):
        file.write("%s\t" % new_peaks.iloc[i][0])
        file.write("%s\t" % new_peaks.iloc[i][1])
        file.write("%s\t" % new_peaks.iloc[i][2])
        file.write("%s\n" % new_peaks.iloc[i][3])

chr4 = new_peaks[new_peaks['chr']==chr_r] 
chr4 = chr4.sort_values("start") # [12:29], [30:42]
chr4 = chr4.reset_index(drop=True)
chr4

new_peaks = pd.read_csv(f'drive/MyDrive/ML/gcPBM/{protein}/{protein}_rep2_peaks.bed', sep='\t', header = None)
new_peaks = new_peaks[[0,1,2,4]]
new_peaks = new_peaks.rename(columns={0: "chr", 1: "start", 2:"end", 4:"score"})
min_val = min(new_peaks["score"])
normalization= max(new_peaks["score"])-min(new_peaks["score"])
new_peaks["score"] = list(((new_peaks["score"]-min_val)/normalization))

with open(f'drive/MyDrive/ML/gcPBM/{protein}/{protein}_rep2_peaks.bedgraph','w') as file:
    for i in range(len(new_peaks)):
        file.write("%s\t" % new_peaks.iloc[i][0])
        file.write("%s\t" % new_peaks.iloc[i][1])
        file.write("%s\t" % new_peaks.iloc[i][2])
        file.write("%s\n" % new_peaks.iloc[i][3])

chr4 = new_peaks[new_peaks['chr']==chr_r] 
chr4 = chr4.sort_values("start") # [12:29], [30:42]
chr4 = chr4.reset_index(drop=True)
chr4





# start_time = time.time()

# genomic_tmer_list = [sequence[k:k+30] for k in np.arange(start_bp-000000,end_bp-000000)]
# genomic_tmer_dict = {}

# for tmer in genomic_tmer_list:
#     if tmer not in genomic_tmer_dict.values():
#         df = pd.DataFrame([[tmer,0.1]], columns=["ID_REF", "VALUE"]) #0.1 is a dummy value
#         tmer_dset = Dataset(protein, df, model_target, randomize_fce, chosen_features, score, selected_tetramers)
#         model.predict(tmer_dset)
#         genomic_tmer_dict[tmer] = model.y_pred[0][0]

# print("It took %s seconds" % (time.time() - start_time))
# with open(f'drive/MyDrive/ML/gcPBM/{protein}/{protein}_dictionary_chr5.txt','w') as file:
#     for tmer in list(genomic_tmer_list):
#         file.write("%s\t" % tmer)
#         file.write("%s\n" % genomic_tmer_dict[tmer])

# print('The dictionary has been created')

# data = pd.read_csv(f'drive/MyDrive/ML/gcPBM/{protein}/{protein}_dictionary_long.txt', sep='\t', header=None)
# data = pd.read_csv(f'drive/MyDrive/ML/gcPBM/{protein}/{protein}_dictionary_chr4_400k-550k.txt', sep='\t', header=None)
# data = pd.read_csv(f'drive/MyDrive/ML/gcPBM/{protein}/{protein}_dictionary_chr5_100k-200k.txt', sep='\t', header=None)

#data = pd.read_csv(f'drive/MyDrive/ML/gcPBM/{protein}/{protein}_dictionary_chr6_{number}.txt', sep='\t', header=None)
#data = pd.read_csv(f'drive/MyDrive/ML/gcPBM/{protein}/{protein}_dictionary_chr10_{number}.txt', sep='\t', header=None)
#data = pd.read_csv(f'drive/MyDrive/ML/gcPBM/{protein}/{protein}_dictionary_chr10_relaxed.txt', sep='\t', header=None)
data = pd.read_csv(f'drive/MyDrive/ML/gcPBM/{protein}/{protein}_dictionary_{chr}_{number}.txt', sep='\t', header=None)

genomic_tmer_dict = dict(data.values.tolist())



# make sure you're using the correct dict
print(list(genomic_tmer_dict.keys())[0])
print(sequence[start_bp-start_bp:start_bp-start_bp+30])



# once we have the dictionary, time to plot stuff

# our affinities:
profile_matrix = []
affinities = []

for i in np.arange(start_bp-start_bp,end_bp-29-start_bp):
  start = i
  end = start + 30
  seq = sequence[start:end]
  # if i>start_bp-250000+49720-25 and i<start_bp-250000+49720+5:
  #   print(genomic_tmer_dict[seq], seq)
  profile_matrix.append([seq, genomic_tmer_dict[seq], start+start_bp, end+start_bp])
  aff = 0
  if genomic_tmer_dict[seq] > 0.3:
    aff = genomic_tmer_dict[seq]
  affinities.append(aff)

with open(f'drive/MyDrive/ML/gcPBM/{protein}/affinities_{protein}.bedgraph','w') as file:
  for prof in profile_matrix:
    file.write(f"{chr}\t")
    file.write("%s\t" % prof[2])
    file.write("%s\t" % prof[3])
    file.write("%s\n" % prof[1])

plt.plot(affinities)



experiment='chip'

# Read in Chip-seq peaks from protein
chip_peaks = pd.read_csv(f'drive/MyDrive/ML/gcPBM/{protein}/{protein}_{experiment}-exo.bedgraph', sep='\t', header = None)
chip_peaks.rename({0: "chr",1: "start",2: "end",3: "score"},axis = "columns", inplace = True)

peaks_chr_this = chip_peaks[chip_peaks['chr']==chr] 
filtered_peaks = peaks_chr_this[(peaks_chr_this["start"] > start_bp) & (peaks_chr_this["end"] <= end_bp)]
peaks = filtered_peaks.copy()
peaks = peaks.sort_values("start")
peaks = peaks.reset_index(drop=True)

# Set indices starting at our 0
peaks['start'] = peaks['start'] - start_bp
peaks['end'] = peaks['end'] - start_bp
peaks['length'] = peaks['end']-peaks['start']
peaks = peaks[['start', 'end', 'length','score']]

# create step plot from chip seq data
chip_step_chip = [0] * len(affinities)
i = 0
for j in range(len(affinities)):
  if i < len(peaks):
    if j == peaks["start"].values[i]:
      start = peaks["start"].values[i]
      if (peaks["end"].values[i] >= len(affinities)):
        end = len(affinities)
      else:
        end = peaks["end"].values[i]

      # chip_step[start:end] = [np.max(affinities)] * (end - start) # will input max for length of peak
      chip_step_chip[start:end] = [peaks["score"].values[i]] * (end - start) # will input 1 for length of peak
      i += 1
  else:
    break

peaks

experiment='pb'

# Read in Chip-seq peaks from protein
chip_peaks = pd.read_csv(f'drive/MyDrive/ML/gcPBM/{protein}/{protein}_{experiment}-exo.bedgraph', sep='\t', header = None)
chip_peaks.rename({0: "chr",1: "start",2: "end",3: "score"},axis = "columns", inplace = True)

peaks_chr_this = chip_peaks[chip_peaks['chr']==chr] 
filtered_peaks = peaks_chr_this[(peaks_chr_this["start"] > start_bp) & (peaks_chr_this["end"] <= end_bp)]
peaks = filtered_peaks.copy()
peaks = peaks.sort_values("start")
peaks = peaks.reset_index(drop=True)

# Set indices starting at our 0
peaks['start'] = peaks['start'] - start_bp
peaks['end'] = peaks['end'] - start_bp
peaks['length'] = peaks['end']-peaks['start']
peaks = peaks[['start', 'end', 'length','score']]

# create step plot from chip seq data
chip_step_pb = [0] * len(affinities)
i = 0
for j in range(len(affinities)):
  if i < len(peaks):
    if j == peaks["start"].values[i]:
      start = peaks["start"].values[i]
      if (peaks["end"].values[i] >= len(affinities)):
        end = len(affinities)
      else:
        end = peaks["end"].values[i]

      # chip_step[start:end] = [np.max(affinities)] * (end - start) # will input max for length of peak
      chip_step_pb[start:end] = [peaks["score"].values[i]] * (end - start) # will input 1 for length of peak
      i += 1
  else:
    break



# new peaks

# Read in Chip-seq peaks from protein
chip_peaks = pd.read_csv(f'drive/MyDrive/ML/gcPBM/{protein}/{protein}_rep1_peaks.bedgraph', sep='\t', header = None)
chip_peaks.rename({0: "chr",1: "start",2: "end",3: "score"},axis = "columns", inplace = True)

peaks_chr_this = chip_peaks[chip_peaks['chr']==chr_r] 
filtered_peaks = peaks_chr_this[(peaks_chr_this["start"] > start_bp) & (peaks_chr_this["end"] <= end_bp)]
peaks = filtered_peaks.copy()
peaks = peaks.sort_values("start")
peaks = peaks.reset_index(drop=True)

# Set indices starting at our 0
peaks['start'] = peaks['start'] - start_bp
peaks['end'] = peaks['end'] - start_bp
peaks['length'] = peaks['end']-peaks['start']
peaks = peaks[['start', 'end', 'length','score']]

# create step plot from chip seq data
chip_step_new = [0] * len(affinities)
i = 0
for j in range(len(affinities)):
  if i < len(peaks):
    if j == peaks["start"].values[i]:
      start = peaks["start"].values[i]
      if (peaks["end"].values[i] >= len(affinities)):
        end = len(affinities)
      else:
        end = peaks["end"].values[i]

      # chip_step[start:end] = [np.max(affinities)] * (end - start) # will input max for length of peak
      chip_step_new[start:end] = [peaks["score"].values[i]] * (end - start) # will input 1 for length of peak
      i += 1
  else:
    break

# new peaks second replicate

# Read in Chip-seq peaks from protein
chip_peaks = pd.read_csv(f'drive/MyDrive/ML/gcPBM/{protein}/{protein}_rep2_peaks.bedgraph', sep='\t', header = None)
chip_peaks.rename({0: "chr",1: "start",2: "end",3: "score"},axis = "columns", inplace = True)

peaks_chr_this = chip_peaks[chip_peaks['chr']==chr_r] 
filtered_peaks = peaks_chr_this[(peaks_chr_this["start"] > start_bp) & (peaks_chr_this["end"] <= end_bp)]
peaks = filtered_peaks.copy()
peaks = peaks.sort_values("start")
peaks = peaks.reset_index(drop=True)

# Set indices starting at our 0
peaks['start'] = peaks['start'] - start_bp
peaks['end'] = peaks['end'] - start_bp
peaks['length'] = peaks['end']-peaks['start']
peaks = peaks[['start', 'end', 'length','score']]

# create step plot from chip seq data
chip_step_new_2 = [0] * len(affinities)
i = 0
for j in range(len(affinities)):
  if i < len(peaks):
    if j == peaks["start"].values[i]:
      start = peaks["start"].values[i]
      if (peaks["end"].values[i] >= len(affinities)):
        end = len(affinities)
      else:
        end = peaks["end"].values[i]

      # chip_step[start:end] = [np.max(affinities)] * (end - start) # will input max for length of peak
      chip_step_new_2[start:end] = [peaks["score"].values[i]] * (end - start) # will input 1 for length of peak
      i += 1
  else:
    break

peaks

# nucleosomes

# nuc_peaks = pd.read_csv(f'drive/MyDrive/ML/gcPBM/yeast_nucleosomes.gff', sep='\t', header = None)
# nuc_peaks.rename({0: "chr",3: "start",4: "end",5: "score"},axis = "columns", inplace = True)
# nuc_chr_this = nuc_peaks.loc[nuc_peaks['chr'] == 'chrX']
# filtered_nuc = nuc_chr_this[(nuc_chr_this["start"] > start_bp) & (nuc_chr_this["end"] <= end_bp)]
# nucleosomes = filtered_nuc.copy()

# Set indices starting at our 0
# nucleosomes['start'] = nucleosomes['start'] - start_bp
# nucleosomes['end'] = nucleosomes['end'] - start_bp
# nucleosomes['length'] = nucleosomes['end']-nucleosomes['start']
# nucleosomes = nucleosomes[['start', 'end', 'length','score']]

# create step plot from nuc data
# nuc_step = [0] * len(affinities)
# i = 0
# for j in range(len(affinities)):
#   if i < len(nucleosomes):
#     if j == nucleosomes["start"].values[i]:
#       start = nucleosomes["start"].values[i]
#       if (nucleosomes["end"].values[i] >= len(affinities)):
#         end = len(affinities)
#       else:
#         end = nucleosomes["end"].values[i]

      # nuc_step[start:end] = [np.max(affinities)] * (end - start) # will input max for length of peak
#       nuc_step[start:end] = [1] * (end - start) # will input 1 for length of peak
#       i += 1
#   else:
#     break

# new nucleosomes

nuc_peaks = pd.read_csv(f'drive/MyDrive/ML/gcPBM/avg_nuc_calls.txt', sep=',')
nuc_peaks = nuc_peaks[['chr', 'start', 'end', 'score']]
nuc_chr_this = nuc_peaks.loc[nuc_peaks['chr'] == chr_r]
filtered_nuc = nuc_chr_this[(nuc_chr_this["start"] > start_bp) & (nuc_chr_this["end"] <= end_bp)]
nucleosomes = filtered_nuc.copy()

# Set indices starting at our 0
nucleosomes['start'] = nucleosomes['start'] - start_bp
nucleosomes['end'] = nucleosomes['end'] - start_bp
nucleosomes['length'] = nucleosomes['end']-nucleosomes['start']
nucleosomes = nucleosomes[['start', 'end', 'length','score']]

# create step plot from nuc data
nuc_step = [0] * len(affinities)
i = 0
for j in range(len(affinities)):
  if i < len(nucleosomes):
    if j == nucleosomes["start"].values[i]:
      start = nucleosomes["start"].values[i]
      if (nucleosomes["end"].values[i] >= len(affinities)):
        end = len(affinities)
      else:
        end = nucleosomes["end"].values[i]

      # nuc_step[start:end] = [np.max(affinities)] * (end - start) # will input max for length of peak
      nuc_step[start:end] = [nucleosomes['score'].values[i]] * (end - start) # will input 1 for length of peak
      i += 1
  else:
    break

# polymerase

pol_peaks_1 = pd.read_csv(f'drive/MyDrive/ML/poly/MACS2_ph0_R1_narrowPeaks.bed', sep='\t', header=None)
pol_peaks_2 = pd.read_csv(f'drive/MyDrive/ML/poly/MACS2_ph0_R2_narrowPeaks.bed', sep='\t', header=None)
pol_peaks_1.rename({0: "chr",1: "start",2: "end",4: "score"},axis = "columns", inplace = True)
pol_peaks_2.rename({0: "chr",1: "start",2: "end",4: "score"},axis = "columns", inplace = True)
pol_peaks_1 = pol_peaks_1[['chr', 'start', 'end', 'score']]
pol_peaks_2 = pol_peaks_2[['chr', 'start', 'end', 'score']]

pol_peaks_this_1 = pol_peaks_1[pol_peaks_1['chr']==chr_r] 
filtered_peaks_1 = pol_peaks_this_1[(pol_peaks_this_1["start"] > start_bp) & (pol_peaks_this_1["end"] <= end_bp)]
peaks_1 = filtered_peaks_1.copy()
peaks_1 = peaks_1.sort_values("start")
peaks_1 = peaks_1.reset_index(drop=True)

pol_peaks_this_2 = pol_peaks_2[pol_peaks_2['chr']==chr_r] 
filtered_peaks_2 = pol_peaks_this_2[(pol_peaks_this_2["start"] > start_bp) & (pol_peaks_this_2["end"] <= end_bp)]
peaks_2 = filtered_peaks_2.copy()
peaks_2 = peaks_2.sort_values("start")
peaks_2 = peaks_2.reset_index(drop=True)

# Set indices starting at our 0
peaks_1['start'] = peaks_1['start'] - start_bp
peaks_1['end'] = peaks_1['end'] - start_bp
peaks_1['length'] = peaks_1['end']-peaks_1['start']
peaks_1 = peaks_1[['start', 'end', 'length','score']]
normalization= max(peaks_1["score"])-min(peaks_1["score"])
min_val = min(peaks_1["score"])
peaks_1["score"] = list(((peaks_1["score"]-min_val)/normalization))

peaks_2['start'] = peaks_2['start'] - start_bp
peaks_2['end'] = peaks_2['end'] - start_bp
peaks_2['length'] = peaks_2['end']-peaks_2['start']
peaks_2 = peaks_2[['start', 'end', 'length','score']]
normalization= max(peaks_2["score"])-min(peaks_2["score"])
min_val = min(peaks_2["score"])
peaks_2["score"] = list(((peaks_2["score"]-min_val)/normalization))

# create step plot from nuc data
peak_step_1 = [0] * len(affinities)
i = 0
for j in range(len(affinities)):
  if i < len(peaks_1):
    if j == peaks_1["start"].values[i]:
      start = peaks_1["start"].values[i]
      if (peaks_1["end"].values[i] >= len(affinities)):
        end = len(affinities)
      else:
        end = peaks_1["end"].values[i]

      # nuc_step[start:end] = [np.max(affinities)] * (end - start) # will input max for length of peak
      peak_step_1[start:end] = [peaks_1['score'].values[i]] * (end - start) # will input 1 for length of peak
      i += 1
  else:
    break

# create step plot from nuc data
peak_step_2 = [0] * len(affinities)
i = 0
for j in range(len(affinities)):
  if i < len(peaks_2):
    if j == peaks_2["start"].values[i]:
      start = peaks_2["start"].values[i]
      if (peaks_2["end"].values[i] >= len(affinities)):
        end = len(affinities)
      else:
        end = peaks_2["end"].values[i]

      # nuc_step[start:end] = [np.max(affinities)] * (end - start) # will input max for length of peak
      peak_step_2[start:end] = [peaks_2['score'].values[i]] * (end - start) # will input 1 for length of peak
      i += 1
  else:
    break





window_start = 56000
window_end =   59000

fig = plt.gcf()
fig.set_size_inches(24, 10)
plt.plot(nuc_step[window_start:window_end], color='red')
plt.plot(affinities[window_start:window_end], color='blue')
plt.plot(chip_step_pb[window_start:window_end],color='green')
plt.plot(chip_step_chip[window_start:window_end],color='orange')
plt.plot(chip_step_new[window_start:window_end],color='purple')
plt.plot(chip_step_new_2[window_start:window_end],color='hotpink')

plt.plot(peak_step_1[window_start:window_end],color='black')
plt.plot(peak_step_2[window_start:window_end],color='darkgray')

plt.xlabel('Base')
plt.ylabel('Relative Intensity')
plt.title(f'{chr} - {number-1}00k-{number}00k subframe')
labels = [ "Nucleosomes", 'ML Prediction','PB-exo', 'ChIP-exo', 'New data - Rep1', 'New data - Rep2', 'Poly 1', 'Poly 2']
plt.legend(labels, loc="best")



tipo = 'PB'
window_start = 56000
window_end =   59000
width = 24
height = 4

fig = plt.gcf()
fig.set_size_inches(width, height)
plt.plot(chip_step_pb[window_start:window_end],color='green')
plt.plot(chip_step_chip[window_start:window_end],color='black')
plt.plot(peak_step_1[window_start:window_end],color='darkorange')
plt.plot(chip_step_new[window_start:window_end],color='black')
plt.plot(chip_step_new_2[window_start:window_end],color='black')
plt.plot(peak_step_2[window_start:window_end],color='darkorange')
plt.ylabel('Relative Intensity')
plt.title(f'{chr_r} sequence')
labels = ['PB-exo', 'ChIP-exo', 'Polymerase']
plt.legend(labels, loc="best", prop={'size': 20})
plt.show()
fig.savefig(f'drive/MyDrive/ML/gcPBM/figures/{tipo}/pbchippoly.png')
plt.clf()


fig = plt.gcf()
fig.set_size_inches(width, height)
plt.plot(affinities[window_start:window_end], color='blue')
plt.ylabel('Relative Intensity')
plt.title(f'{chr_r} sequence')
labels = ['ML Prediction']
plt.legend(labels, loc="best", prop={'size': 20})
plt.show()
fig.savefig(f'drive/MyDrive/ML/gcPBM/figures/{tipo}/prediction.png')
plt.clf()


fig = plt.gcf()
fig.set_size_inches(width, height)
plt.plot(nuc_step[window_start:window_end], color='red')
plt.ylabel('Relative Intensity')
plt.title(f'{chr_r} sequence')
labels = [ "Nucleosomes"]
plt.legend(labels, loc="best", prop={'size': 20})
plt.show()
fig.savefig(f'drive/MyDrive/ML/gcPBM/figures/{tipo}/nucleosomes.png')

window_start = 36 * 1000
window_end =   37 * 1000


fig = plt.gcf()
fig.set_size_inches(16, 8)
plt.plot(affinities[window_start:window_end])
plt.plot(nuc_step[window_start:window_end], color='red')
plt.plot(chip_step_pb[window_start:window_end],color='green')
plt.plot(chip_step_chip[window_start:window_end],color='orange')
plt.plot(chip_step_new[window_start:window_end],color='purple')
plt.plot(chip_step_new_2[window_start:window_end],color='hotpink')
labels = ['ML Prediction', 'Nucleosomes', 'PB-exo', 'ChIP-exo (1st Study)', 'ChIP-exo (2nd Study)- Rep1', 'ChIP-exo (2nd Study)- Rep2']
plt.legend(labels, loc="best")

counts = pd.read_csv(f'drive/MyDrive/ML/gcPBM/{protein}/genome_counts.txt', sep='\t')
counts.sum()

121/69

tfbs = sequence[267900-start_bp+0:267900-start_bp+200]

scores = []
for i in np.arange(len(tfbs)-5):
  score = 0
  print(i,tfbs[i:i+6])
  for j in np.arange(0,len(freq_matrix)):
      a = translate[tfbs[i+j]]
      score = score + freq_matrix[j][a]
  scores.append(score)

[(tfbs[i:i+6],i,scores[i]) for i in range(len(scores)) if scores[i]>=4]

top_feat = pd.read_csv(f'drive/MyDrive/ML/top_features.txt', sep='\t', header = None)
freqs = top_feat[0].value_counts()

fig = plt.gcf()
fig.set_size_inches(16, 8)
plt.yticks(range(0,230,20))
freqs.plot(kind='barh')
print(freqs)

top_feat = pd.read_csv(f'drive/MyDrive/ML/top_features.txt', sep='\t', header = None)
top_feat_no_pres = pd.read_csv(f'drive/MyDrive/ML/top_features_no_pres.txt', sep='\t', header = None)
top_feat_no_pres_rand = pd.read_csv(f'drive/MyDrive/ML/top_features_rand_no_pres.txt', sep='\t', header = None)

top_feat.rename({0: "Presence"},axis = "columns", inplace = True)
freqs = top_feat['Presence'].value_counts()
top_feat_no_pres.rename({0: "No presence"},axis = "columns", inplace = True)
freqs_np = top_feat_no_pres["No presence"].value_counts()
top_feat_no_pres_rand.rename({0: "No presence - Randomized"},axis = "columns", inplace = True)
freqs_npr = top_feat_no_pres_rand["No presence - Randomized"].value_counts()

df = pd.concat([freqs,freqs_np,freqs_npr],axis=1)
df = df.drop(['Presence'])
print(df)

ax =df.plot(kind='line',xticks = range(len(df)),figsize=(16,8),rot=90)
ax.set_xticklabels(list(df.index))
ax.set_ylabel('Frequency')

shape_feats = ['SHIFT', 'SLIDE', 'RISE', 'TILT', 'ROLL', 'TWIST']

base_features = pd.read_csv(f'drive/MyDrive/ML/gcPBM/feature_study/all_shape_feats.txt', sep='\t', header = None)
no_SHIFT = pd.read_csv(f'drive/MyDrive/ML/gcPBM/feature_study/SHIFT_feature.txt', sep='\t', header = None)
no_SLIDE = pd.read_csv(f'drive/MyDrive/ML/gcPBM/feature_study/SLIDE_feature.txt', sep='\t', header = None)
no_RISE = pd.read_csv(f'drive/MyDrive/ML/gcPBM/feature_study/RISE_feature.txt', sep='\t', header = None)
no_TILT = pd.read_csv(f'drive/MyDrive/ML/gcPBM/feature_study/TILT_feature.txt', sep='\t', header = None)
no_ROLL = pd.read_csv(f'drive/MyDrive/ML/gcPBM/feature_study/ROLL_feature.txt', sep='\t', header = None)
no_TWIST = pd.read_csv(f'drive/MyDrive/ML/gcPBM/feature_study/TWIST_feature.txt', sep='\t', header = None)

freqs_base = base_features[0].value_counts(normalize=True)
freqs_SHIFT = no_SHIFT[0].value_counts(normalize=True)
freqs_SLIDE = no_SLIDE[0].value_counts(normalize=True)
freqs_RISE = no_RISE[0].value_counts(normalize=True)
freqs_TILT = no_TILT[0].value_counts(normalize=True)
freqs_ROLL = no_ROLL[0].value_counts(normalize=True)
freqs_TWIST = no_TWIST[0].value_counts(normalize=True)

freqs = []
freqs.append(freqs_SHIFT)
freqs.append(freqs_SLIDE)
freqs.append(freqs_RISE)
freqs.append(freqs_TILT)
freqs.append(freqs_ROLL)
freqs.append(freqs_TWIST)

df = pd.DataFrame(index=(freqs_base.index.sort_values()))
for (freq,feat) in zip(freqs,shape_feats):
  df['No '+feat] = (freq-freqs_base)/freqs_base*100

df

df = df[['No RISE','No ROLL', 'No SHIFT', 'No SLIDE', 'No TILT', 'No TWIST']]
base_features

import seaborn as sns; sns.set()
#max_val = max(abs(df.max().max()), abs(df.min().min()))
t = df.transpose()
f, ax = plt.subplots(figsize=(24, 10))
ax = sns.heatmap(t,  annot=True, cmap="RdBu", center=0)



# correlations

# feature-dependent correlations
shape_df = pd.read_csv(f'drive/MyDrive/ML/gcPBM/feature_study/correlations.txt', sep='\t',header=None)
prots = [shape_df[0][j] for j in range(len(shape_df)) if shape_df[1][j] == 'SHIFT']
SHIFT_df = [shape_df[2][j] for j in range(len(shape_df)) if shape_df[1][j] == 'SHIFT']
SLIDE_df = [shape_df[2][j] for j in range(len(shape_df)) if shape_df[1][j] == 'SLIDE']
RISE_df = [shape_df[2][j] for j in range(len(shape_df)) if shape_df[1][j] == 'RISE']
TILT_df = [shape_df[2][j] for j in range(len(shape_df)) if shape_df[1][j] == 'TILT']
ROLL_df = [shape_df[2][j] for j in range(len(shape_df)) if shape_df[1][j] == 'ROLL']
TWIST_df = [shape_df[2][j] for j in range(len(shape_df)) if shape_df[1][j] == 'TWIST']

df = pd.DataFrame(list(zip(prots,SHIFT_df, SLIDE_df, RISE_df, TILT_df, ROLL_df, TWIST_df)))
df = df.set_index(0)

fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(20, 20))
#medianprops = dict(linestyle='-', linewidth=4, color='red')

#t = df.transpose()
#labels = list(t.index)
df.rename({1:'No SHIFT', 2:'No SLIDE', 3:'No RISE', 4:'No TILT', 5:'No ROLL', 6:'No TWIST'},axis = "columns", inplace = True)

#labels = ['No SHIFT', 'No SLIDE', 'No RISE', 'No TILT', 'No ROLL', 'No TWIST']
#green_diamond = dict(markerfacecolor='g', marker='D')
#bplot = ax.boxplot(t,  flierprops=green_diamond,   patch_artist=True, medianprops=medianprops,labels=labels)

#ax.set_yticks(np.arange(-0.7,1,0.02))
#ax.set_title('Correlations')
#ax.yaxis.grid(True, color='lightgrey')
#ax.set_facecolor('white')

#print(df.mean())
print(df)
vals, names, xs = [],[],[]
for i, col in enumerate(df.columns):
    vals.append(df[col].values)
    names.append(col)
    xs.append(np.random.normal(i + 1, 0.04, df[col].values.shape[0]))  # adds jitter to the data points - can be adjusted

#labeling = list(df.index)
df.boxplot(ax = ax,showfliers=False, showmeans=True,fontsize = 40,rot=45,medianprops = dict(linestyle='-', linewidth=3, color='black'), 
                   meanprops={"marker":"o","markerfacecolor":"black", "markeredgecolor":"black","markersize":"10"})
palette = ['r', 'g', 'b', 'y', 'hotpink', 'darkorange']
for x, val, c in zip(xs, vals, palette):
    plt.scatter(x, val, alpha=0.8, color=c)

ax.set_yticks(np.arange(-0.8,1.1,0.1))
ax.set_ylabel('$R^2$',size=50)
ax.yaxis.grid(True, color='lightgrey')
ax.set_facecolor('white')



pres_corr = pd.read_csv(f'drive/MyDrive/ML/electro_corr.txt', sep='\t')
pres_corr = pres_corr.set_index('protein')
pres_corr = pres_corr[['0vs5']]
pres_corr.rename({"0vs5": "Presence"},axis = "columns", inplace = True)

np_corr = pd.read_csv(f'drive/MyDrive/ML/no_pres_corr.txt', sep='\t',header=None)
np_corr.rename({0: "protein", 2:"0vs5"},axis = "columns", inplace = True)
np_corr = np_corr.set_index('protein')
np_corr = np_corr[['0vs5']]
np_corr.rename({"0vs5": "No presence"},axis = "columns", inplace = True)

npr_corr = pd.read_csv(f'drive/MyDrive/ML/rand_no_pres_corr.txt', sep='\t',header=None)
npr_corr.rename({0: "protein", 2:"0vs5"},axis = "columns", inplace = True)
npr_corr = npr_corr.set_index('protein')
npr_corr = npr_corr[['0vs5']]
npr_corr.rename({"0vs5": "No presence-Randomized"},axis = "columns", inplace = True)


df = pd.concat([pres_corr,np_corr,npr_corr],axis=1)

ax = df.plot(kind='line',xticks = range(len(df)), yticks=np.arange(0,1,0.1),figsize=(20,8))
ax.set_xticklabels(list(df.index))
ax.set_ylabel('R²')

fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10, 10))

t = df[['Presence', 'No presence']]

vals, names, xs = [],[],[]
for i, col in enumerate(t.columns):
    vals.append(t[col].values)
    names.append(col)
    xs.append(np.random.normal(i + 1, 0.04, t[col].values.shape[0]))  # adds jitter to the data points - can be adjusted



labeling = list(t.index)
t.boxplot(ax = ax,showfliers=False, showmeans=True,fontsize = 30, rot=45,medianprops = dict(linestyle='-', linewidth=3, color='black'), 
                   meanprops={"marker":"o","markerfacecolor":"black", "markeredgecolor":"black","markersize":"10"})
palette = ['r', 'b', 'b', 'y']
for x, val, c in zip(xs, vals, palette):
    plt.scatter(x, val, alpha=0.8, color=c)

ax.set_yticks(np.arange(-0.2,1.1,0.1))

ax.yaxis.grid(True, color='lightgrey')
ax.set_ylabel('HT-SELEX $R^2$',size=40)
ax.set_facecolor('white')

pres_corr = pd.read_csv(f'drive/MyDrive/ML/electro_corr_linear.txt', sep=' ',header=None)
pres_corr.rename({0: "protein", 2:"0vs5"},axis = "columns", inplace = True)
pres_corr = pres_corr.set_index('protein')
pres_corr = pres_corr[['0vs5']]
pres_corr.rename({"0vs5": "Presence"},axis = "columns", inplace = True)

np_corr = pd.read_csv(f'drive/MyDrive/ML/no_pres_corr_linear.txt', sep='\t',header=None)
np_corr.rename({0: "protein", 2:"0vs5"},axis = "columns", inplace = True)
np_corr = np_corr.set_index('protein')
np_corr = np_corr[['0vs5']]
np_corr.rename({"0vs5": "No presence"},axis = "columns", inplace = True)

npr_corr = pd.read_csv(f'drive/MyDrive/ML/rand_no_pres_corr_linear.txt', sep='\t',header=None)
npr_corr.rename({0: "protein", 2:"0vs5"},axis = "columns", inplace = True)
npr_corr = npr_corr.set_index('protein')
npr_corr = npr_corr[['0vs5']]
npr_corr.rename({"0vs5": "No presence-Randomized"},axis = "columns", inplace = True)

df = pd.concat([pres_corr,np_corr,npr_corr],axis=1)

ax = df.plot(kind='line',xticks = range(len(df)), yticks=np.arange(-0.6,1,0.1),figsize=(20,8))
ax.set_xticklabels(list(df.index))
ax.set_ylabel('R²')

fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10, 8))
medianprops = dict(linestyle='-', linewidth=4, color='red')

t = df.transpose()
labels = list(t.index)
green_diamond = dict(markerfacecolor='g', marker='D')
bplot = ax.boxplot(t,  flierprops=green_diamond,   patch_artist=True, medianprops=medianprops,labels=labels)

ax.set_title('Correlations')
ax.yaxis.grid(True, color='lightgrey')
ax.set_facecolor('white')



selex_1 = pd.read_csv(f'drive/MyDrive/ML/correlations/electro_corr_clear.txt', sep=' ')
selex_2 = pd.read_csv(f'drive/MyDrive/ML/correlations/electro_corr_2.txt', sep='\t')
upbm = pd.read_csv(f'drive/MyDrive/ML/correlations/corr_ws_2.txt', sep='\t')
new_upbm = pd.read_csv(f'drive/MyDrive/ML/correlations/new_upbm_correlations.txt', sep='\t')
ratios  = pd.read_csv(f'drive/MyDrive/ML/correlations/ratios.txt', sep=' ',header=None)

upbm = pd.concat([upbm,ratios],axis=1)
l = []
for i in range(len(upbm[0])):
    if upbm[0][i] > 0.1:
        l.append(i)
upbm.index = [i for i in range(len(upbm))]
u_pbm = upbm.drop(l)[['PROTEINS','FAMILY', 'PRESENCE', 'AVG+DIAG', 'P+A+D']]

new_upbm = pd.concat([new_upbm,ratios],axis=1)
ll = []
for i in range(len(new_upbm[0])):
    if new_upbm[0][i] > 0.1:
        ll.append(i)
new_upbm.index =  [i for i in range(len(new_upbm))]
new_u_pbm = new_upbm.drop(ll)[['protein', 'correlation']]

#new_u_pbm = new_upbm.drop(l)[['protein', 'presence','electro','shape']]

u_pbm.rename({'PROTEINS': "protein", 'P+A+D':'uPBM'},axis = "columns", inplace = True)

selex_1.index = selex_1['protein']
selex_2.index = selex_2['protein']
selex_1 = selex_1['0vs3']
selex_2 = selex_2['0vs5']
selex_full = pd.concat([selex_1,selex_2])

u_pbm.index= u_pbm['protein']
new_u_pbm.index= u_pbm['protein']
u_pbm = u_pbm['uPBM']

gc_pbm = {'gcPBM': [0.905, 0.951, 0.922], 'protein': ['myc', 'mad1', 'max']}
gc_pbm = pd.DataFrame(data=gc_pbm)
gc_pbm.index = gc_pbm['protein']
gc_pbm = gc_pbm['gcPBM']

t = pd.concat([selex_full,u_pbm,gc_pbm],axis=1)
t.rename({0: "HT-SELEX"},axis = "columns", inplace = True)
t = t.sort_index()

import seaborn as sns; sns.set()
#t = df.transpose()
f, ax = plt.subplots(figsize=(5, 10))
ax = sns.heatmap(t[0:30],  annot=True, cmap="RdBu", center=0)
f, ax = plt.subplots(figsize=(5, 10))
ax = sns.heatmap(t[30:60],  annot=True, cmap="RdBu", center=0)
f, ax = plt.subplots(figsize=(5, 10))
ax = sns.heatmap(t[60:90],  annot=True, cmap="RdBu", center=0)
f, ax = plt.subplots(figsize=(5, 10))
ax = sns.heatmap(t[90:],  annot=True, cmap="RdBu", center=0)

t.to_csv(r'drive/MyDrive/ML/correlations/all_datasets.txt', sep='\t', mode='w')

(t['uPBM'].sort_values()[:34]).to_csv(r'drive/MyDrive/ML/correlations/upbm_filtered.txt', sep='\t', mode='w')

t.mean(), t.std()

upbm = pd.read_csv(f'drive/MyDrive/ML/correlations/corr_ws_2.txt', sep='\t')
new_upbm = pd.read_csv(f'drive/MyDrive/ML/correlations/new_upbm_correlations.txt', sep='\t')
ratios  = pd.read_csv(f'drive/MyDrive/ML/correlations/ratios.txt', sep=' ',header=None)
features = pd.read_csv(f'drive/MyDrive/ML/correlations/upbm_features.txt', sep='\t')

upbm = pd.concat([upbm[['PROTEINS','P+A+D']],ratios[0]],axis=1)
upbm.index = upbm['PROTEINS']
new_upbm.index = new_upbm['protein']
features.index = features['protein']
new_upbm = new_upbm.drop('sp140')
features = features.drop('sp140')

DF = pd.concat([upbm,new_upbm,features],axis=1)

l = []
for i in range(len(DF[0])):
    if DF[0][i] > 0.1:
        l.append(i)
DF.index = [i for i in range(len(DF))]
df = DF.drop(l)[['PROTEINS','P+A+D','correlation', 0, 'presence','electro','shape']]
f_upbm = list(df.mean()[['presence','electro','shape']])

gcpbm_feat = pd.read_csv(f'drive/MyDrive/ML/gcPBM/feature_study/gcpbm_features.txt', sep='\t', header= None)
gcpbm_feat.rename({1: "presence", 2:'electro', 3:'shape'},axis = "columns", inplace = True)
f_gcpbm = list(gcpbm_feat.mean())

selex_feat = pd.read_csv(f'drive/MyDrive/ML/correlations/htselex_features.txt', sep='\t', header= None)
selex_feat.rename({1: "presence", 2:'electro', 3:'shape'},axis = "columns", inplace = True)
f_selex = list(selex_feat.mean())





len(selex_feat)

print(f_upbm,f_gcpbm,f_selex)

labels = ["Presence", "Electrostatic", "Shape"]
#colors = ['red','darkorange','yellow','hotpink','green','blue','cyan','lightgreen','darkgray','purple','lightblue','olive']
colors = ['red', 'cyan', 'purple']

fig1, ax1 = plt.subplots()
fig1.set_size_inches(15, 15)
ax1.pie(f_upbm, labels = labels,startangle=0,  colors=colors,  autopct='%1.1f%%',textprops={'fontsize': 25})
ax1.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
plt.show()

fig1, ax1 = plt.subplots()
fig1.set_size_inches(15, 15)
ax1.pie(f_gcpbm, labels = labels,startangle=0,  colors=colors,  autopct='%1.1f%%',textprops={'fontsize': 25})
ax1.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
plt.show()

fig1, ax1 = plt.subplots()
fig1.set_size_inches(15, 15)
ax1.pie(f_selex, labels = labels,startangle=0,  colors=colors,  autopct='%1.1f%%',textprops={'fontsize': 25})
ax1.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
plt.show()

(146-20-17-15)/1326

selex_1 = pd.read_csv(f'drive/MyDrive/ML/correlations/electro_corr.txt', sep=' ')
selex_2 = pd.read_csv(f'drive/MyDrive/ML/correlations/electro_corr_2.txt', sep='\t')
upbm = pd.read_csv(f'drive/MyDrive/ML/correlations/corr_ws.txt', sep='\t')
ratios  = pd.read_csv(f'drive/MyDrive/ML/correlations/ratios.txt', sep=' ',header=None)

upbm = pd.concat([upbm,ratios],axis=1)
l = []
for i in range(len(upbm[0])):
    if upbm[0][i] > 0.1:
        l.append(i)
upbm.index = [i for i in range(len(upbm))]
u_pbm = upbm.drop(l)[['PROTEINS','FAMILY', 'PRESENCE', 'AVG+DIAG', 'P+A+D']]

selex_1 = selex_1['0vs3']
selex_2 = selex_2['0vs5']
u_pbm = u_pbm['P+A+D']
gc_pbm = {'gcPBM': [0.905, 0.951, 0.922]}

gc_pbm = pd.DataFrame(data=gc_pbm)
selex_1.mean(), selex_2.mean(), u_pbm.mean(), gc_pbm.mean()

selex_2.mean(),selex_2.std()

t

selex = pd.concat([selex_1,selex_2],axis=0)
#change name col
selex = selex.reset_index(drop=True)
t = pd.concat([selex,u_pbm,gc_pbm,filtered['Filtered']],axis=1)
t.rename({0: "HT-SELEX", 'P+A+D': 'uPBM','Filtered':'HT-SELEX Filtered'},axis = "columns", inplace = True)
t.dropna()

fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(24, 20))

vals, names, xs = [],[],[]
for i, col in enumerate(t.columns):
    vals.append(t[col].values)
    names.append(col)
    xs.append(np.random.normal(i + 1, 0.04, t[col].values.shape[0]))  # adds jitter to the data points - can be adjusted

#t_1 = (t[['selex']].dropna()).transpose()
#t_2 = (t[['P+A+D']].dropna()).transpose()
#t_3 = (t[['gcPBM']].dropna()).transpose()
#positions = np.arange(3)+1
labeling = ['HT-SELEX', 'uPBM', 'gcPBM','HT-SELEX Filtered']
bplot = t.boxplot(labels=labeling, ax = ax,showfliers=False, showmeans=True,fontsize = 40, rot=45,medianprops = dict(linestyle='-', linewidth=3, color='black'), 
                   meanprops={"marker":"o","markerfacecolor":"black", "markeredgecolor":"black","markersize":"10"})
#bplot = ax.boxplot(t_1, showfliers=False, showmeans=True,medianprops = dict(linestyle='-', linewidth=3, color='black'), 
#                   meanprops={"marker":"o","markerfacecolor":"black", "markeredgecolor":"black","markersize":"10"})
#bplot = ax.boxplot(t_2, showfliers=False,showmeans=True,medianprops = dict(linestyle='-', linewidth=3, color='black'), 
#                   meanprops={"marker":"o","markerfacecolor":"black", "markeredgecolor":"black","markersize":"10"})
#bplot = ax.boxplot(t_3, showfliers=False,showmeans=True,medianprops = dict(linestyle='-', linewidth=3, color='black'), 
#                   meanprops={"marker":"o","markerfacecolor":"black", "markeredgecolor":"black","markersize":"10"})

palette = ['r', 'g', 'b', 'y']
for x, val, c in zip(xs, vals, palette):
    plt.scatter(x, val, alpha=0.8, color=c)

ax.set_yticks(np.arange(-0.1,1.01,0.1))
ax.yaxis.grid(True, color='lightgrey')
ax.set_ylabel('$R^2$',size=40)
ax.set_facecolor('white')



# only rohs's
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(24, 20))

vals, names, xs = [],[],[]
rohs_x, rohs_y = [], []
for i, col in enumerate(t.columns):
    vals.append(t[col].values)
    names.append(col)
    xs.append(np.random.normal(i + 1, 0.04, t[col].values.shape[0]))  # adds jitter to the data points - can be adjusted
    rohs_x.append(np.arange(i+1-0.23,i+1+0.24,0.02))

rohs_y.append([0.67]*len(rohs_x[0]))
rohs_y.append([0.46]*len(rohs_x[0]))
rohs_y.append([0.93]*len(rohs_x[0]))
rohs_y.append([0.67]*len(rohs_x[0]))

dl_y = [0.64]*len(rohs_x[0])


#t_1 = (t[['selex']].dropna()).transpose()
#t_2 = (t[['P+A+D']].dropna()).transpose()
#t_3 = (t[['gcPBM']].dropna()).transpose()
#positions = np.arange(3)+1
labeling = ['HT-SELEX', 'uPBM', 'gcPBM','HT-SELEX Filtered']
bplot = t.boxplot(labels=labeling, ax = ax,showfliers=False, showmeans=True,fontsize = 40, rot=45,medianprops = dict(linestyle='-', linewidth=3, color='black'), 
                   meanprops={"marker":"o","markerfacecolor":"black", "markeredgecolor":"black","markersize":"10"})
#bplot = ax.boxplot(t_1, showfliers=False, showmeans=True,medianprops = dict(linestyle='-', linewidth=3, color='black'), 
#                   meanprops={"marker":"o","markerfacecolor":"black", "markeredgecolor":"black","markersize":"10"})
#bplot = ax.boxplot(t_2, showfliers=False,showmeans=True,medianprops = dict(linestyle='-', linewidth=3, color='black'), 
#                   meanprops={"marker":"o","markerfacecolor":"black", "markeredgecolor":"black","markersize":"10"})
#bplot = ax.boxplot(t_3, showfliers=False,showmeans=True,medianprops = dict(linestyle='-', linewidth=3, color='black'), 
#                   meanprops={"marker":"o","markerfacecolor":"black", "markeredgecolor":"black","markersize":"10"})

palette = ['r', 'g', 'b', 'y']
for x, val, c, r_x, r_y in zip(xs, vals, palette, rohs_x, rohs_y):
    plt.scatter(x, val, alpha=0.8, color=c)
    plt.scatter(r_x, r_y, color = 'darkorange',marker='_')

plt.scatter(rohs_x[1],dl_y,color='navy',marker='_')

ax.set_yticks(np.arange(-0.1,1.01,0.1))
ax.yaxis.grid(True, color='lightgrey')
ax.set_ylabel('$R^2$',size=50)
ax.set_facecolor('white')

t = pd.read_csv(f'drive/MyDrive/ML/correlations/upbm_dl.txt', sep='\t', header= None)
t.rename({0: "protein"},axis = "columns", inplace = True)
t.index = t['protein']
t = t[[1,2,3,4,5,6,7]]
t = t.dropna()
len(t)

t.columns

# deep learning folks
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(26, 20))

vals, names, xs = [],[],[]
#rohs_x, rohs_y = [], []
for i, col in enumerate(t.columns):
    vals.append(t[col].values)
    names.append(col)
    xs.append(np.random.normal(i + 1, 0.04, t[col].values.shape[0]))  # adds jitter to the data points - can be adjusted
#    rohs_x.append(np.arange(i+1-0.15,i+1+0.17,0.02))

#rohs_y.append([0.67]*len(rohs_x[0]))
#rohs_y.append([0.46]*len(rohs_x[0]))
#rohs_y.append([0.93]*len(rohs_x[0]))

#dl_y = [0.64]*len(rohs_x[0])

#t_1 = (t[['selex']].dropna()).transpose()
#t_2 = (t[['P+A+D']].dropna()).transpose()
#t_3 = (t[['gcPBM']].dropna()).transpose()
#positions = np.arange(3)+1
labeling = ['K-spectrum+shape', 'Dimismatch+shape', 'Deepbind', 'DLBSS', 'CRPT', 'CRPTS', 'DNAffinity']

for k in np.arange(1,8,1):
  t.rename({k: labeling[k-1]},axis = "columns", inplace = True)

bplot = t.boxplot(labels=labeling, ax = ax,showfliers=False, showmeans=True,fontsize = 40, rot=45,medianprops = dict(linestyle='-', linewidth=3, color='black'), 
                   meanprops={"marker":"o","markerfacecolor":"black", "markeredgecolor":"black","markersize":"10"})
#bplot = ax.boxplot(t_1, showfliers=False, showmeans=True,medianprops = dict(linestyle='-', linewidth=3, color='black'), 
#                   meanprops={"marker":"o","markerfacecolor":"black", "markeredgecolor":"black","markersize":"10"})
#bplot = ax.boxplot(t_2, showfliers=False,showmeans=True,medianprops = dict(linestyle='-', linewidth=3, color='black'), 
#                   meanprops={"marker":"o","markerfacecolor":"black", "markeredgecolor":"black","markersize":"10"})
#bplot = ax.boxplot(t_3, showfliers=False,showmeans=True,medianprops = dict(linestyle='-', linewidth=3, color='black'), 
#                   meanprops={"marker":"o","markerfacecolor":"black", "markeredgecolor":"black","markersize":"10"})

palette = ['r', 'g', 'b', 'y', 'darkorange','navy', 'olive']
for x, val, c in zip(xs, vals, palette):
    plt.scatter(x, val, alpha=0.8, color=c)

ax.set_yticks(np.arange(0,1.01,0.1))
ax.yaxis.grid(True, color='lightgrey')
ax.set_ylabel('$R^2$',size=50)
ax.set_facecolor('white')

microarray = t.dropna()

fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(20, 20))

plt.scatter(microarray['DNAffinity'],microarray['K-spectrum+shape'], alpha=0.8, marker='.', s=1000)
#plt.scatter(microarray['DNAffinity'],microarray['Dimismatch+shape'], alpha=0.8, marker='.', s=1000)
#plt.scatter(microarray['DNAffinity'],microarray['Deepbind'], alpha=0.8, marker='.', s=1000)
#plt.scatter(microarray['DNAffinity'],microarray['DLBSS'], alpha=0.8, marker='.', s=1000)
#plt.scatter(microarray['DNAffinity'],microarray['CRPT'], alpha=0.8, marker='.', s=1000)
#plt.scatter(microarray['DNAffinity'],microarray['CRPTS'], alpha=0.8, marker='.', s=1000)

plt.plot(np.arange(0,1.01,1), np.arange(0,1.01,1), color = 'red',linewidth=1)

ax.set_ylabel('K-spectrum+shape',size=50)
#ax.set_ylabel('Dimismatch+shape',size=50)
#ax.set_ylabel('Deepbind',size=50)
#ax.set_ylabel('DLBSS',size=50)
#ax.set_ylabel('CRPT',size=50)
#ax.set_ylabel('CRPTS',size=50)



ax.set_xlabel('DNAffinity',size=50)
ax.set_yticks(np.arange(0,1.01,0.1))
ax.set_xticks(np.arange(0,1.01,0.1))
ax.yaxis.grid(True, color='lightgrey')
ax.xaxis.grid(True, color='lightgrey')
ax.set_facecolor('white')







###### new shape feats
AVG = pd.read_csv(f'drive/MyDrive/ML/input/avg_tetramer_copy.dat', sep=' ')
propeller = pd.read_excel(f'drive/MyDrive/ML/input/propeller_ok.xlsx')
propeller.index = propeller['trimer']
AVG.index = AVG['TETRA']
propeller = propeller['Propeller_ok']

AVG['TWIST_1'] = np.nan
AVG['TWIST_2'] = np.nan

for tetra in list(AVG['TETRA']):
  tri_1 = tetra[:3]
  tri_2 = tetra[1:]
  if tetra[2] == 'g' or tetra[2] == 'J':
    AVG['TWIST_1'][tetra] = np.nan
    AVG['TWIST_2'][tetra] = np.nan
  else:
    AVG['TWIST_1'][tetra] = propeller[tri_1]
    AVG['TWIST_2'][tetra] = propeller[tri_2]

(AVG[['SHIFT','SLIDE','RISE','TILT','ROLL','TWIST','TWIST_1','TWIST_2']]).to_csv(r'drive/MyDrive/ML/input/new_avg.dat', sep=' ', mode='w')

only_shape_corr = pd.read_csv(f'drive/MyDrive/ML/correlations/only_shape_corr.txt', sep=' ')
only_shape_corr.index = only_shape_corr['protein']
only_shape_corr = only_shape_corr['correlation']

new_shape_corr = pd.read_csv(f'drive/MyDrive/ML/correlations/new_avg_corr_2.txt', sep='\t')
new_shape_corr.index = new_shape_corr['protein']
new_shape_corr = new_shape_corr['correlation']

fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(20, 20))

plt.scatter(only_shape_corr, new_shape_corr, alpha=0.8, color='darkorange', marker='.', s=1000)
plt.plot(np.arange(0,1.01,1), np.arange(0,1.01,1), color = 'red',linewidth=1)
ax.set_ylabel('New Propeller Twist',size=30)
ax.set_xlabel('Old model',size=30)
ax.set_yticks(np.arange(0,1.01,0.1))
ax.set_xticks(np.arange(0,1.01,0.1))
ax.yaxis.grid(True, color='lightgrey')
ax.xaxis.grid(True, color='lightgrey')
ax.set_facecolor('white')

fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(15, 10))
plt.plot(new_shape_corr-only_shape_corr)
plt.plot(np.arange(0,19.01,19),[0,0],  color = 'red',linewidth=1)

print(new_shape_corr.mean(), only_shape_corr.mean())
 print((new_shape_corr.drop('tbx15')).mean(), (only_shape_corr.drop('tbx15')).mean())

improv = list(new_shape_corr-only_shape_corr)
print('We improve corr on',sum([1 for k in range(19) if improv[k]>0]))
print('same corr on',sum([1 for k in range(19) if improv[k]==0.0]))
print('We do not improve corr on',sum([1 for k in range(19) if improv[k]<0]))





##pca
AVG = pd.read_csv(f'drive/MyDrive/ML/input/avg_tetramer_copy.dat', sep=' ')
AVG = AVG.drop(np.arange(256,290,1))
AVG.index = AVG['TETRA']
AVG = AVG[['SHIFT','SLIDE','RISE','TILT','ROLL','TWIST']]

tetra_fce = {line.split()[0] : np.array([float(x) for x in line.split()[1:]])  for line in open(f'drive/MyDrive/ML/input/fce_tetramer.dat') if 'SHIFT' not in line}
diag = pd.DataFrame({tt: tetra_fce[tt][list(range(0,36,7))] for tt in tetra_fce.keys()})
diag = diag.transpose()

df = pd.concat([AVG,diag],axis=1)
from sklearn.preprocessing import StandardScaler
scaler = StandardScaler()
shape_matrix = scaler.fit_transform(df)
from sklearn.decomposition import PCA
pca = PCA()
# shape_matrix should have 256 rows (all 4-mers) & 12 columns (all features)
transformed_shape = pca.fit_transform(shape_matrix)  
# transformed_shape will have the same shape as shape_matrix
transformed_shape[:,0] # <- this will be the "main" stiffness vector, combined from many correlated ones, think of it as the 256 stiffness values e.g. for average shift
transformed_shape[:,1] # <- this should be the next one importance-wise etc.

print(pca.explained_variance_ratio_)
print(pca.singular_values_)

df

np.sum(list(abs(pca.components_[0])**2))

pca.components_[1]

pca.components_[2]

pca.components_[3]

sum = 0
for j in range(len(pca.components_)):
    sum = sum + abs(pca.components_[j])*pca.explained_variance_ratio_[j]
    print(j)
sum

shape = ['SHIFT','SLIDE','RISE','TILT','ROLL','TWIST']
average = ['AVG '+feat for feat in shape]
diag = ['DIAG '+feat for feat in shape]
labeling = average+diag

colors = ['red','darkorange','yellow','hotpink','green','blue','cyan','lightgreen','darkgray','purple','lightblue','olive']

fig1, ax1 = plt.subplots()
fig1.set_size_inches(15, 15)
ax1.pie(sum, labels = labeling,startangle=0,  colors=colors,  autopct='%1.1f%%',textprops={'fontsize': 25})
ax1.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
plt.show()

np.sum(sum**2)

(20 + 17) / 146









exclude = []
conf_filter = pd.read_csv(f'drive/MyDrive/ML/conf_filter.txt', sep=',')
rohs_protz = pd.read_csv(f'drive/MyDrive/ML/correlations/electro_corr.txt', sep=' ')
for i in range(len(conf_filter)):
  if conf_filter['confidence'][i] == 0:
    exclude.append(conf_filter['protein'][i])

rohs_protz.index = rohs_protz['protein']
filtered = rohs_protz.drop(exclude)

rohs_protz.rename({'0vs3':'All proteins'},axis = "columns", inplace = True)
filtered.rename({'0vs3':'Filtered'},axis = "columns", inplace = True)
t = pd.DataFrame([rohs_protz['All proteins'],filtered['Filtered']])
t = t.transpose()

fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10, 10))

vals, names, xs = [],[],[]
for i, col in enumerate(t.columns):
    vals.append(t[col].values)
    names.append(col)
    xs.append(np.random.normal(i + 1, 0.04, t[col].values.shape[0]))  # adds jitter to the data points - can be adjusted

#t_1 = (t[['selex']].dropna()).transpose()
#t_2 = (t[['P+A+D']].dropna()).transpose()
#t_3 = (t[['gcPBM']].dropna()).transpose()
#positions = np.arange(3)+1
labeling = ['All proteins', 'Filtered']
bplot = t.boxplot(labels=labeling, ax = ax,showfliers=False, showmeans=True,fontsize = 40, rot=45,medianprops = dict(linestyle='-', linewidth=3, color='black'), 
                   meanprops={"marker":"o","markerfacecolor":"black", "markeredgecolor":"black","markersize":"10"})
#bplot = ax.boxplot(t_1, showfliers=False, showmeans=True,medianprops = dict(linestyle='-', linewidth=3, color='black'), 
#                   meanprops={"marker":"o","markerfacecolor":"black", "markeredgecolor":"black","markersize":"10"})
#bplot = ax.boxplot(t_2, showfliers=False,showmeans=True,medianprops = dict(linestyle='-', linewidth=3, color='black'), 
#                   meanprops={"marker":"o","markerfacecolor":"black", "markeredgecolor":"black","markersize":"10"})
#bplot = ax.boxplot(t_3, showfliers=False,showmeans=True,medianprops = dict(linestyle='-', linewidth=3, color='black'), 
#                   meanprops={"marker":"o","markerfacecolor":"black", "markeredgecolor":"black","markersize":"10"})

palette = ['r', 'g', 'b', 'y']
for x, val, c in zip(xs, vals, palette):
    plt.scatter(x, val, alpha=0.8, color=c)

ax.set_yticks(np.arange(-0.1,1.01,0.1))
ax.yaxis.grid(True, color='lightgrey')
ax.set_ylabel('HT-SELEX $R^2$',size=50)
ax.set_facecolor('white')

t.transpose()

conf_filter

filtered

