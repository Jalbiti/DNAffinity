######################## MACHINE LEARNING ROUTINE ############################################
# from the matrix_aligned file we perform machine learning algorithm                       #
# first we undersample by removing points in given intervals                                 #
#     - it is important to establish the number of divisions L                               #
#     - and the largest value we want to undersample (say, leave the top X points)           #
# secondly we oversample by adding all values larger than the boundary                       #
#     notice that these data may be really similar in affinity and changing a letter or two  #
# finally we run the code choosing features using a dictionary                               #
# the final correlations are stored in
# an array after being averaged                         #
# this is done to avoid the randomness of the random forest                                  #
##############################################################################################

# pick a protein from the following:         
#'Tbx2','Foxc2','Pou1f1', 'Zscan10', 'Sox10', 'Nr5a2', 'Rora', 'Xbp1', 'Sp140', 'Sdccag8', 'Gata4', 'Ar', 'Foxo6', 'Nfil3', 'Zfp202', 'Zfp263', 'Egr3', 'Mlx', 'Tfec', 'Snai1'
#'Ahcftf1', 'Atf3', 'Atf4', 'Dbp', 'Dmrtc2', 'Esrrb', 'Esrrg', 'Klf8', 'Klf9', 'Klf12', 'Mybl2', 'Mzf1', 'Nhlh2', 'Nkx2-9', 'Nr2e1', 'Nr2f1', 'Nr2f6', 'Nr4a2', 'P42pop', 'Pit1', 'Prdm11', 'Rarg', 'Rfx7', 'Rorb', 'Sox3', 'Srebf1', 'Tbx1', 'Tbx4', 'Tbx5', 'Tbx20', 'Zbtb1', 'Zfp300', 'Zfp3', 'Zfp637', 'Zfx', 'Zic5', 'Zkscan1', 'Zkscan5', 'Znf740', 'Foxg1'

#we got it

import sys
import os
from itertools import product
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split, cross_val_score
from sklearn.ensemble import RandomForestRegressor, RandomForestClassifier
from sklearn.preprocessing import MinMaxScaler
from sklearn import metrics
from sklearn.linear_model import LinearRegression
from sklearn.linear_model import Ridge
from sklearn.svm import SVR
from sklearn.model_selection import KFold
import pickle
from collections import defaultdict
import random

np.random.seed(20)

class Dataset:
    def __init__(self, protein, fitxer, model_target='tetramers', randomize_fce=False,
                 chosen_features=['full_fce', 'avg'], score='Median_intensity',
                 selected_tetramers = [0, 1, 2, 3, 4, 5, 6, 7, 8], weighted=False):

        self.protein = protein
        self.pbm35 = fitxer

        if weighted:
            self.pbm35 = self.pbm35[["ID_REF", "VALUE", 'WEIGHT']]
            self.pbm35.columns = ["ID_REF", "VALUE", 'WEIGHT']
        else:
            self.pbm35 = self.pbm35[["ID_REF", "VALUE"]]
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
            tetra_avg = {line.split()[0] : np.array([float(x) for x in line.split()[1:]])  for line in open(f'input/avg_tetramer.dat') if 'SHIFT' not in line}
            try:
                print(name_file)
                sublist = [np.delete(k, int(num_shp_params_avg)) for k in tetra_avg.values()]
            except:
                #print("Not working jeje")
                sublist = tetra_avg.values() 
            tetra_avg = dict(zip(tetra_avg.keys(), sublist))
            self.tetra_avg = np.array([np.concatenate([tetra_avg[otmer[i:i+4]] for i in self.selected_tetramers]) for otmer in self.sequences])
        if 'full_fce' in self.feats or 'diagonal_fce' in self.feats:
            tetra_fce = {line.split()[0] : np.array([float(x) for x in line.split()[1:]])  for line in open(f'input/fce_tetramer.dat') if 'SHIFT' not in line}
            if self.randomize:
                keys = list(tetra_fce.keys())
                permut_keys = np.random.permutation(keys)
                tetra_fce = {key: tetra_fce[val] for key, val in zip(keys, permut_keys)}
            self.tetra_fce = np.array([np.concatenate([tetra_fce[otmer[i:i+4]] for i in self.selected_tetramers]) for otmer in self.sequences])
        # We might want to scramble the matchings to do 'negative control', i.e., remove all physical information from the dataset
        # Let's also keep track of the reduced matrix
        if 'diagonal_fce' in self.feats:
            tetra_fce_reduced = {tt: tetra_fce[tt][list(range(0,36,7))] for tt in tetra_fce.keys()}
            try:
                print(name_file)
                sublist = [np.delete(k, int(num_shp_params_diag)) for k in tetra_fce_reduced.values()]
            except:
                #print("Not working jeje")
                sublist = tetra_fce_reduced.values()
            tetra_fce_reduced = dict(zip(tetra_fce_reduced.keys(), sublist))
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
            self.mgw = {line.split()[0] : [float(x) for x in line.split()[1:]]  for line in open(f'input/mgw_rohs.txt') if 'SHIFT' not in line}
            self.mgw = np.array([np.concatenate([self.mgw[otmer[i:i+4]] for i in self.selected_tetramers]) for otmer in self.sequences])

        if 'electrostatic' in self.feats:
            self.electrostatic = {line.split()[0] : [x for x in line.split()[1:]]  for line in open(f'input/electrostatic.txt') if 'SHIFT' not in line}
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
    regressors = {'random_forest': RandomForestRegressor(n_estimators=10, random_state=None),
                  'linear': LinearRegression(), 'ridge': Ridge(alpha=0.0001), 'svr': SVR()}

    def __init__(self, dataset, len_aln, regressor='random_forest', training_set_size=0.8, weighted=False):
        plt.clf()
        self.data = dataset
        self.X = self.data.features
        self.y = self.data.scores
        if weighted:
            self.weights = self.data.pbm35['WEIGHT']
        else:
            self.weights = None
        self.features = []
        self.l = []
        self.len_aln = len_aln
        self.training_set_size = training_set_size
        self.y_pred, self.r2, self.mse = 3 * [None]
        self.X_train, self.X_test, self.y_train, self.y_test, self.w_train, self.w_test = 6 * [None]

        print("features have shape: ", self.X.shape)
        print("target values have shape: ", self.y.shape)
        self.regressor = Model.regressors[regressor]
        self.train()

    def train(self):
        """
        Splits the data from the Dataset and trains the regressor, also extracting data
        on feature importances
        :return: None
        """
        # splits w/o randomness for training and testing
        if self.weights is not None:
            self.X_train, self.X_test, self.y_train, self.y_test, \
                self.w_train, self.w_test = train_test_split(self.X, self.y.ravel(), self.weights,
                                                             train_size=self.training_set_size,
                                                             random_state=0, shuffle=False)
            self.regressor.fit(self.X_train, self.y_train, self.w_train)
        else:
            self.X_train, self.X_test, self.y_train, \
                self.y_test = train_test_split(self.X, self.y.ravel(), train_size=self.training_set_size,
                                               random_state=0, shuffle=False)
            self.regressor.fit(self.X_train, self.y_train)


    def predict(self, testing_dataset=None):
        """
        Performs the prediction on the test data from the internal (default) or external dataset;
        sets the predictions and corresponding errors/correlations as attributes
        :param testing_dataset: Dataset instance, test data to supply for cross-evaluation (optional)
        :return: None
        """
        if testing_dataset is not None:
            self.X_test, self.y_test = testing_dataset.features, testing_dataset.scores
        self.y_pred = self.regressor.predict(self.X_test).reshape(-1, 1)
        self.y_test = self.y_test.reshape(-1, 1)
        self.r2 = np.round(metrics.r2_score(self.y_test, self.y_pred), 3)
        self.mse = np.round(metrics.mean_absolute_error(self.y_test, self.y_pred), 6)    

    def plot(self):
        """
        Plots the prediction vs ground truth
        :return: None
        """
        plt.rc('xtick', labelsize=15)
        plt.rc('ytick', labelsize=15)
        plt.xlabel('Testing')
        plt.ylabel('Predicted')
        # plt.title(f'MSE = {self.mse}, r$^2$ = {self.r2}')
        plt.xlim = (0, 1)
        plt.ylim = (0, 1)
        plt.scatter(self.y_test, self.y_pred, label=f'r$^2$ = {self.r2}')

    @staticmethod
    def show():
        """
        To allow for multiple overlaid scatter plots, show() can be called separately
        :return: None
        """
        plt.legend()
        plt.show()

