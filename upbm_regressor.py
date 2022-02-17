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


import sys
import os
from itertools import product
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestRegressor, RandomForestClassifier
from sklearn.preprocessing import MinMaxScaler
from sklearn import metrics
from sklearn.linear_model import LinearRegression
from sklearn.linear_model import Ridge
from sklearn.svm import SVR

from collections import defaultdict
import random


protein = sys.argv[1]


# first read the processed data -- the 12mers and their corresponding affinity
# the _over file has some extra 12mers per sequence, namely those with affinity
# larger than max_aff/2

deca = pd.read_csv(f'proteins/{protein}/{protein}_matrix_aligned.txt', sep='\t')

strings = list(deca["ID_REF"])
values = list(deca["VALUE"])
weights = list(deca["WEIGHTEDSCORE"])  # pick either WEIGHT or WEIGHTEDSCORE

values = [(value-min(values))/(max(values)-min(values)) for value in values]

dictionary = dict(zip(values, strings))
weight_dict = dict(zip(values, weights))

len_aln = len(strings[0])

dictionary[1] # quick check that everything is alright

#### undersampling

# number of divisions we are using
L = 25000

# affinitiy boundaries
sort_aff = sorted(values)
# pick top points ~ proportion 1:500 before undersampling i.e. we pick 50 top points not to undersample
upper_bound = sort_aff[-1000]

# interval dictionary --
interval_dict = defaultdict(list)
for value in values:
    interval_dict[int(value*L)].append(value)

# creates a matrix of values: we pick one value for the intervals where i<L/2
# and up to two for i>L/2

training_ordered = []
start = 1 # i.e. starts picking values larger than start/L from the normalized affinites

for i in np.arange(start,L+1):
    if len(interval_dict[i])>0:
        if i <= L*upper_bound:
            random_values = [random.choice(interval_dict[i]) for j in np.arange(0,20)]
            random_values = set(random_values)
            value = random.choice(interval_dict[i])
            training_ordered.append([dictionary[value], value, weight_dict[value]])
        if i > L*upper_bound:
            for value in interval_dict[i]:
                training_ordered.append([dictionary[value], value, weight_dict[value]])

# writes the file w/o randomization (not used but in case we need it)

with open(f'proteins/{protein}/{protein}_training_ordered.txt','w') as file:
    file.write('ID_REF\tVALUE\tWEIGHT\n')

with open(f'proteins/{protein}/{protein}_training_ordered.txt','a') as file:
    for vector in training_ordered:
        file.write("%s\t" % vector[0])
        file.write("%s\t" % vector[1])
        file.write("%s\n" % vector[2])


# writes the randomized file ready for training


np.random.shuffle(training_ordered)

with open(f'proteins/{protein}/{protein}_training.txt','w') as file:
    file.write('ID_REF\tVALUE\tWEIGHT\n')
    
with open(f'proteins/{protein}/{protein}_training.txt','a') as file:
    for vector in training_ordered:
        file.write("%s\t" % vector[0])
        file.write("%s\t" % vector[1])
        file.write("%s\n" % vector[2])


# classes

l=[]
features=[]

class Dataset:
    def __init__(self, protein, fitxer, model_target='tetramers', randomize_fce=False,
                 chosen_features=['full_fce', 'avg'], score='Median_intensity',
                 selected_tetramers = [0, 1, 2, 3, 4, 5, 6, 7, 8], pbm = 10):

        self.protein = protein
        self.pbm35 = fitxer

        ##self.pbm35 = pd.read_csv(f'drive/My Drive/ML/{protein}/{fitxer}.txt', sep='\t')
        ##self.pbm35 = pd.read_csv(f'fitxer, sep='\t')

        self.pbm35 = self.pbm35[["ID_REF","VALUE", 'WEIGHT']]
        self.pbm35.columns = ["ID_REF", "VALUE", 'WEIGHT']
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
            tetra_avg = {line.split()[0] : np.array([float(x) for x in line.split()[1:]])  for line in open(f'/orozco/projects/proteinBinding/PBM_SELEX/input/avg_tetramer.dat') if 'SHIFT' not in line}
            self.tetra_avg = np.array([np.concatenate([tetra_avg[otmer[i:i+4]] for i in self.selected_tetramers]) for otmer in self.sequences])
        if 'full_fce' in self.feats or 'diagonal_fce' in self.feats:
            tetra_fce = {line.split()[0] : np.array([float(x) for x in line.split()[1:]])  for line in open(f'/orozco/projects/proteinBinding/PBM_SELEX/input/fce_tetramer.dat') if 'SHIFT' not in line}
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
            self.mgw = {line.split()[0] : [float(x) for x in line.split()[1:]]  for line in open(f'/orozco/projects/proteinBinding/PBM_SELEX/input/mgw_rohs.txt') if 'SHIFT' not in line}
            self.mgw = np.array([np.concatenate([self.mgw[otmer[i:i+4]] for i in self.selected_tetramers]) for otmer in self.sequences])

        if 'electrostatic' in self.feats:
            self.electrostatic = {line.split()[0] : [x for x in line.split()[1:]]  for line in open(f'/orozco/projects/proteinBinding/PBM_SELEX/input/electrostatic.txt') if 'SHIFT' not in line}
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
        self.weights = self.data.pbm35['WEIGHT']

        print("features have shape: ", self.X.shape)
        print("target values have shape: ", self.y.shape)
        self.regressor = Model.regressors[regressor]
        # splits w/o randomness for training and testing
        self.X_train, self.X_test, self.y_train, self.y_test, self.w_train, self.w_test = train_test_split(self.X, self.y, self.weights, train_size=training_set_size, random_state=0, shuffle=False)
        self.regressor.fit(self.X_train, self.y_train, self.w_train)

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
    
    def predict(self, testing_dataset=None):
        if testing_dataset is not None:
            self.X_test, self.y_test = testing_dataset.features, testing_dataset.scores
        self.y_pred = self.regressor.predict(self.X_test).reshape(-1,1)
        self.y_test = self.y_test.reshape(-1,1)
        self.r2 = np.round(metrics.r2_score(self.y_test, self.y_pred, self.w_test), 3) # with weights
        #self.r2 = np.round(metrics.r2_score(self.y_test, self.y_pred), 3)           # w/o weights
        self.mse = np.round(metrics.mean_absolute_error(self.y_test, self.y_pred), 6)

    def plot(self):
        plt.rc('xtick', labelsize=15)
        plt.rc('ytick', labelsize=15)
        plt.xlabel('Testing')
        plt.ylabel('Predicted')
        #plt.title(f'MSE = {self.mse}, r$^2$ = {self.r2}')
        plt.xlim = (0,1)
        plt.ylim = (0,1)
        plt.scatter(self.y_test, self.y_pred, label=f'r$^2$ = {self.r2}')

    def show(self):
        plt.legend()
        plt.show()


the_features = defaultdict(list)
the_features = {0: ['diagonal_fce'], 1: ['presence_tetramer'], 2: ['avg'], 3: ['presence_tetramer', 'avg', 'diagonal_fce'], 4: ['presence_tetramer', 'avg', 'diagonal_fce','electrostatic'], 5: ['avg', 'diagonal_fce']}
len_aln = len(strings[0])

df_train = pd.read_csv(f'proteins/{protein}/{protein}_training.txt', sep='\t')

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

print('the uPBM model has been trained')
print("It took %s seconds" % (time.time() - start_time))

model.predict()
print(model.y_test.shape, model.y_pred.shape)
print('The correlation is ', model.r2)

plt.savefig(f'output_upbm/{protein}/{protein}_corr.png')

with open(f'output_upbm/uupbm_correlations.txt','a') as file:
   file.write("%s\t" % protein)
   file.write("%s\n" % model.r2)

# features

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
plt.savefig(f'output_upbm/{protein}/features.png')

with open(f'output_upbm/upbm_features.txt','a') as file:
   file.write("%s\t" % protein)
   file.write("%s\t" % p)
   file.write("%s\t" % e)
   file.write("%s\n" % shape)

print(df[0:30])


