# -*- coding: utf-8 -*-
# !pip install biopython
# !pip install folium
# !curl -O https://raw.githubusercontent.com/jperkel/example_notebook/master/NC_005816.gb

import sys, os
import numpy as np
import time
import pandas as pd
import matplotlib.pyplot as plt
import pickle

from collections import defaultdict
from itertools import product

from regressor_classes import *

if len(sys.argv) > 1:
    protein = sys.argv[1]
    cycle = sys.argv[2]
    data_dir = f'proteins/{protein}'
else:
    protein = 'Gata4'
    cycle = '4'
    data_dir = f'test_data/{protein}'

#####################################
# files required to run:
# '{data_dir}/{cycle}.txt'
#####################################

deca = pd.read_csv(f'{data_dir}/{cycle}.txt', sep='\t')

strings = list(deca["Kmer"])
values = list(deca["Affinity"])
# counts = list(deca["ExpectedCount"])
counts = list(deca["ObservedCount"])
prob = list(deca["Probability"])

# undersampling

prob_bound = sorted(prob)[int(len(prob)*0.9)]

mat = []
for i in np.arange(0, len(prob)):
    if prob[i] > prob_bound:
        mat.append([strings[i], values[i]])

# ordered

if not os.path.isdir('output_selex'):
    os.mkdir('output_selex')
if not os.path.isdir(f'output_selex/{protein}'):
    os.mkdir(f'output_selex/{protein}')


with open(f'output_selex/{protein}/{protein}_training_ordered.txt', 'w') as file:
    file.write('ID_REF\tVALUE\n')
    for vector in mat:
        file.write("%s\t" % vector[0])
        file.write("%s\n" % vector[1])

# writes the randomized file ready for training

np.random.shuffle(mat)

with open(f'output_selex/{protein}/SELEX_training.txt', 'w') as file:
    file.write('ID_REF\tVALUE\n')
    for vector in mat:
        file.write("%s\t" % vector[0])
        file.write("%s\n" % vector[1])


the_features = {0: ['diagonal_fce'], 1: ['presence_tetramer'], 2: ['avg'],
                3: ['presence_tetramer', 'avg', 'diagonal_fce'],
                4: ['presence_tetramer', 'avg', 'diagonal_fce', 'electrostatic'],
                5: ['avg', 'diagonal_fce']}
len_aln = len(strings[0])

df_train = pd.read_csv(f'output_selex/{protein}/SELEX_training.txt', sep='\t')


chosen_features = the_features[3]
model_target = 'octamers'                                 
regressor = 'random_forest'                                
training_set_size = 0.9
randomize_fce = False                                     
score = 'Median_intensity'
selected_tetramers = list(np.arange(0, len_aln-3))

# time
start_time = time.time()
dset_train = Dataset(protein, df_train, model_target, randomize_fce, chosen_features, score, selected_tetramers)
model = Model(dset_train, len_aln, regressor, training_set_size)
print('the SELEX model has been trained')
print("It took %s seconds" % (time.time() - start_time))

model.predict()
print(model.y_test.shape, model.y_pred.shape)
print('The r2 is ', model.r2)

# save the model to disk
filename = f'trained_models/{protein}_finalized_model.sav'
pickle.dump(model, open(filename, 'wb'), protocol=4)

