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
# 'Tbx2','Foxc2','Pou1f1', 'Zscan10', 'Sox10', 'Nr5a2', 'Rora', 'Xbp1', 'Sp140', 'Sdccag8', 'Gata4', 'Ar', 'Foxo6',
# 'Nfil3', 'Zfp202', 'Zfp263', 'Egr3', 'Mlx', 'Tfec', 'Snai1', 'Ahcftf1', 'Atf3', 'Atf4', 'Dbp', 'Dmrtc2', 'Esrrb',
# 'Esrrg', 'Klf8', 'Klf9', 'Klf12', 'Mybl2', 'Mzf1', 'Nhlh2', 'Nkx2-9', 'Nr2e1', 'Nr2f1', 'Nr2f6', 'Nr4a2', 'P42pop',
# 'Pit1', 'Prdm11', 'Rarg', 'Rfx7', 'Rorb', 'Sox3', 'Srebf1', 'Tbx1', 'Tbx4', 'Tbx5', 'Tbx20', 'Zbtb1', 'Zfp300',
# 'Zfp3', 'Zfp637', 'Zfx', 'Zic5', 'Zkscan1', 'Zkscan5', 'Znf740', 'Foxg1'
import os.path
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import time
import pickle

from collections import defaultdict
import random

from regressor_classes import *

if len(sys.argv) > 1:
    protein = sys.argv[1]
    data_dir = f'proteins/{protein}'
else:
    protein = 'Gata4'
    data_dir = f'test_data/{protein}'

#####################################
# files required to run:
# '{data_dir}/{protein}_training.txt'
#####################################
np.random.seed(20)
protein = sys.argv[1]

the_features = {0: ['diagonal_fce'], 1: ['presence_tetramer'], 2: ['avg'],
                3: ['presence_tetramer', 'avg', 'diagonal_fce'],
                4: ['presence_tetramer', 'avg', 'diagonal_fce', 'electrostatic'],
                5: ['avg', 'diagonal_fce']}

df_train = pd.read_csv(f'{data_dir}/{protein}_training.txt', sep='\t')
strings = list(df_train["ID_REF"])
len_aln = len(strings[0])

chosen_features = the_features[3]
model_target = 'octamers'
regressor = 'random_forest'
training_set_size = 0.9
randomize_fce = False
score = 'Median_intensity'
selected_tetramers = list(np.arange(0,len_aln-3))

# timings and the actual training
start_time = time.time()
dset_train = Dataset(protein, df_train, model_target, randomize_fce, chosen_features, score, selected_tetramers, weighted=True)
model = Model(dset_train, len_aln, regressor, training_set_size, weighted=True)
print('the uPBM model has been trained')
print("It took %s seconds" % (time.time() - start_time))

model.predict()
print(model.y_test.shape, model.y_pred.shape)
print('The r2 is ', model.r2)

# save the model to disk
filename = f'trained_models/{protein}_finalized_model.sav'
pickle.dump(model, open(filename, 'wb'))

