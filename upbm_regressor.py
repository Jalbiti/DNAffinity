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

from collections import defaultdict
import random

from dataset import Dataset
from model import Model

if len(sys.argv) > 1:
    protein = sys.argv[1]
    data_dir = f'proteins/{protein}'
else:
    protein = 'Gata4'
    data_dir = f'test_data'

#####################################
# files required to run:
# '{data_dir}/{protein}_matrix_aligned.txt'
#####################################

# first read the processed data -- the 12mers and their corresponding affinity
# the _over file has some extra 12mers per sequence, namely those with affinity
# larger than max_aff/2

deca = pd.read_csv(f'{data_dir}/{protein}_matrix_aligned.txt', sep='\t')

strings = list(deca["ID_REF"])
values = list(deca["VALUE"])
weights = list(deca["WEIGHTEDSCORE"])  # pick either WEIGHT or WEIGHTEDSCORE

values = [(value-min(values))/(max(values)-min(values)) for value in values]

dictionary = dict(zip(values, strings))
weight_dict = dict(zip(values, weights))

len_aln = len(strings[0])

print(dictionary[1])  # quick check that everything is alright

# undersampling

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
start = 1  # i.e. starts picking values larger than start/L from the normalized affinites

for i in np.arange(start, L+1):
    if len(interval_dict[i]) > 0:
        if i <= L*upper_bound:
            random_values = [random.choice(interval_dict[i]) for j in np.arange(0, 20)]
            random_values = set(random_values)
            value = random.choice(interval_dict[i])
            training_ordered.append([dictionary[value], value, weight_dict[value]])
        if i > L*upper_bound:
            for value in interval_dict[i]:
                training_ordered.append([dictionary[value], value, weight_dict[value]])

# writes the file w/o randomization (not used but in case we need it)
with open(f'{data_dir}/{protein}_training_ordered.txt', 'w') as file:
    file.write('ID_REF\tVALUE\tWEIGHT\n')
    for vector in training_ordered:
        file.write("%s\t" % vector[0])
        file.write("%s\t" % vector[1])
        file.write("%s\n" % vector[2])

# writes the randomized file ready for training
np.random.shuffle(training_ordered)

with open(f'{data_dir}/{protein}_training.txt', 'w') as file:
    file.write('ID_REF\tVALUE\tWEIGHT\n')
    for vector in training_ordered:
        file.write("%s\t" % vector[0])
        file.write("%s\t" % vector[1])
        file.write("%s\n" % vector[2])


the_features = {0: ['diagonal_fce'], 1: ['presence_tetramer'], 2: ['avg'],
                3: ['presence_tetramer', 'avg', 'diagonal_fce'],
                4: ['presence_tetramer', 'avg', 'diagonal_fce', 'electrostatic'],
                5: ['avg', 'diagonal_fce']}
len_aln = len(strings[0])

df_train = pd.read_csv(f'{data_dir}/{protein}_training.txt', sep='\t')

chosen_features = the_features[4]
model_target = 'octamers'
regressor = 'random_forest'
training_set_size = 0.9
randomize_fce = False
score = 'Median_intensity'
selected_tetramers = list(np.arange(0, len_aln-3))

# timings and the actual training
start_time = time.time()
dset_train = Dataset(protein, df_train, model_target, randomize_fce, chosen_features, score, selected_tetramers,
                     weighted=True)
model = Model(dset_train, len_aln, regressor, training_set_size, weighted=True)
print('the uPBM model has been trained')
print("It took %s seconds" % (time.time() - start_time))

model.predict()
print(model.y_test.shape, model.y_pred.shape)
print('The correlation is ', model.r2)

if not os.path.isdir('output_upbm'):
    os.mkdir('output_upbm')
if not os.path.isdir(f'output_upbm/{protein}'):
    os.mkdir(f'output_upbm/{protein}')

plt.savefig(f'output_upbm/{protein}/{protein}_corr.png')

with open(f'output_upbm/uupbm_correlations.txt', 'a') as file:
    file.write("%s\t" % protein)
    file.write("%s\n" % model.r2)

# features

df = pd.DataFrame(model.features)

index = [k for k in range(model.X.shape[1]) if df[0][k] == 'Presence']
p = sum([df[1][k] for k in index])

index = [k for k in range(model.X.shape[1]) if df[0][k] == 'Electro']
e = sum([df[1][k] for k in index])

shape = 1-e-p

y = np.array([p, e, shape])

fig = plt.gcf()
fig.set_size_inches(8, 8)
labels = ["Presence", "Electrostatic", "Shape"]
patches, texts = plt.pie(y, startangle=0)
plt.legend(patches, labels, loc="best")
plt.savefig(f'output_upbm/{protein}/features.png')

with open(f'output_upbm/upbm_features.txt', 'a') as file:
    file.write("%s\t" % protein)
    file.write("%s\t" % p)
    file.write("%s\t" % e)
    file.write("%s\n" % shape)

print(df[0:30])
