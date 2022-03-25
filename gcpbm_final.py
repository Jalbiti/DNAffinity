# -*- coding: utf-8 -*-
### conda code v1

!pip install biopython
!pip install folium
!curl -O https://raw.githubusercontent.com/jperkel/example_notebook/master/NC_005816.gb

##################

import time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from dataset import Dataset
from model import Model

protein = 'cbf1'
concentration = '100' # '100' or '200'
data_dir = f'drive/MyDrive/ML/gcPBM/{protein}'

#####################################
# files required to run:
# '{data_dir}/{protein}_{concentration}.txt'
# '{data_dir}/cbf1_freq_matrix_6.txt'
#####################################

# yeast

raw_data = pd.read_csv(f'{data_dir}/{protein}_{concentration}.txt', sep='\t')
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
print(strings[0][12:18])



freq_matrix = pd.read_csv(f'{data_dir}/cbf1_freq_matrix_6.txt', delim_whitespace=True)
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

discarded = list(discarded)
print('We have discarded',len(discarded)/len(strings)*100,'% of the data')


## options:

# no smart under, no weighting
df_train = pd.DataFrame({'ID_REF': strings, 'VALUE': values})
df_train = df_train.drop(discarded)
df_train = df_train.reset_index(drop=True)
# df_train = df_train.sample(frac=1).reset_index(drop=True)[0:15000]
len(df_train)
# # smart under, no weighting
# df_train = pd.read_csv(f'{data_dir}/{protein}_training.txt', sep='\t')

# # no smart under, weighting
# df_train = pd.DataFrame({'ID_REF': strings, 'VALUE': values, 'WEIGHT':weights})
# df_train = df_train.sample(frac=1).reset_index(drop=True)[0:15000]



the_features = {0: ['electrostatic'], 1: ['presence_tetramer'], 2: ['avg'], 3: ['diagonal_fce'],
                4: ['presence_tetramer', 'avg', 'diagonal_fce','electrostatic'], 5: ['avg', 'diagonal_fce']}
len_aln = len(strings[0])

chosen_features = the_features[4]
model_target = 'octamers'
regressor = 'random_forest'
training_set_size = 0.9
randomize_fce = False
score = 'Median_intensity'
selected_tetramers = list(np.arange(0,len_aln-3))

## time

start_time = time.time()
dset_train = Dataset(protein, df_train, model_target, randomize_fce, chosen_features, score, selected_tetramers)
model = Model(dset_train, len_aln, regressor, training_set_size)
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

df = pd.DataFrame(model.features)

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

print(df[0:10])
