import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from Bio import SeqIO
from Bio import pairwise2

from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

protein = 'cbf1'
concentration = '100'  # '100' or '200'
experiment = 'pb'  # 'chip' or 'pb' for in vivo or in vitro resp.

# genomic

raw_peaks = pd.read_csv(f'drive/MyDrive/ML/gcPBM/{protein}/{protein}_{experiment}-exo.gff', sep='\t', header=None)
proc_peaks = raw_peaks[[0, 3, 4, 5, 8]]
proc_peaks = proc_peaks.rename(columns={0: "chr", 3: "start", 4: "end", 5: "score", 8: "distance"})
proc_peaks["distance"] = proc_peaks["distance"].str[12:].astype(int)
proc_peaks["start"] = proc_peaks["start"] - proc_peaks["distance"]
for i in range(len(proc_peaks)):
    if proc_peaks["start"][i] < 0:
        proc_peaks["start"][i] = 0
proc_peaks["end"] = proc_peaks["end"] + proc_peaks["distance"]
proc_peaks = proc_peaks[["chr", "start", "end", "score"]]
min_val = min(proc_peaks["score"])
normalization = max(proc_peaks["score"]) - min(proc_peaks["score"])
proc_peaks["score"] = list(((proc_peaks["score"] - min_val) / normalization))

with open(f'drive/MyDrive/ML/gcPBM/{protein}/{protein}_{experiment}-exo.bedgraph', 'w') as file:
    for i in range(len(proc_peaks)):
        file.write("%s\t" % proc_peaks.iloc[i][0])
        file.write("%s\t" % proc_peaks.iloc[i][1])
        file.write("%s\t" % proc_peaks.iloc[i][2])
        file.write("%s\n" % proc_peaks.iloc[i][3])

for k in range(1, 17):
    chrom = 'chr' + str(k)
    print(chrom)
    chrsome = proc_peaks[proc_peaks['chr'] == chrom]
    chrsome = chrsome.sort_values("start")
    chrsome = chrsome.reset_index(drop=True)
    print(len(chrsome) / (chrsome['end'][len(chrsome) - 1] - chrsome['start'][0]) * 100000)

# intervals: 25-26 (1), 27-28 (2), 29-30 (3), 32-33 (4), 34-35 (5)

chrsome = 'chr14'
chr_r = 'chrXIV'
number = 2
start_bp = 100000 * (number - 1) + 1
end_bp = 100000 * number
# end_bp = 230218 # chr1
# end_bp = 813184 # chr2
# end_bp = 316620 # chr3
# end_bp = 1531933 # chr4
# end_bp = 576874 # chr5
# end_bp = 270161 # chr6
# end_bp = 1090940 # chr7
# end_bp = 562643 # chr8
# end_bp = 439888 # chr9
# end_bp = 745752 # chr10
# end_bp = 666816 # chr11
# end_bp = 1078177 # chr12
# end_bp = 924431 # chr13
# end_bp = 784333 # chr14
# end_bp = 1091291 # chr15
# end_bp =  948066 # chr16

print(start_bp, end_bp, end_bp - start_bp)

fasta = SeqIO.read(f"drive/MyDrive/ML/gcPBM/{chrsome}_{number}.fa", "fasta")
sequence = str(fasta.seq)

experiment = 'chip'

# genomic

raw_peaks = pd.read_csv(f'drive/MyDrive/ML/gcPBM/{protein}/{protein}_{experiment}-exo.gff', sep='\t', header=None)
proc_peaks = raw_peaks[[0, 3, 4, 5, 8]]
proc_peaks = proc_peaks.rename(columns={0: "chr", 3: "start", 4: "end", 5: "score", 8: "distance"})
proc_peaks["distance"] = proc_peaks["distance"].str[12:].astype(int)
proc_peaks["start"] = proc_peaks["start"] - proc_peaks["distance"]
for i in range(len(proc_peaks)):
    if proc_peaks["start"][i] < 0:
        proc_peaks["start"][i] = 0
proc_peaks["end"] = proc_peaks["end"] + proc_peaks["distance"]
proc_peaks = proc_peaks[["chr", "start", "end", "score"]]
min_val = min(proc_peaks["score"])
normalization = max(proc_peaks["score"]) - min(proc_peaks["score"])
proc_peaks["score"] = list(((proc_peaks["score"] - min_val) / normalization))

with open(f'drive/MyDrive/ML/gcPBM/{protein}/{protein}_{experiment}-exo.bedgraph', 'w') as file:
    for i in range(len(proc_peaks)):
        file.write("%s\t" % proc_peaks.iloc[i][0])
        file.write("%s\t" % proc_peaks.iloc[i][1])
        file.write("%s\t" % proc_peaks.iloc[i][2])
        file.write("%s\n" % proc_peaks.iloc[i][3])

chr4 = proc_peaks[proc_peaks['chr'] == chrsome]
chr4 = chr4.sort_values("start")  # [1:12], [13:32]
chr4 = chr4.reset_index(drop=True)
print(chr4)

experiment = 'pb'

## genomic

raw_peaks = pd.read_csv(f'drive/MyDrive/ML/gcPBM/{protein}/{protein}_{experiment}-exo.gff', sep='\t', header=None)
proc_peaks = raw_peaks[[0, 3, 4, 5, 8]]
proc_peaks = proc_peaks.rename(columns={0: "chr", 3: "start", 4: "end", 5: "score", 8: "distance"})
proc_peaks["distance"] = proc_peaks["distance"].str[12:].astype(int)
proc_peaks["start"] = proc_peaks["start"] - proc_peaks["distance"]
for i in range(len(proc_peaks)):
    if proc_peaks["start"][i] < 0:
        proc_peaks["start"][i] = 0
proc_peaks["end"] = proc_peaks["end"] + proc_peaks["distance"]
proc_peaks = proc_peaks[["chr", "start", "end", "score"]]
min_val = min(proc_peaks["score"])
normalization = max(proc_peaks["score"]) - min(proc_peaks["score"])
proc_peaks["score"] = list(((proc_peaks["score"] - min_val) / normalization))

with open(f'drive/MyDrive/ML/gcPBM/{protein}/{protein}_{experiment}-exo.bedgraph', 'w') as file:
    for i in range(len(proc_peaks)):
        file.write("%s\t" % proc_peaks.iloc[i][0])
        file.write("%s\t" % proc_peaks.iloc[i][1])
        file.write("%s\t" % proc_peaks.iloc[i][2])
        file.write("%s\n" % proc_peaks.iloc[i][3])

chr4 = proc_peaks[proc_peaks['chr'] == chrsome]
chr4 = chr4.sort_values("start")[34:49]  # [12:29], [30:42]
chr4 = chr4.reset_index(drop=True)
print(chr4)

new_peaks = pd.read_csv(f'drive/MyDrive/ML/gcPBM/{protein}/{protein}_rep1_peaks.bed', sep='\t', header=None)
new_peaks = new_peaks[[0, 1, 2, 4]]
new_peaks = new_peaks.rename(columns={0: "chr", 1: "start", 2: "end", 4: "score"})
min_val = min(new_peaks["score"])
normalization = max(new_peaks["score"]) - min(new_peaks["score"])
new_peaks["score"] = list(((new_peaks["score"] - min_val) / normalization))

with open(f'drive/MyDrive/ML/gcPBM/{protein}/{protein}_rep1_peaks.bedgraph', 'w') as file:
    for i in range(len(new_peaks)):
        file.write("%s\t" % new_peaks.iloc[i][0])
        file.write("%s\t" % new_peaks.iloc[i][1])
        file.write("%s\t" % new_peaks.iloc[i][2])
        file.write("%s\n" % new_peaks.iloc[i][3])

chr4 = new_peaks[new_peaks['chr'] == chr_r]
chr4 = chr4.sort_values("start")  # [12:29], [30:42]
chr4 = chr4.reset_index(drop=True)
print(chr4)

new_peaks = pd.read_csv(f'drive/MyDrive/ML/gcPBM/{protein}/{protein}_rep2_peaks.bed', sep='\t', header=None)
new_peaks = new_peaks[[0, 1, 2, 4]]
new_peaks = new_peaks.rename(columns={0: "chr", 1: "start", 2: "end", 4: "score"})
min_val = min(new_peaks["score"])
normalization = max(new_peaks["score"]) - min(new_peaks["score"])
new_peaks["score"] = list(((new_peaks["score"] - min_val) / normalization))

with open(f'drive/MyDrive/ML/gcPBM/{protein}/{protein}_rep2_peaks.bedgraph', 'w') as file:
    for i in range(len(new_peaks)):
        file.write("%s\t" % new_peaks.iloc[i][0])
        file.write("%s\t" % new_peaks.iloc[i][1])
        file.write("%s\t" % new_peaks.iloc[i][2])
        file.write("%s\n" % new_peaks.iloc[i][3])

chr4 = new_peaks[new_peaks['chr'] == chr_r]
chr4 = chr4.sort_values("start")  # [12:29], [30:42]
chr4 = chr4.reset_index(drop=True)
print(chr4)

data = pd.read_csv(f'drive/MyDrive/ML/gcPBM/{protein}/{protein}_dictionary_{chrsome}_{number}.txt', sep='\t',
                   header=None)

genomic_tmer_dict = dict(data.values.tolist())

# make sure you're using the correct dict
print(list(genomic_tmer_dict.keys())[0])
print(sequence[start_bp - start_bp:start_bp - start_bp + 30])

# once we have the dictionary, time to plot stuff

# our affinities:
profile_matrix = []
affinities = []

for i in np.arange(start_bp - start_bp, end_bp - 29 - start_bp):
    start = i
    end = start + 30
    seq = sequence[start:end]
    # if i>start_bp-250000+49720-25 and i<start_bp-250000+49720+5:
    #   print(genomic_tmer_dict[seq], seq)
    profile_matrix.append([seq, genomic_tmer_dict[seq], start + start_bp, end + start_bp])
    aff = 0
    if genomic_tmer_dict[seq] > 0.3:
        aff = genomic_tmer_dict[seq]
    affinities.append(aff)

with open(f'drive/MyDrive/ML/gcPBM/{protein}/affinities_{protein}.bedgraph', 'w') as file:
    for prof in profile_matrix:
        file.write(f"{chrsome}\t")
        file.write("%s\t" % prof[2])
        file.write("%s\t" % prof[3])
        file.write("%s\n" % prof[1])

experiment = 'chip'

# Read in Chip-seq peaks from protein
chip_peaks = pd.read_csv(f'drive/MyDrive/ML/gcPBM/{protein}/{protein}_{experiment}-exo.bedgraph', sep='\t', header=None)
chip_peaks.rename({0: "chr", 1: "start", 2: "end", 3: "score"}, axis="columns", inplace=True)

peaks_chr_this = chip_peaks[chip_peaks['chr'] == chrsome]
filtered_peaks = peaks_chr_this[(peaks_chr_this["start"] > start_bp) & (peaks_chr_this["end"] <= end_bp)]
peaks = filtered_peaks.copy()
peaks = peaks.sort_values("start")
peaks = peaks.reset_index(drop=True)

# Set indices starting at our 0
peaks['start'] = peaks['start'] - start_bp
peaks['end'] = peaks['end'] - start_bp
peaks['length'] = peaks['end'] - peaks['start']
peaks = peaks[['start', 'end', 'length', 'score']]

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
            chip_step_chip[start:end] = [peaks["score"].values[i]] * (end - start)  # will input 1 for length of peak
            i += 1
    else:
        break

print(peaks)

experiment = 'pb'

# Read in Chip-seq peaks from protein
chip_peaks = pd.read_csv(f'drive/MyDrive/ML/gcPBM/{protein}/{protein}_{experiment}-exo.bedgraph', sep='\t', header=None)
chip_peaks.rename({0: "chr", 1: "start", 2: "end", 3: "score"}, axis="columns", inplace=True)

peaks_chr_this = chip_peaks[chip_peaks['chr'] == chrsome]
filtered_peaks = peaks_chr_this[(peaks_chr_this["start"] > start_bp) & (peaks_chr_this["end"] <= end_bp)]
peaks = filtered_peaks.copy()
peaks = peaks.sort_values("start")
peaks = peaks.reset_index(drop=True)

# Set indices starting at our 0
peaks['start'] = peaks['start'] - start_bp
peaks['end'] = peaks['end'] - start_bp
peaks['length'] = peaks['end'] - peaks['start']
peaks = peaks[['start', 'end', 'length', 'score']]

# create step plot from chip seq data
chip_step_pb = [0] * len(affinities)
i = 0
for j in range(len(affinities)):
    if i < len(peaks):
        if j == peaks["start"].values[i]:
            start = peaks["start"].values[i]
            if peaks["end"].values[i] >= len(affinities):
                end = len(affinities)
            else:
                end = peaks["end"].values[i]

            # chip_step[start:end] = [np.max(affinities)] * (end - start) # will input max for length of peak
            chip_step_pb[start:end] = [peaks["score"].values[i]] * (end - start)  # will input 1 for length of peak
            i += 1
    else:
        break

# new peaks

# Read in Chip-seq peaks from protein
chip_peaks = pd.read_csv(f'drive/MyDrive/ML/gcPBM/{protein}/{protein}_rep1_peaks.bedgraph', sep='\t', header=None)
chip_peaks.rename({0: "chr", 1: "start", 2: "end", 3: "score"}, axis="columns", inplace=True)

peaks_chr_this = chip_peaks[chip_peaks['chr'] == chr_r]
filtered_peaks = peaks_chr_this[(peaks_chr_this["start"] > start_bp) & (peaks_chr_this["end"] <= end_bp)]
peaks = filtered_peaks.copy()
peaks = peaks.sort_values("start")
peaks = peaks.reset_index(drop=True)

# Set indices starting at our 0
peaks['start'] = peaks['start'] - start_bp
peaks['end'] = peaks['end'] - start_bp
peaks['length'] = peaks['end'] - peaks['start']
peaks = peaks[['start', 'end', 'length', 'score']]

# create step plot from chip seq data
chip_step_new = [0] * len(affinities)
i = 0
for j in range(len(affinities)):
    if i < len(peaks):
        if j == peaks["start"].values[i]:
            start = peaks["start"].values[i]
            if peaks["end"].values[i] >= len(affinities):
                end = len(affinities)
            else:
                end = peaks["end"].values[i]

            # chip_step[start:end] = [np.max(affinities)] * (end - start) # will input max for length of peak
            chip_step_new[start:end] = [peaks["score"].values[i]] * (end - start)  # will input 1 for length of peak
            i += 1
    else:
        break

# new peaks second replicate

# Read in Chip-seq peaks from protein
chip_peaks = pd.read_csv(f'drive/MyDrive/ML/gcPBM/{protein}/{protein}_rep2_peaks.bedgraph', sep='\t', header=None)
chip_peaks.rename({0: "chr", 1: "start", 2: "end", 3: "score"}, axis="columns", inplace=True)

peaks_chr_this = chip_peaks[chip_peaks['chr'] == chr_r]
filtered_peaks = peaks_chr_this[(peaks_chr_this["start"] > start_bp) & (peaks_chr_this["end"] <= end_bp)]
peaks = filtered_peaks.copy()
peaks = peaks.sort_values("start")
peaks = peaks.reset_index(drop=True)

# Set indices starting at our 0
peaks['start'] = peaks['start'] - start_bp
peaks['end'] = peaks['end'] - start_bp
peaks['length'] = peaks['end'] - peaks['start']
peaks = peaks[['start', 'end', 'length', 'score']]

# create step plot from chip seq data
chip_step_new_2 = [0] * len(affinities)
i = 0
for j in range(len(affinities)):
    if i < len(peaks):
        if j == peaks["start"].values[i]:
            start = peaks["start"].values[i]
            if peaks["end"].values[i] >= len(affinities):
                end = len(affinities)
            else:
                end = peaks["end"].values[i]

            # chip_step[start:end] = [np.max(affinities)] * (end - start) # will input max for length of peak
            chip_step_new_2[start:end] = [peaks["score"].values[i]] * (end - start)  # will input 1 for length of peak
            i += 1
    else:
        break

print(peaks)

nuc_peaks = pd.read_csv(f'drive/MyDrive/ML/gcPBM/avg_nuc_calls.txt', sep=',')
nuc_peaks = nuc_peaks[['chr', 'start', 'end', 'score']]
nuc_chr_this = nuc_peaks.loc[nuc_peaks['chr'] == chr_r]
filtered_nuc = nuc_chr_this[(nuc_chr_this["start"] > start_bp) & (nuc_chr_this["end"] <= end_bp)]
nucleosomes = filtered_nuc.copy()

# Set indices starting at our 0
nucleosomes['start'] = nucleosomes['start'] - start_bp
nucleosomes['end'] = nucleosomes['end'] - start_bp
nucleosomes['length'] = nucleosomes['end'] - nucleosomes['start']
nucleosomes = nucleosomes[['start', 'end', 'length', 'score']]

# create step plot from nuc data
nuc_step = [0] * len(affinities)
i = 0
for j in range(len(affinities)):
    if i < len(nucleosomes):
        if j == nucleosomes["start"].values[i]:
            start = nucleosomes["start"].values[i]
            if nucleosomes["end"].values[i] >= len(affinities):
                end = len(affinities)
            else:
                end = nucleosomes["end"].values[i]

            # nuc_step[start:end] = [np.max(affinities)] * (end - start) # will input max for length of peak
            nuc_step[start:end] = [nucleosomes['score'].values[i]] * (end - start)  # will input 1 for length of peak
            i += 1
    else:
        break

# polymerase

pol_peaks_1 = pd.read_csv(f'drive/MyDrive/ML/poly/MACS2_ph0_R1_narrowPeaks.bed', sep='\t', header=None)
pol_peaks_2 = pd.read_csv(f'drive/MyDrive/ML/poly/MACS2_ph0_R2_narrowPeaks.bed', sep='\t', header=None)
pol_peaks_1.rename({0: "chr", 1: "start", 2: "end", 4: "score"}, axis="columns", inplace=True)
pol_peaks_2.rename({0: "chr", 1: "start", 2: "end", 4: "score"}, axis="columns", inplace=True)
pol_peaks_1 = pol_peaks_1[['chr', 'start', 'end', 'score']]
pol_peaks_2 = pol_peaks_2[['chr', 'start', 'end', 'score']]

pol_peaks_this_1 = pol_peaks_1[pol_peaks_1['chr'] == chr_r]
filtered_peaks_1 = pol_peaks_this_1[(pol_peaks_this_1["start"] > start_bp) & (pol_peaks_this_1["end"] <= end_bp)]
peaks_1 = filtered_peaks_1.copy()
peaks_1 = peaks_1.sort_values("start")
peaks_1 = peaks_1.reset_index(drop=True)

pol_peaks_this_2 = pol_peaks_2[pol_peaks_2['chr'] == chr_r]
filtered_peaks_2 = pol_peaks_this_2[(pol_peaks_this_2["start"] > start_bp) & (pol_peaks_this_2["end"] <= end_bp)]
peaks_2 = filtered_peaks_2.copy()
peaks_2 = peaks_2.sort_values("start")
peaks_2 = peaks_2.reset_index(drop=True)

# Set indices starting at our 0
peaks_1['start'] = peaks_1['start'] - start_bp
peaks_1['end'] = peaks_1['end'] - start_bp
peaks_1['length'] = peaks_1['end'] - peaks_1['start']
peaks_1 = peaks_1[['start', 'end', 'length', 'score']]
normalization = max(peaks_1["score"]) - min(peaks_1["score"])
min_val = min(peaks_1["score"])
peaks_1["score"] = list(((peaks_1["score"] - min_val) / normalization))

peaks_2['start'] = peaks_2['start'] - start_bp
peaks_2['end'] = peaks_2['end'] - start_bp
peaks_2['length'] = peaks_2['end'] - peaks_2['start']
peaks_2 = peaks_2[['start', 'end', 'length', 'score']]
normalization = max(peaks_2["score"]) - min(peaks_2["score"])
min_val = min(peaks_2["score"])
peaks_2["score"] = list(((peaks_2["score"] - min_val) / normalization))

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
            peak_step_1[start:end] = [peaks_1['score'].values[i]] * (end - start)  # will input 1 for length of peak
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
            if peaks_2["end"].values[i] >= len(affinities):
                end = len(affinities)
            else:
                end = peaks_2["end"].values[i]

            # nuc_step[start:end] = [np.max(affinities)] * (end - start) # will input max for length of peak
            peak_step_2[start:end] = [peaks_2['score'].values[i]] * (end - start)  # will input 1 for length of peak
            i += 1
    else:
        break

window_start = 56000
window_end = 59000

fig = plt.gcf()
fig.set_size_inches(24, 10)
plt.plot(nuc_step[window_start:window_end], color='red')
plt.plot(affinities[window_start:window_end], color='blue')
plt.plot(chip_step_pb[window_start:window_end], color='green')
plt.plot(chip_step_chip[window_start:window_end], color='orange')
plt.plot(chip_step_new[window_start:window_end], color='purple')
plt.plot(chip_step_new_2[window_start:window_end], color='hotpink')

plt.plot(peak_step_1[window_start:window_end], color='black')
plt.plot(peak_step_2[window_start:window_end], color='darkgray')

plt.xlabel('Base')
plt.ylabel('Relative Intensity')
plt.title(f'{chrsome} - {number - 1}00k-{number}00k subframe')
labels = ["Nucleosomes", 'ML Prediction', 'PB-exo', 'ChIP-exo', 'New data - Rep1', 'New data - Rep2', 'Poly 1',
          'Poly 2']
plt.legend(labels, loc="best")

tipo = 'PB'
window_start = 56000
window_end = 59000
width = 24
height = 4

fig = plt.gcf()
fig.set_size_inches(width, height)
plt.plot(chip_step_pb[window_start:window_end], color='green')
plt.plot(chip_step_chip[window_start:window_end], color='black')
plt.plot(peak_step_1[window_start:window_end], color='darkorange')
plt.plot(chip_step_new[window_start:window_end], color='black')
plt.plot(chip_step_new_2[window_start:window_end], color='black')
plt.plot(peak_step_2[window_start:window_end], color='darkorange')
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
labels = ["Nucleosomes"]
plt.legend(labels, loc="best", prop={'size': 20})
plt.show()
fig.savefig(f'drive/MyDrive/ML/gcPBM/figures/{tipo}/nucleosomes.png')

window_start = 36 * 1000
window_end = 37 * 1000

fig = plt.gcf()
fig.set_size_inches(16, 8)
plt.plot(affinities[window_start:window_end])
plt.plot(nuc_step[window_start:window_end], color='red')
plt.plot(chip_step_pb[window_start:window_end], color='green')
plt.plot(chip_step_chip[window_start:window_end], color='orange')
plt.plot(chip_step_new[window_start:window_end], color='purple')
plt.plot(chip_step_new_2[window_start:window_end], color='hotpink')
labels = ['ML Prediction', 'Nucleosomes', 'PB-exo', 'ChIP-exo (1st Study)', 'ChIP-exo (2nd Study)- Rep1',
          'ChIP-exo (2nd Study)- Rep2']
plt.legend(labels, loc="best")

counts = pd.read_csv(f'drive/MyDrive/ML/gcPBM/{protein}/genome_counts.txt', sep='\t')
counts.sum()

# 121 / 69

tfbs = sequence[267900 - start_bp + 0:267900 - start_bp + 200]

# TODO find out what these are
scores = []
for i in np.arange(len(tfbs) - 5):
    score = 0
    print(i, tfbs[i:i + 6])
    for j in np.arange(0, len(freq_matrix)):
        a = translate[tfbs[i + j]]
        score = score + freq_matrix[j][a]
    scores.append(score)

print([(tfbs[i:i + 6], i, scores[i]) for i in range(len(scores)) if scores[i] >= 4])

top_feat = pd.read_csv(f'drive/MyDrive/ML/top_features.txt', sep='\t', header=None)
freqs = top_feat[0].value_counts()

fig = plt.gcf()
fig.set_size_inches(16, 8)
plt.yticks(range(0, 230, 20))
freqs.plot(kind='barh')
print(freqs)

top_feat = pd.read_csv(f'drive/MyDrive/ML/top_features.txt', sep='\t', header=None)
top_feat_no_pres = pd.read_csv(f'drive/MyDrive/ML/top_features_no_pres.txt', sep='\t', header=None)
top_feat_no_pres_rand = pd.read_csv(f'drive/MyDrive/ML/top_features_rand_no_pres.txt', sep='\t', header=None)

top_feat.rename({0: "Presence"}, axis="columns", inplace=True)
freqs = top_feat['Presence'].value_counts()
top_feat_no_pres.rename({0: "No presence"}, axis="columns", inplace=True)
freqs_np = top_feat_no_pres["No presence"].value_counts()
top_feat_no_pres_rand.rename({0: "No presence - Randomized"}, axis="columns", inplace=True)
freqs_npr = top_feat_no_pres_rand["No presence - Randomized"].value_counts()

df = pd.concat([freqs, freqs_np, freqs_npr], axis=1)
df = df.drop(['Presence'])
print(df)

ax = df.plot(kind='line', xticks=range(len(df)), figsize=(16, 8), rot=90)
ax.set_xticklabels(list(df.index))
ax.set_ylabel('Frequency')

shape_feats = ['SHIFT', 'SLIDE', 'RISE', 'TILT', 'ROLL', 'TWIST']

base_features = pd.read_csv(f'drive/MyDrive/ML/gcPBM/feature_study/all_shape_feats.txt', sep='\t', header=None)
no_SHIFT = pd.read_csv(f'drive/MyDrive/ML/gcPBM/feature_study/SHIFT_feature.txt', sep='\t', header=None)
no_SLIDE = pd.read_csv(f'drive/MyDrive/ML/gcPBM/feature_study/SLIDE_feature.txt', sep='\t', header=None)
no_RISE = pd.read_csv(f'drive/MyDrive/ML/gcPBM/feature_study/RISE_feature.txt', sep='\t', header=None)
no_TILT = pd.read_csv(f'drive/MyDrive/ML/gcPBM/feature_study/TILT_feature.txt', sep='\t', header=None)
no_ROLL = pd.read_csv(f'drive/MyDrive/ML/gcPBM/feature_study/ROLL_feature.txt', sep='\t', header=None)
no_TWIST = pd.read_csv(f'drive/MyDrive/ML/gcPBM/feature_study/TWIST_feature.txt', sep='\t', header=None)

freqs_base = base_features[0].value_counts(normalize=True)
freqs_SHIFT = no_SHIFT[0].value_counts(normalize=True)
freqs_SLIDE = no_SLIDE[0].value_counts(normalize=True)
freqs_RISE = no_RISE[0].value_counts(normalize=True)
freqs_TILT = no_TILT[0].value_counts(normalize=True)
freqs_ROLL = no_ROLL[0].value_counts(normalize=True)
freqs_TWIST = no_TWIST[0].value_counts(normalize=True)

freqs = [freqs_SHIFT, freqs_SLIDE, freqs_RISE, freqs_TILT, freqs_ROLL, freqs_TWIST]

df = pd.DataFrame(index=(freqs_base.index.sort_values()))
for (freq, feat) in zip(freqs, shape_feats):
    df['No ' + feat] = (freq - freqs_base) / freqs_base * 100

print(df)

df = df[['No RISE', 'No ROLL', 'No SHIFT', 'No SLIDE', 'No TILT', 'No TWIST']]
print(base_features)


sns.set()
# max_val = max(abs(df.max().max()), abs(df.min().min()))
t = df.transpose()
f, ax = plt.subplots(figsize=(24, 10))
ax = sns.heatmap(t, annot=True, cmap="RdBu", center=0)

# correlations

# feature-dependent correlations
shape_df = pd.read_csv(f'drive/MyDrive/ML/gcPBM/feature_study/correlations.txt', sep='\t', header=None)
prots = [shape_df[0][j] for j in range(len(shape_df)) if shape_df[1][j] == 'SHIFT']
SHIFT_df = [shape_df[2][j] for j in range(len(shape_df)) if shape_df[1][j] == 'SHIFT']
SLIDE_df = [shape_df[2][j] for j in range(len(shape_df)) if shape_df[1][j] == 'SLIDE']
RISE_df = [shape_df[2][j] for j in range(len(shape_df)) if shape_df[1][j] == 'RISE']
TILT_df = [shape_df[2][j] for j in range(len(shape_df)) if shape_df[1][j] == 'TILT']
ROLL_df = [shape_df[2][j] for j in range(len(shape_df)) if shape_df[1][j] == 'ROLL']
TWIST_df = [shape_df[2][j] for j in range(len(shape_df)) if shape_df[1][j] == 'TWIST']

df = pd.DataFrame(list(zip(prots, SHIFT_df, SLIDE_df, RISE_df, TILT_df, ROLL_df, TWIST_df)))
df = df.set_index(0)

fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(20, 20))

df.rename({1: 'No SHIFT', 2: 'No SLIDE', 3: 'No RISE', 4: 'No TILT', 5: 'No ROLL', 6: 'No TWIST'}, axis="columns",
          inplace=True)

# labels = ['No SHIFT', 'No SLIDE', 'No RISE', 'No TILT', 'No ROLL', 'No TWIST']
# green_diamond = dict(markerfacecolor='g', marker='D')
# bplot = ax.boxplot(t,  flierprops=green_diamond,   patch_artist=True, medianprops=medianprops,labels=labels)

# ax.set_yticks(np.arange(-0.7,1,0.02))
# ax.set_title('Correlations')
# ax.yaxis.grid(True, color='lightgrey')
# ax.set_facecolor('white')

# print(df.mean())
print(df)
vals, names, xs = [], [], []
for i, col in enumerate(df.columns):
    vals.append(df[col].values)
    names.append(col)
    xs.append(
        np.random.normal(i + 1, 0.04, df[col].values.shape[0]))  # adds jitter to the data points - can be adjusted

# labeling = list(df.index)
df.boxplot(ax=ax, showfliers=False, showmeans=True, fontsize=40, rot=45,
           medianprops=dict(linestyle='-', linewidth=3, color='black'),
           meanprops={"marker": "o", "markerfacecolor": "black", "markeredgecolor": "black", "markersize": "10"})
palette = ['r', 'g', 'b', 'y', 'hotpink', 'darkorange']
for x, val, c in zip(xs, vals, palette):
    plt.scatter(x, val, alpha=0.8, color=c)

ax.set_yticks(np.arange(-0.8, 1.1, 0.1))
ax.set_ylabel('$R^2$', size=50)
ax.yaxis.grid(True, color='lightgrey')
ax.set_facecolor('white')

pres_corr = pd.read_csv(f'drive/MyDrive/ML/electro_corr.txt', sep='\t')
pres_corr = pres_corr.set_index('protein')
pres_corr = pres_corr[['0vs5']]
pres_corr.rename({"0vs5": "Presence"}, axis="columns", inplace=True)

np_corr = pd.read_csv(f'drive/MyDrive/ML/no_pres_corr.txt', sep='\t', header=None)
np_corr.rename({0: "protein", 2: "0vs5"}, axis="columns", inplace=True)
np_corr = np_corr.set_index('protein')
np_corr = np_corr[['0vs5']]
np_corr.rename({"0vs5": "No presence"}, axis="columns", inplace=True)

npr_corr = pd.read_csv(f'drive/MyDrive/ML/rand_no_pres_corr.txt', sep='\t', header=None)
npr_corr.rename({0: "protein", 2: "0vs5"}, axis="columns", inplace=True)
npr_corr = npr_corr.set_index('protein')
npr_corr = npr_corr[['0vs5']]
npr_corr.rename({"0vs5": "No presence-Randomized"}, axis="columns", inplace=True)

df = pd.concat([pres_corr, np_corr, npr_corr], axis=1)

ax = df.plot(kind='line', xticks=range(len(df)), yticks=np.arange(0, 1, 0.1), figsize=(20, 8))
ax.set_xticklabels(list(df.index))
ax.set_ylabel('R²')

fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10, 10))

t = df[['Presence', 'No presence']]

vals, names, xs = [], [], []
for i, col in enumerate(t.columns):
    vals.append(t[col].values)
    names.append(col)
    xs.append(np.random.normal(i + 1, 0.04, t[col].values.shape[0]))  # adds jitter to the data points - can be adjusted

labeling = list(t.index)
t.boxplot(ax=ax, showfliers=False, showmeans=True, fontsize=30, rot=45,
          medianprops=dict(linestyle='-', linewidth=3, color='black'),
          meanprops={"marker": "o", "markerfacecolor": "black", "markeredgecolor": "black", "markersize": "10"})
palette = ['r', 'b', 'b', 'y']
for x, val, c in zip(xs, vals, palette):
    plt.scatter(x, val, alpha=0.8, color=c)

ax.set_yticks(np.arange(-0.2, 1.1, 0.1))

ax.yaxis.grid(True, color='lightgrey')
ax.set_ylabel('HT-SELEX $R^2$', size=40)
ax.set_facecolor('white')

pres_corr = pd.read_csv(f'drive/MyDrive/ML/electro_corr_linear.txt', sep=' ', header=None)
pres_corr.rename({0: "protein", 2: "0vs5"}, axis="columns", inplace=True)
pres_corr = pres_corr.set_index('protein')
pres_corr = pres_corr[['0vs5']]
pres_corr.rename({"0vs5": "Presence"}, axis="columns", inplace=True)

np_corr = pd.read_csv(f'drive/MyDrive/ML/no_pres_corr_linear.txt', sep='\t', header=None)
np_corr.rename({0: "protein", 2: "0vs5"}, axis="columns", inplace=True)
np_corr = np_corr.set_index('protein')
np_corr = np_corr[['0vs5']]
np_corr.rename({"0vs5": "No presence"}, axis="columns", inplace=True)

npr_corr = pd.read_csv(f'drive/MyDrive/ML/rand_no_pres_corr_linear.txt', sep='\t', header=None)
npr_corr.rename({0: "protein", 2: "0vs5"}, axis="columns", inplace=True)
npr_corr = npr_corr.set_index('protein')
npr_corr = npr_corr[['0vs5']]
npr_corr.rename({"0vs5": "No presence-Randomized"}, axis="columns", inplace=True)

df = pd.concat([pres_corr, np_corr, npr_corr], axis=1)

ax = df.plot(kind='line', xticks=range(len(df)), yticks=np.arange(-0.6, 1, 0.1), figsize=(20, 8))
ax.set_xticklabels(list(df.index))
ax.set_ylabel('R²')

fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10, 8))
medianprops = dict(linestyle='-', linewidth=4, color='red')

t = df.transpose()
labels = list(t.index)
green_diamond = dict(markerfacecolor='g', marker='D')
bplot = ax.boxplot(t, flierprops=green_diamond, patch_artist=True, medianprops=medianprops, labels=labels)

ax.set_title('Correlations')
ax.yaxis.grid(True, color='lightgrey')
ax.set_facecolor('white')

selex_1 = pd.read_csv(f'drive/MyDrive/ML/correlations/electro_corr_clear.txt', sep=' ')
selex_2 = pd.read_csv(f'drive/MyDrive/ML/correlations/electro_corr_2.txt', sep='\t')
upbm = pd.read_csv(f'drive/MyDrive/ML/correlations/corr_ws_2.txt', sep='\t')
new_upbm = pd.read_csv(f'drive/MyDrive/ML/correlations/new_upbm_correlations.txt', sep='\t')
ratios = pd.read_csv(f'drive/MyDrive/ML/correlations/ratios.txt', sep=' ', header=None)

upbm = pd.concat([upbm, ratios], axis=1)
l = []
for i in range(len(upbm[0])):
    if upbm[0][i] > 0.1:
        l.append(i)
upbm.index = [i for i in range(len(upbm))]
u_pbm = upbm.drop(l)[['PROTEINS', 'FAMILY', 'PRESENCE', 'AVG+DIAG', 'P+A+D']]

new_upbm = pd.concat([new_upbm, ratios], axis=1)
ll = []
for i in range(len(new_upbm[0])):
    if new_upbm[0][i] > 0.1:
        ll.append(i)
new_upbm.index = [i for i in range(len(new_upbm))]
new_u_pbm = new_upbm.drop(ll)[['protein', 'correlation']]

# new_u_pbm = new_upbm.drop(l)[['protein', 'presence','electro','shape']]

u_pbm.rename({'PROTEINS': "protein", 'P+A+D': 'uPBM'}, axis="columns", inplace=True)

selex_1.index = selex_1['protein']
selex_2.index = selex_2['protein']
selex_1 = selex_1['0vs3']
selex_2 = selex_2['0vs5']
selex_full = pd.concat([selex_1, selex_2])

u_pbm.index = u_pbm['protein']
new_u_pbm.index = u_pbm['protein']
u_pbm = u_pbm['uPBM']

gc_pbm = {'gcPBM': [0.905, 0.951, 0.922], 'protein': ['myc', 'mad1', 'max']}
gc_pbm = pd.DataFrame(data=gc_pbm)
gc_pbm.index = gc_pbm['protein']
gc_pbm = gc_pbm['gcPBM']

t = pd.concat([selex_full, u_pbm, gc_pbm], axis=1)
t.rename({0: "HT-SELEX"}, axis="columns", inplace=True)
t = t.sort_index()

sns.set()
# t = df.transpose()
f, ax = plt.subplots(figsize=(5, 10))
ax = sns.heatmap(t[0:30], annot=True, cmap="RdBu", center=0)
f, ax = plt.subplots(figsize=(5, 10))
ax = sns.heatmap(t[30:60], annot=True, cmap="RdBu", center=0)
f, ax = plt.subplots(figsize=(5, 10))
ax = sns.heatmap(t[60:90], annot=True, cmap="RdBu", center=0)
f, ax = plt.subplots(figsize=(5, 10))
ax = sns.heatmap(t[90:], annot=True, cmap="RdBu", center=0)

t.to_csv(r'drive/MyDrive/ML/correlations/all_datasets.txt', sep='\t', mode='w')

(t['uPBM'].sort_values()[:34]).to_csv(r'drive/MyDrive/ML/correlations/upbm_filtered.txt', sep='\t', mode='w')

t.mean(), t.std()

upbm = pd.read_csv(f'drive/MyDrive/ML/correlations/corr_ws_2.txt', sep='\t')
new_upbm = pd.read_csv(f'drive/MyDrive/ML/correlations/new_upbm_correlations.txt', sep='\t')
ratios = pd.read_csv(f'drive/MyDrive/ML/correlations/ratios.txt', sep=' ', header=None)
features = pd.read_csv(f'drive/MyDrive/ML/correlations/upbm_features.txt', sep='\t')

upbm = pd.concat([upbm[['PROTEINS', 'P+A+D']], ratios[0]], axis=1)
upbm.index = upbm['PROTEINS']
new_upbm.index = new_upbm['protein']
features.index = features['protein']
new_upbm = new_upbm.drop('sp140')
features = features.drop('sp140')

DF = pd.concat([upbm, new_upbm, features], axis=1)

l = []
for i in range(len(DF[0])):
    if DF[0][i] > 0.1:
        l.append(i)

DF.index = [i for i in range(len(DF))]
df = DF.drop(l)[['PROTEINS', 'P+A+D', 'correlation', 0, 'presence', 'electro', 'shape']]
f_upbm = list(df.mean()[['presence', 'electro', 'shape']])

gcpbm_feat = pd.read_csv(f'drive/MyDrive/ML/gcPBM/feature_study/gcpbm_features.txt', sep='\t', header=None)
gcpbm_feat.rename({1: "presence", 2: 'electro', 3: 'shape'}, axis="columns", inplace=True)
f_gcpbm = list(gcpbm_feat.mean())

selex_feat = pd.read_csv(f'drive/MyDrive/ML/correlations/htselex_features.txt', sep='\t', header=None)
selex_feat.rename({1: "presence", 2: 'electro', 3: 'shape'}, axis="columns", inplace=True)
f_selex = list(selex_feat.mean())

len(selex_feat)

print(f_upbm, f_gcpbm, f_selex)

labels = ["Presence", "Electrostatic", "Shape"]
# colors = ['red','darkorange','yellow','hotpink','green','blue','cyan','lightgreen','darkgray','purple',
# 'lightblue','olive']
colors = ['red', 'cyan', 'purple']

fig1, ax1 = plt.subplots()
fig1.set_size_inches(15, 15)
ax1.pie(f_upbm, labels=labels, startangle=0, colors=colors, autopct='%1.1f%%', textprops={'fontsize': 25})
ax1.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
plt.show()

fig1, ax1 = plt.subplots()
fig1.set_size_inches(15, 15)
ax1.pie(f_gcpbm, labels=labels, startangle=0, colors=colors, autopct='%1.1f%%', textprops={'fontsize': 25})
ax1.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
plt.show()

fig1, ax1 = plt.subplots()
fig1.set_size_inches(15, 15)
ax1.pie(f_selex, labels=labels, startangle=0, colors=colors, autopct='%1.1f%%', textprops={'fontsize': 25})
ax1.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
plt.show()

# (146 - 20 - 17 - 15) / 1326

selex_1 = pd.read_csv(f'drive/MyDrive/ML/correlations/electro_corr.txt', sep=' ')
selex_2 = pd.read_csv(f'drive/MyDrive/ML/correlations/electro_corr_2.txt', sep='\t')
upbm = pd.read_csv(f'drive/MyDrive/ML/correlations/corr_ws.txt', sep='\t')
ratios = pd.read_csv(f'drive/MyDrive/ML/correlations/ratios.txt', sep=' ', header=None)

upbm = pd.concat([upbm, ratios], axis=1)
l = []
for i in range(len(upbm[0])):
    if upbm[0][i] > 0.1:
        l.append(i)
upbm.index = [i for i in range(len(upbm))]
u_pbm = upbm.drop(l)[['PROTEINS', 'FAMILY', 'PRESENCE', 'AVG+DIAG', 'P+A+D']]

selex_1 = selex_1['0vs3']
selex_2 = selex_2['0vs5']
u_pbm = u_pbm['P+A+D']
gc_pbm = {'gcPBM': [0.905, 0.951, 0.922]}

gc_pbm = pd.DataFrame(data=gc_pbm)
selex_1.mean(), selex_2.mean(), u_pbm.mean(), gc_pbm.mean()

selex_2.mean(), selex_2.std()

print(t)

selex = pd.concat([selex_1, selex_2], axis=0)
# change name col
selex = selex.reset_index(drop=True)
# TODO what is `filtered` and should it be here?
t = pd.concat([selex, u_pbm, gc_pbm, filtered['Filtered']], axis=1)
t.rename({0: "HT-SELEX", 'P+A+D': 'uPBM', 'Filtered': 'HT-SELEX Filtered'}, axis="columns", inplace=True)
t.dropna()

fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(24, 20))

vals, names, xs = [], [], []
for i, col in enumerate(t.columns):
    vals.append(t[col].values)
    names.append(col)
    xs.append(np.random.normal(i + 1, 0.04, t[col].values.shape[0]))  # adds jitter to the data points - can be adjusted

labeling = ['HT-SELEX', 'uPBM', 'gcPBM', 'HT-SELEX Filtered']
bplot = t.boxplot(labels=labeling, ax=ax, showfliers=False, showmeans=True, fontsize=40, rot=45,
                  medianprops=dict(linestyle='-', linewidth=3, color='black'),
                  meanprops={"marker": "o", "markerfacecolor": "black", "markeredgecolor": "black", "markersize": "10"})
# bplot = ax.boxplot(t_1, showfliers=False, showmeans=True,medianprops = dict(linestyle='-', linewidth=3,color='black'),
#                   meanprops={"marker":"o","markerfacecolor":"black", "markeredgecolor":"black","markersize":"10"})
# bplot = ax.boxplot(t_2, showfliers=False,showmeans=True,medianprops = dict(linestyle='-', linewidth=3, color='black'),
#                   meanprops={"marker":"o","markerfacecolor":"black", "markeredgecolor":"black","markersize":"10"})
# bplot = ax.boxplot(t_3, showfliers=False,showmeans=True,medianprops = dict(linestyle='-', linewidth=3, color='black'),
#                   meanprops={"marker":"o","markerfacecolor":"black", "markeredgecolor":"black","markersize":"10"})

palette = ['r', 'g', 'b', 'y']
for x, val, c in zip(xs, vals, palette):
    plt.scatter(x, val, alpha=0.8, color=c)

ax.set_yticks(np.arange(-0.1, 1.01, 0.1))
ax.yaxis.grid(True, color='lightgrey')
ax.set_ylabel('$R^2$', size=40)
ax.set_facecolor('white')

# only rohs's
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(24, 20))

vals, names, xs = [], [], []
rohs_x, rohs_y = [], []
for i, col in enumerate(t.columns):
    vals.append(t[col].values)
    names.append(col)
    xs.append(np.random.normal(i + 1, 0.04, t[col].values.shape[0]))  # adds jitter to the data points - can be adjusted
    rohs_x.append(np.arange(i + 1 - 0.23, i + 1 + 0.24, 0.02))

rohs_y.append([0.67] * len(rohs_x[0]))
rohs_y.append([0.46] * len(rohs_x[0]))
rohs_y.append([0.93] * len(rohs_x[0]))
rohs_y.append([0.67] * len(rohs_x[0]))

dl_y = [0.64] * len(rohs_x[0])

# t_1 = (t[['selex']].dropna()).transpose()
# t_2 = (t[['P+A+D']].dropna()).transpose()
# t_3 = (t[['gcPBM']].dropna()).transpose()
# positions = np.arange(3)+1
labeling = ['HT-SELEX', 'uPBM', 'gcPBM', 'HT-SELEX Filtered']
bplot = t.boxplot(labels=labeling, ax=ax, showfliers=False, showmeans=True, fontsize=40, rot=45,
                  medianprops=dict(linestyle='-', linewidth=3, color='black'),
                  meanprops={"marker": "o", "markerfacecolor": "black", "markeredgecolor": "black", "markersize": "10"})
# bplot = ax.boxplot(t_1, showfliers=False, showmeans=True,medianprops = dict(linestyle='-', linewidth=3,color='black'),
#                   meanprops={"marker":"o","markerfacecolor":"black", "markeredgecolor":"black","markersize":"10"})
# bplot = ax.boxplot(t_2, showfliers=False,showmeans=True,medianprops = dict(linestyle='-', linewidth=3, color='black'),
#                   meanprops={"marker":"o","markerfacecolor":"black", "markeredgecolor":"black","markersize":"10"})
# bplot = ax.boxplot(t_3, showfliers=False,showmeans=True,medianprops = dict(linestyle='-', linewidth=3, color='black'),
#                   meanprops={"marker":"o","markerfacecolor":"black", "markeredgecolor":"black","markersize":"10"})

palette = ['r', 'g', 'b', 'y']
for x, val, c, r_x, r_y in zip(xs, vals, palette, rohs_x, rohs_y):
    plt.scatter(x, val, alpha=0.8, color=c)
    plt.scatter(r_x, r_y, color='darkorange', marker='_')

plt.scatter(rohs_x[1], dl_y, color='navy', marker='_')

ax.set_yticks(np.arange(-0.1, 1.01, 0.1))
ax.yaxis.grid(True, color='lightgrey')
ax.set_ylabel('$R^2$', size=50)
ax.set_facecolor('white')

t = pd.read_csv(f'drive/MyDrive/ML/correlations/upbm_dl.txt', sep='\t', header=None)
t.rename({0: "protein"}, axis="columns", inplace=True)
t.index = t['protein']
t = t[[1, 2, 3, 4, 5, 6, 7]]
t = t.dropna()
len(t)

print(t.columns)

# deep learning folks
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(26, 20))

vals, names, xs = [], [], []

for i, col in enumerate(t.columns):
    vals.append(t[col].values)
    names.append(col)
    xs.append(np.random.normal(i + 1, 0.04, t[col].values.shape[0]))  # adds jitter to the data points - can be adjusted

labeling = ['K-spectrum+shape', 'Dimismatch+shape', 'Deepbind', 'DLBSS', 'CRPT', 'CRPTS', 'DNAffinity']

for k in np.arange(1, 8, 1):
    t.rename({k: labeling[k - 1]}, axis="columns", inplace=True)

bplot = t.boxplot(labels=labeling, ax=ax, showfliers=False, showmeans=True, fontsize=40, rot=45,
                  medianprops=dict(linestyle='-', linewidth=3, color='black'),
                  meanprops={"marker": "o", "markerfacecolor": "black", "markeredgecolor": "black", "markersize": "10"})
# bplot = ax.boxplot(t_1, showfliers=False, showmeans=True,medianprops = dict(linestyle='-', linewidth=3,color='black'),
#                   meanprops={"marker":"o","markerfacecolor":"black", "markeredgecolor":"black","markersize":"10"})
# bplot = ax.boxplot(t_2, showfliers=False,showmeans=True,medianprops = dict(linestyle='-', linewidth=3, color='black'),
#                   meanprops={"marker":"o","markerfacecolor":"black", "markeredgecolor":"black","markersize":"10"})
# bplot = ax.boxplot(t_3, showfliers=False,showmeans=True,medianprops = dict(linestyle='-', linewidth=3, color='black'),
#                   meanprops={"marker":"o","markerfacecolor":"black", "markeredgecolor":"black","markersize":"10"})

palette = ['r', 'g', 'b', 'y', 'darkorange', 'navy', 'olive']
for x, val, c in zip(xs, vals, palette):
    plt.scatter(x, val, alpha=0.8, color=c)

ax.set_yticks(np.arange(0, 1.01, 0.1))
ax.yaxis.grid(True, color='lightgrey')
ax.set_ylabel('$R^2$', size=50)
ax.set_facecolor('white')

microarray = t.dropna()

fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(20, 20))

plt.scatter(microarray['DNAffinity'], microarray['K-spectrum+shape'], alpha=0.8, marker='.', s=1000)
# plt.scatter(microarray['DNAffinity'],microarray['Dimismatch+shape'], alpha=0.8, marker='.', s=1000)
# plt.scatter(microarray['DNAffinity'],microarray['Deepbind'], alpha=0.8, marker='.', s=1000)
# plt.scatter(microarray['DNAffinity'],microarray['DLBSS'], alpha=0.8, marker='.', s=1000)
# plt.scatter(microarray['DNAffinity'],microarray['CRPT'], alpha=0.8, marker='.', s=1000)
# plt.scatter(microarray['DNAffinity'],microarray['CRPTS'], alpha=0.8, marker='.', s=1000)

plt.plot(np.arange(0, 1.01, 1), np.arange(0, 1.01, 1), color='red', linewidth=1)

ax.set_ylabel('K-spectrum+shape', size=50)
# ax.set_ylabel('Dimismatch+shape',size=50)
# ax.set_ylabel('Deepbind',size=50)
# ax.set_ylabel('DLBSS',size=50)
# ax.set_ylabel('CRPT',size=50)
# ax.set_ylabel('CRPTS',size=50)


ax.set_xlabel('DNAffinity', size=50)
ax.set_yticks(np.arange(0, 1.01, 0.1))
ax.set_xticks(np.arange(0, 1.01, 0.1))
ax.yaxis.grid(True, color='lightgrey')
ax.xaxis.grid(True, color='lightgrey')
ax.set_facecolor('white')

# new shape feats
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

(AVG[['SHIFT', 'SLIDE', 'RISE', 'TILT', 'ROLL', 'TWIST', 'TWIST_1', 'TWIST_2']]).to_csv(
    r'drive/MyDrive/ML/input/new_avg.dat', sep=' ', mode='w')

only_shape_corr = pd.read_csv(f'drive/MyDrive/ML/correlations/only_shape_corr.txt', sep=' ')
only_shape_corr.index = only_shape_corr['protein']
only_shape_corr = only_shape_corr['correlation']

new_shape_corr = pd.read_csv(f'drive/MyDrive/ML/correlations/new_avg_corr_2.txt', sep='\t')
new_shape_corr.index = new_shape_corr['protein']
new_shape_corr = new_shape_corr['correlation']

fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(20, 20))

plt.scatter(only_shape_corr, new_shape_corr, alpha=0.8, color='darkorange', marker='.', s=1000)
plt.plot(np.arange(0, 1.01, 1), np.arange(0, 1.01, 1), color='red', linewidth=1)
ax.set_ylabel('New Propeller Twist', size=30)
ax.set_xlabel('Old model', size=30)
ax.set_yticks(np.arange(0, 1.01, 0.1))
ax.set_xticks(np.arange(0, 1.01, 0.1))
ax.yaxis.grid(True, color='lightgrey')
ax.xaxis.grid(True, color='lightgrey')
ax.set_facecolor('white')

fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(15, 10))
plt.plot(new_shape_corr - only_shape_corr)
plt.plot(np.arange(0, 19.01, 19), [0, 0], color='red', linewidth=1)

print(new_shape_corr.mean(), only_shape_corr.mean())
print((new_shape_corr.drop('tbx15')).mean(), (only_shape_corr.drop('tbx15')).mean())

improv = list(new_shape_corr - only_shape_corr)
print('We improve corr on', sum([1 for k in range(19) if improv[k] > 0]))
print('same corr on', sum([1 for k in range(19) if improv[k] == 0.0]))
print('We do not improve corr on', sum([1 for k in range(19) if improv[k] < 0]))

# pca
AVG = pd.read_csv(f'drive/MyDrive/ML/input/avg_tetramer_copy.dat', sep=' ')
AVG = AVG.drop(np.arange(256, 290, 1))
AVG.index = AVG['TETRA']
AVG = AVG[['SHIFT', 'SLIDE', 'RISE', 'TILT', 'ROLL', 'TWIST']]

tetra_fce = {line.split()[0]: np.array([float(x) for x in line.split()[1:]]) for line in
             open(f'drive/MyDrive/ML/input/fce_tetramer.dat') if 'SHIFT' not in line}
diag = pd.DataFrame({tt: tetra_fce[tt][list(range(0, 36, 7))] for tt in tetra_fce.keys()})
diag = diag.transpose()

df = pd.concat([AVG, diag], axis=1)

scaler = StandardScaler()
shape_matrix = scaler.fit_transform(df)

pca = PCA()
# shape_matrix should have 256 rows (all 4-mers) & 12 columns (all features)
transformed_shape = pca.fit_transform(shape_matrix)
# transformed_shape will have the same shape as shape_matrix
print(transformed_shape[:, 0])  # <- this will be the "main" stiffness vector,
# combined from many correlated ones, think of it as the 256 stiffness values e.g. for average shift
print(transformed_shape[:, 1])  # <- this should be the next one importance-wise etc.

print(pca.explained_variance_ratio_)
print(pca.singular_values_)

print(df)

np.sum(list(abs(pca.components_[0]) ** 2))

print(pca.components_[1])

print(pca.components_[2])

print(pca.components_[3])

pcsum = 0
for j in range(len(pca.components_)):
    pcsum = pcsum + abs(pca.components_[j]) * pca.explained_variance_ratio_[j]
    print(j)
print(pcsum)

shape = ['SHIFT', 'SLIDE', 'RISE', 'TILT', 'ROLL', 'TWIST']
average = ['AVG ' + feat for feat in shape]
diag = ['DIAG ' + feat for feat in shape]
labeling = average + diag

colors = ['red', 'darkorange', 'yellow', 'hotpink', 'green', 'blue', 'cyan', 'lightgreen', 'darkgray', 'purple',
          'lightblue', 'olive']

fig1, ax1 = plt.subplots()
fig1.set_size_inches(15, 15)
ax1.pie(pcsum, labels=labeling, startangle=0, colors=colors, autopct='%1.1f%%', textprops={'fontsize': 25})
ax1.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
plt.show()

np.sum(pcsum ** 2)

# (20 + 17) / 146

exclude = []
conf_filter = pd.read_csv(f'drive/MyDrive/ML/conf_filter.txt', sep=',')
rohs_protz = pd.read_csv(f'drive/MyDrive/ML/correlations/electro_corr.txt', sep=' ')
for i in range(len(conf_filter)):
    if conf_filter['confidence'][i] == 0:
        exclude.append(conf_filter['protein'][i])

rohs_protz.index = rohs_protz['protein']
filtered = rohs_protz.drop(exclude)

rohs_protz.rename({'0vs3': 'All proteins'}, axis="columns", inplace=True)
filtered.rename({'0vs3': 'Filtered'}, axis="columns", inplace=True)
t = pd.DataFrame([rohs_protz['All proteins'], filtered['Filtered']])
t = t.transpose()

fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10, 10))

vals, names, xs = [], [], []
for i, col in enumerate(t.columns):
    vals.append(t[col].values)
    names.append(col)
    xs.append(np.random.normal(i + 1, 0.04, t[col].values.shape[0]))  # adds jitter to the data points - can be adjusted

labeling = ['All proteins', 'Filtered']
bplot = t.boxplot(labels=labeling, ax=ax, showfliers=False, showmeans=True, fontsize=40, rot=45,
                  medianprops=dict(linestyle='-', linewidth=3, color='black'),
                  meanprops={"marker": "o", "markerfacecolor": "black", "markeredgecolor": "black", "markersize": "10"})
# bplot = ax.boxplot(t_1, showfliers=False, showmeans=True,medianprops = dict(linestyle='-', linewidth=3,color='black'),
#                   meanprops={"marker":"o","markerfacecolor":"black", "markeredgecolor":"black","markersize":"10"})
# bplot = ax.boxplot(t_2, showfliers=False,showmeans=True,medianprops = dict(linestyle='-', linewidth=3, color='black'),
#                   meanprops={"marker":"o","markerfacecolor":"black", "markeredgecolor":"black","markersize":"10"})
# bplot = ax.boxplot(t_3, showfliers=False,showmeans=True,medianprops = dict(linestyle='-', linewidth=3, color='black'),
#                   meanprops={"marker":"o","markerfacecolor":"black", "markeredgecolor":"black","markersize":"10"})

palette = ['r', 'g', 'b', 'y']
for x, val, c in zip(xs, vals, palette):
    plt.scatter(x, val, alpha=0.8, color=c)

ax.set_yticks(np.arange(-0.1, 1.01, 0.1))
ax.yaxis.grid(True, color='lightgrey')
ax.set_ylabel('HT-SELEX $R^2$', size=50)
ax.set_facecolor('white')

t.transpose()

print(conf_filter)

print(filtered)
