######################################################################
# Analysis of the SNP matrix that might inform us on how to parse the 
# matrix to define haplotypes in the most parsimonious way. 
# PSB
######################################################################

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Info on haplotypes and dataset
snp_matrix_path = '/home/ada/Desktop/Shraiman_lab/9_clade_analysis/haplotyping/fixed_9_clade_snp_matrix.csv'
clade_9_strains = ['PB87','PB40','GGACTCCT-AGAGTAGA','GGACTCCT-GTAAGGAG','GGACTCCT-CTAAGCCT','TCCTGAGC-TAGATCGC','TCCTGAGC-CTCTCTAT','PB73','PB80']

# Read in the data matrix
snp_matrix = pd.read_csv(snp_matrix_path, index_col=0)
columns = snp_matrix.columns.tolist()
strains = snp_matrix.index.tolist()


## What is the overall distance distribution between the columns of the matrix? 
all_distances = []

for i in range(0, len(columns)-1):
    dist = int(float(columns[i+1])) - int(float(columns[i]))
    all_distances.append(dist)

"""
plt.hist(all_distances, bins=range(0,max(all_distances)+100, 100))
plt.yscale('log')
plt.xlabel('distance (bp)')
plt.ylabel('counts (log)')
plt.title('What is the overall distribution of distances between the consecutive columns of the PSB SNP matrix?')
plt.show()
"""

## What is the distance between each green dot in the SNP matrix? (use all strains except for the 9 clade)
distances_green_green = []

for i in range(1, len(columns)-1):
    for s in range(1, len(strains)):
        if strains[s]  not in clade_9_strains:
            current_snp = snp_matrix.iloc[s, i]
            if current_snp  == snp_matrix.iloc[s, i+1]: 
                dist = int(float(columns[i+1])) - int(float(columns[i]))
                distances_green_green.append(dist)

"""
plt.hist(distances, bins=range(0,max(distances)+100, 100))
plt.yscale('log')
plt.xlabel('distance (bp)')
plt.ylabel('counts (log)')
plt.title('What is the distribution of genomic distances between consecutive green loci? (green -> green)')
plt.show()
"""

# Calculate the cumulant for the green -> white and white -> green distances of F clade
sizes, bins, _ = plt.hist(distances_green_green, bins=range(0,40000, 100), edgecolor='black')
all_sizes, bins, _ = plt.hist(all_distances, bins=range(0,40000, 100), edgecolor='black')
sizes = np.array(sizes)
total = float(sum(sizes))

fractions_gg = [0.0]
current_total = 0.0
for i in sizes:
    current_total += (i/total)
    fractions_gg.append(current_total)

bins=range(100,40100, 100)

## What about the distance between WT/green transitions ? (use all strains except for the 9 clade)
transition_distances = []

for i in range(1, len(columns)-1):
    for s in range(1, len(strains)):
        if strains[s]  not in clade_9_strains:
            current_snp = snp_matrix.iloc[s, i]
            if (current_snp == 1  and snp_matrix.iloc[s, i+1] != 1) or (current_snp != 1  and snp_matrix.iloc[s, i+1] == 1): 
                dist = int(float(columns[i+1])) - int(float(columns[i]))
                transition_distances.append(dist)


# Calculate the cumulant for the green -> white and white -> green distances of F clade
sizes, bins, _ = plt.hist(transition_distances, bins=range(0,40000, 100), edgecolor='black')
sizes = np.array(sizes)
total = float(sum(sizes))

fractions_transitions = [0.0]
current_total = 0.0
for i in sizes:
    current_total += (i/total)
    fractions_transitions.append(current_total)

bins=range(100,40100, 100)
plt.clf()

plt.plot(bins, fractions_gg, color='green', label='same color transitions (g->g, w->w, etc)')

plt.plot(bins, fractions_transitions, color='blue', label='WT/green transitions')

# Calculate the cumulant for all successive loci
all_sizes = np.array(all_sizes)
all_total = float(sum(all_sizes))

all_fractions = [0.0]
all_current_total = 0.0
for i in all_sizes:
    all_current_total += (i/all_total)
    all_fractions.append(all_current_total)


plt.plot(bins, all_fractions, color='red', label='all transitions')
plt.xlabel('distance between columns (bp) (log)')
plt.ylabel('fraction of column pairs')
plt.title('Is there a significant difference between same color  \ntransitions vs all  vs WT/green transitions distances?')
plt.legend()
plt.xscale('log')
plt.show()

"""
## How does the distance across the chromosome change between SNP columns?
# Make a list of all the pairs that have a very large real genomic distance between them.
distances = []
current_chr = snp_matrix.iloc[0, 0]   # [row_position, column_position]
pairs = []

for i in range(0, len(columns)-1):
    if current_chr == snp_matrix.iloc[0, i+1]:
        dist = int(float(columns[i+1])) - int(float(columns[i]))
        distances.append(dist)
        if dist > 5000 and dist < 7000:
            pairs.append([columns[i], columns[i+1], current_chr])

    current_chr = snp_matrix.iloc[0, i+1]

for i in pairs:
    print(i)

plt.plot(range(0,201-1), distances)
plt.xlabel('Column number of the SNP matrix')
plt.ylabel('Distance between two SNP columns (bp)')
plt.title('How does the distance across the chromosome change between SNP columns?')
#plt.show()


#for i in pairs:
#    print(i)

## Look at all pairs of successive loci containing a transition from green to white or 
# white to green and plot a cumulative histogram of inter-loci distance in this set.

f_matrix = snp_matrix.loc[f_clade_strains, :]
transition_distances = []
all_distances = []
current_chr = snp_matrix.iloc[0, 0]   # [row_position, column_position]
f_n = len(f_clade_strains)

for i in range(1, len(columns)-1):

    if current_chr == snp_matrix.iloc[0, i+1]:   # check if the next column belong to the same chr
        col1 = np.array(f_matrix.iloc[:, i])
        col2 = np.array(f_matrix.iloc[:, i+1])

        if sum(col1==1) >= f_n/2 and sum(col2==1) >= f_n/2:    # Check if the current columns belong to the F clade fixed matrix

            for s in range(1, len(strains)):                # Go through the whole matrix, not just the F clade strains
                current_snp = snp_matrix.iloc[s, i]
                next_snp = snp_matrix.iloc[s, i+1]

                if (current_snp == -1 and next_snp != -1) or (current_snp != -1 and next_snp == -1):   # Check if we transition from one color to another.
                    dist = int(float(columns[i+1])) - int(float(columns[i]))
                    transition_distances.append(dist)
                elif current_snp == next_snp:
                    dist = int(float(columns[i+1])) - int(float(columns[i]))
                    all_distances.append(dist)

    current_chr = snp_matrix.iloc[0, i+1]

# Compare the distribution of all distances between successive loci vs transitions of green -> white and white -> green
# Make a cumulant graph.

# Calculate the cumulant for the green -> white and white -> green distances of F clade
sizes, bins, _ = plt.hist(transition_distances, bins=range(0,40000, 100), edgecolor='black')
all_sizes, bins, _ = plt.hist(all_distances, bins=range(0,40000, 100), edgecolor='black')
sizes = np.array(sizes)
total = float(sum(sizes))

fractions = [0.0]
current_total = 0.0
for i in sizes:
    current_total += (i/total)
    fractions.append(current_total)

bins=range(100,40100, 100)
plt.clf()
#plt.plot(bins, fractions, color='green', label='WT/green transitions')

# Calculate the cumulant for all successive loci
all_sizes = np.array(all_sizes)
all_total = float(sum(all_sizes))

all_fractions = [0.0]
all_current_total = 0.0
for i in all_sizes:
    all_current_total += (i/all_total)
    all_fractions.append(all_current_total)


plt.plot(bins, all_fractions, color='red', label='same color SNP pairs \n(green->green, wt->wt, etc)')
plt.xlabel('distance between columns (bp) (log)')
plt.ylabel('fraction of column pairs')
plt.title('Is there a significant difference between WT/green  \ntransition genomic distances vs distances that stays the same color in \n F clade columns of the SNP matrix?')
plt.legend()
plt.xscale('log')
#plt.show()

all_distances = all_distances + transition_distances
all_sizes, bins, _ = plt.hist(all_distances, bins=range(0,40000, 100), edgecolor='black')
transition_sizes, bins, _ = plt.hist(transition_distances, bins=range(0,40000, 100), edgecolor='black')

fractions = np.array(transition_sizes)
all_fractions = np.array(all_sizes)
fractions = fractions + 1
all_fractions = all_fractions +1
divided = fractions / all_fractions
one_minus = 1 - divided

neg_log = -1 * np.log(one_minus)

bins=range(100,40000, 100)
plt.clf()
plt.plot(bins, neg_log )
plt.show()
"""