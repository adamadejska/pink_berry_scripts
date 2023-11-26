#################################################################################
# Now that we have identified different horizontal transfer regions, 
# we want to see if there is a correlation between the number of transfers 
# and the "age" of each strain (aka. number of singletons which we consider
# the good indication of some time scale since the strain diverged from the
# origin)
# PSB 
#################################################################################

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import random
from scipy.stats import pearsonr



# Files and variables used
full_matrix_path = '/home/ada/Desktop/Shraiman_lab/data/snp_data_2020.csv'

clade_9_strains = ['PB87','PB40','GGACTCCT-AGAGTAGA','GGACTCCT-GTAAGGAG','GGACTCCT-CTAAGCCT','TCCTGAGC-TAGATCGC','TCCTGAGC-CTCTCTAT','PB73','PB80']
strains = ['PB87','PB40','GGACTCCT-AGAGTAGA','GGACTCCT-GTAAGGAG','GGACTCCT-CTAAGCCT','TCCTGAGC-TAGATCGC','TCCTGAGC-CTCTCTAT','PB73','PB80','AAGAGGCA-TATCCTCT','PB63','AAGAGGCA-TAGATCGC','PB39','PB_8','GCTACGCT-GTAAGGAG','PB24','PB55','AAGAGGCA-GTAAGGAG','AAGAGGCA-AAGGAGTA','GCTACGCT-ACTGCATA','PB31','AAGAGGCA-CTAAGCCT','AAGAGGCA-CTCTCTAT','PB32','PB16','PB47','PB48','PB64','AGGCAGAA-ACTGCATA','AGGCAGAA-CTAAGCCT','TAAGGCGA-TATCCTCT','PB78','PB61','GGACTCCT-ACTGCATA','GGACTCCT-TAGATCGC','PB45','PB77','PB69','PB_5','PB37','PB13','PB85','PB28','PB67','PB59','PB52','PB26','PB34','PB11','CTCTCTAC-AAGGAGTA','GTAGAGGA-CTAAGCCT','CTCTCTAC-CTAAGCCT','CTCTCTAC-CTCTCTAT','CTCTCTAC-TAGATCGC','CTCTCTAC-ACTGCATA','CAGAGAGG-TAGATCGC','CAGAGAGG-CTAAGCCT','CAGAGAGG-GTAAGGAG','CTCTCTAC-TATCCTCT','CAGAGAGG-ACTGCATA','CAGAGAGG-AGAGTAGA','CAGAGAGG-TATCCTCT','CAGAGAGG-AAGGAGTA','CTCTCTAC-AGAGTAGA','CTCTCTAC-GTAAGGAG','PB44','PB36','PB58','PB82','PB20','PB60','PB51','PB_2','PB25','PB41','CGTACTAG-TAGATCGC','PB57','PB43','PB12','PB66','PB90','PB21','GTAGAGGA-CTCTCTAT','PB75']
"""
# Load in the SNP matrix
snp_matrix = pd.read_csv(full_matrix_path, index_col=0)
columns = snp_matrix.columns.tolist()
strains = snp_matrix.index.tolist()

print(strains)

singleton_counts = {}
# First, determine the number of singletons per mixing layer strain
for strain in strains:
    print(strain)
    singleton_counter = 0
    # For each SNP position (column), determine if it is a singleton for that particular strain
    for c in columns:
        snp_col = np.array(snp_matrix.loc[strains,c])
        if snp_matrix.loc[strain, c] == 1 and np.nansum(snp_col) == 1:
            singleton_counter += 1

    singleton_counts[strain] = singleton_counter

print(singleton_counts)
"""
singleton_counts = {'PB63': 2826,'PB28': 1369, 'GGACTCCT-AGAGTAGA': 699, 'AAGAGGCA-AGAGTAGA': 966, 'PB73': 767, 'PB76': 1250, 'AAGAGGCA-ACTGCATA': 1102, 'GGACTCCT-GTAAGGAG': 711, 'GCTACGCT-GTAAGGAG': 931, 'TCCTGAGC-CTCTCTAT': 597, 'PB55': 855, 'PB_5': 653, 'PB31': 997, 'AAGAGGCA-TATCCTCT': 809, 'PB24': 655, 'PB87': 398, 'AAGAGGCA-AAGGAGTA': 813, 'PB84': 782, 'PB40': 366, 'AGGCAGAA-ACTGCATA': 951, 'AGGCAGAA-CTAAGCCT': 886, 'PB67': 772, 'PB_8': 674, 'PB80': 680, 'PB39': 470, 'PB47': 576, 'PB59': 607, 'PB37': 116, 'PB61': 289, 'TCCTGAGC-TAGATCGC': 504, 'PB77': 150, 'PB69': 320, 'TAAGGCGA-AAGGAGTA': 587, 'PB13': 98, 'PB82': 336, 'PB48': 432, 'PB52': 497, 'GGACTCCT-TAGATCGC': 96, 'PB11': 189, 'GGACTCCT-CTAAGCCT': 551, 'PB44': 351, 'PB36': 397, 'AAGAGGCA-GTAAGGAG': 358, 'PB16': 344, 'PB32': 271, 'PB20': 492, 'GCTACGCT-ACTGCATA': 337, 'PB85': 131, 'PB65': 444, 'PB58': 297, 'PB64': 304, 'CTCTCTAC-GTAAGGAG': 266, 'PB60': 279, 'TAAGGCGA-GTAAGGAG': 315, 'AGGCAGAA-GTAAGGAG': 101, 'PB21': 132, 'TAAGGCGA-TATCCTCT': 194, 'AAGAGGCA-CTCTCTAT': 188, 'PB78': 129, 'PB_1': 86, 'PB91': 120, 'AAGAGGCA-TAGATCGC': 180, 'CTCTCTAC-AAGGAGTA': 100, 'GTAGAGGA-ACTGCATA': 96, 'CAGAGAGG-AAGGAGTA': 67, 'CTCTCTAC-AGAGTAGA': 74, 'GTAGAGGA-GTAAGGAG': 96, 'PB43': 93, 'PB34': 179, 'AAGAGGCA-CTAAGCCT': 182, 'PB92': 125, 'CTCTCTAC-TATCCTCT': 76, 'AGGCAGAA-AAGGAGTA': 65, 'PB90': 101, 'PB12': 98, 'CAGAGAGG-TATCCTCT': 239, 'PB18': 87, 'PB26': 139, 'GTAGAGGA-AAGGAGTA': 61, 'PB41': 37, 'AGGCAGAA-AGAGTAGA': 179, 'PB75': 104, 'PB45': 117, 'CAGAGAGG-TAGATCGC': 85, 'TAAGGCGA-AGAGTAGA': 67, 'CAGAGAGG-GTAAGGAG': 81, 'PB17': 55, 'CAGAGAGG-ACTGCATA': 96, 'CTCTCTAC-CTCTCTAT': 62, 'PB51': 60, 'CGTACTAG-AGAGTAGA': 55, 'PB66': 97, 'PB25': 59, 'PB83': 75, 'CGTACTAG-CTAAGCCT': 54, 'CTCTCTAC-TAGATCGC': 77, 'TCCTGAGC-AAGGAGTA': 218, 'CAGAGAGG-AGAGTAGA': 107, 'PB_3': 66, 'AGGCAGAA-CTCTCTAT': 80, 'CAGAGAGG-CTAAGCCT': 63, 'TAAGGCGA-TAGATCGC': 180, 'PB81': 42, 'CGTACTAG-ACTGCATA': 36, 'CGTACTAG-AAGGAGTA': 61, 'PB89': 36, 'TAAGGCGA-ACTGCATA': 17, 'PB74': 26, 'CTCTCTAC-ACTGCATA': 79, 'PB57': 40, 'PB49': 39, 'GTAGAGGA-CTAAGCCT': 89, 'CGTACTAG-CTCTCTAT': 39, 'GTAGAGGA-CTCTCTAT': 54, 'PB_9': 36, 'PB_4': 26, 'GTAGAGGA-AGAGTAGA': 36, 'GTAGAGGA-TATCCTCT': 44, 'PB29': 29, 'PB10': 42, 'PB_2': 29, 'AGGCAGAA-TAGATCGC': 35, 'CGTACTAG-TATCCTCT': 32, 'PB53': 32, 'GGACTCCT-ACTGCATA': 410, 'PB19': 22, 'CTCTCTAC-CTAAGCCT': 67, 'TCCTGAGC-AGAGTAGA': 30, 'CAGAGAGG-CTCTCTAT': 42, 'PB35': 24, 'PB27': 18, 'GTAGAGGA-TAGATCGC': 26, 'CGTACTAG-TAGATCGC': 24, 'TCCTGAGC-CTAAGCCT': 13, 'TAAGGCGA-CTAAGCCT': 30, 'TCCTGAGC-TATCCTCT': 19, 'PB42': 24, 'TCCTGAGC-GTAAGGAG': 17, 'TAAGGCGA-CTCTCTAT': 18, 'TCCTGAGC-ACTGCATA': 23, 'AGGCAGAA-TATCCTCT': 22, 'CGTACTAG-GTAAGGAG': 21}

# Next, read in the file with the HGT blocks for each strain (disctionary-like format)
hgt_blocks_file_path25 = '/home/ada/Desktop/Shraiman_lab/9_clade_analysis/haplotyping/hgt_blocks_lengths_all_strains_2_5.txt'
hgt_blocks_file_path35 = '/home/ada/Desktop/Shraiman_lab/9_clade_analysis/haplotyping/hgt_blocks_lengths_all_strains_3_5.txt'
hgt_blocks_file_path15 = '/home/ada/Desktop/Shraiman_lab/9_clade_analysis/haplotyping/hgt_blocks_lengths_all_strains_1_5.txt'
hgt_blocks_file_path55 = '/home/ada/Desktop/Shraiman_lab/9_clade_analysis/haplotyping/hgt_blocks_lengths_all_strains_5_5.txt'

hgt_green_blocks_files = [hgt_blocks_file_path25, hgt_blocks_file_path35, hgt_blocks_file_path15,hgt_blocks_file_path55]

"""
#fig, axs = plt.subplots(2, 2)
for i in range(0, len(hgt_green_blocks_files)):

    hgt_blocks_green = {}

    with open(hgt_green_blocks_files[i], 'r') as f:
        for line in f:
            line = line.split(':')
            strain, blocks = line[0], [int(i) for i in line[1].split(',')]
            blocks = [i for i in blocks if i != 0]   # clean zeros

            # Get rid of any strains with less than 10 blocks
            if len(blocks) > 10:
                hgt_blocks_green[strain] = blocks


    # Plot the number of HGT blocks against the number of singleton SNPs for each mixing layer strain
    # Differentiate between green and yellow blocks

    x_g, y_g = [], []
    for strain in strains:
        if strain not in clade_9_strains and strain in hgt_blocks_green.keys():
            total = len(hgt_blocks_green[strain])
            if singleton_counts[strain] < 1000 and total < 120:
                x_g.append(total)
                y_g.append(singleton_counts[strain])

    print(hgt_green_blocks_files[i])
    print(len(hgt_blocks_green['AAGAGGCA-TAGATCGC']))
    a, b = np.polyfit(y_g, x_g, 1)

    l = hgt_green_blocks_files[i].split('_')[-1].split('.')[0]
    penalty = hgt_green_blocks_files[i].split('_')[-2]
    label = 'penalty: ' + penalty + ' L: ' + l

    y_g = np.array(y_g)

    #add line of best fit to plot
    plt.plot(y_g, a*y_g+b, label = label)
    #plt.scatter(y_g, x_g)
 
    if i == 1:
        axs[0, 0].scatter(y_g, x_g)
        axs[0, 0].plot(y_g, a*y_g+b, label = label, c='red')
        axs[0, 0].set_title(label)
    elif i == 2:
        axs[0, 1].scatter(y_g, x_g)
        axs[0, 1].plot(y_g, a*y_g+b, label = label, c='red')
        axs[0, 1].set_title(label)
    elif i == 3:
        axs[1, 0].scatter(y_g, x_g)
        axs[1, 0].plot(y_g, a*y_g+b, label = label, c='red')
        axs[1, 0].set_title(label)
    else:
        axs[1,1].scatter(y_g, x_g)
        axs[1,1].plot(y_g, a*y_g+b, label = label, c='red')
        axs[1,1].set_title(label)

plt.title('Number of HGT blocks vs number of singletons in the mixing layer strains (PSB)\n best fit line only')
plt.xlabel('Singleton number')
plt.ylabel('HGT blocks number')
plt.legend()
plt.show()
    

# Which strains show the most variation in their number of HGT blocks given different variations of the penalty or L? 
# We should look at those strains more closely and determine which values to use for further analysis.

block_numbers_per_strain = {}
variance = []

for i in range(0, len(hgt_green_blocks_files)):

    with open(hgt_green_blocks_files[i], 'r') as f:
        for line in f:
            line = line.split(':')
            strain, blocks = line[0], [int(i) for i in line[1].split(',')]
            blocks = [i for i in blocks if i != 0]   # clean zeros
            if strain not in block_numbers_per_strain.keys():
                block_numbers_per_strain[strain] = [len(blocks)]
            else:
                block_numbers_per_strain[strain].append(len(blocks))

print(block_numbers_per_strain)
# Calculate variance of the number of blocks per strain
for s in block_numbers_per_strain.keys():
    variance.append(np.std(np.array(block_numbers_per_strain[s])))

plt.bar(list(block_numbers_per_strain.keys()), variance)
plt.xticks(rotation=90)
plt.ylabel('Standard deviation')
plt.show()
"""
# Check the pearson's correlation for singletons vs HGT block number
# # Use only one file 

hgt_blocks_green = {}

with open(hgt_blocks_file_path15, 'r') as f:
    for line in f:
        line = line.split(':')
        strain, blocks = line[0], [int(i) for i in line[1].split(',')]
        blocks = [i for i in blocks if i != 0]   # clean zeros
        # Get rid of any strains with less than 10 blocks
        #if len(blocks) > 10:
        hgt_blocks_green[strain] = blocks




# Plot the number of HGT blocks against the number of singleton SNPs for each mixing layer strain
# Differentiate between green and yellow blocks

x_g, y_g = [], []
for strain in strains:
    if strain not in clade_9_strains and strain in hgt_blocks_green.keys():
        total = len(hgt_blocks_green[strain])
        if singleton_counts[strain] < 1500:
            x_g.append(total)
            y_g.append(singleton_counts[strain])
    
y_g = np.array(y_g)
a, b = np.polyfit(y_g, x_g, 1)
plt.plot(y_g, a*y_g+b)
plt.scatter(y_g, x_g) 
plt.suptitle('Number of HGT blocks vs number of singletons in the mixing layer strains (PSB): \n varying penalty variables + best fit line')
plt.xlabel('Singleton number')
plt.ylabel('HGT blocks number')
plt.legend()
plt.show()

corr,_ = pearsonr(y_g, x_g)
#corr = (np.cov(f_freq,ab_freq)/(np.std(f_freq)*np.std(ab_freq)))
print('Correlation between the experimental frequencies is: %.3f' %corr)

random_correlations = []
# Check if the correlation is significant by resampling Y axis coordinates (AB_frequencies) a bunch of times
for i in range(0,100000):
	x_g_random = x_g[:]
	random.shuffle(x_g_random)
	# Calculate Pearson's correlation
	corr,_ = pearsonr(y_g, x_g_random)
	random_correlations.append(corr)

print(np.mean(random_correlations))
print(np.median(random_correlations))
print(np.std(random_correlations))
bins = list(np.arange(-13,16)/20)

# Get the heights of each bin in the histogram (for probability calculations)
counts, bin_edges = np.histogram(random_correlations, bins=bins)
print(bin_edges)
print(counts)

# Make the actual histogram
plt.hist(random_correlations, bins=bins, edgecolor='black')
plt.title('Pearson\'s correlation coefficient calculated 100k times by resampling of HGT block numbers')
plt.xlabel('Pearson\'s correlation coefficient')
plt.ylabel('Frequency')
plt.show()
