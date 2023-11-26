#################################################################################
# Now that we have identified different horizontal transfer regions, 
# we want to see if there is a correlation between the number of transfers 
# and the "age" of each strain (aka. number of singletons which we consider
# the good indication of some time scale since the strain diverged from the
# origin)
#################################################################################

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import random
from scipy.stats import pearsonr



# Files and variables used
full_matrix_path = '/home/ada/Desktop/Shraiman_lab/srb_data/srb_snp_data_2021_chr.csv'

basal_AB_strains = ['CGTACTAG-TATCCTCT','CGTACTAG-GTAAGGAG','PB50','AGGCAGAA-TATCCTCT','GTAGAGGA-TAGATCGC','TAAGGCGA-CTAAGCCT','PB82','PB90','PB27','AGGCAGAA-GTAAGGAG','CGTACTAG-CTAAGCCT','TAAGGCGA-TATCCTCT','PB18','CGTACTAG-AAGGAGTA','GTAGAGGA-AGAGTAGA','PB10','PB34','PB76']
f_clade_strains = ['CAGAGAGG-AAGGAGTA','CAGAGAGG-GTAAGGAG','CTCTCTAC-TATCCTCT','AAGAGGCA-CTCTCTAT','AAGAGGCA-TATCCTCT','AAGAGGCA-CTAAGCCT','GCTACGCT-ACTGCATA','PB39','PB32','PB40','PB88','PB16','PB79','PB87','AAGAGGCA-TAGATCGC','AAGAGGCA-GTAAGGAG','PB64']
topAB_strains = ['TAAGGCGA-CTCTCTAT','TAAGGCGA-ACTGCATA','PB83','AGGCAGAA-CTCTCTAT','AGGCAGAA-TAGATCGC','PB53','GTAGAGGA-CTCTCTAT','TCCTGAGC-ACTGCATA','PB58','TCCTGAGC-AGAGTAGA','PB33','PB_4','PB84','TCCTGAGC-GTAAGGAG','PB41','PB25','GTAGAGGA-TATCCTCT']
mixing_layer_strains = ['CAGAGAGG-TAGATCGC','CAGAGAGG-AGAGTAGA','CAGAGAGG-ACTGCATA','CTCTCTAC-CTAAGCCT','CTCTCTAC-ACTGCATA','PB37','PB_5','CAGAGAGG-TATCCTCT','CTCTCTAC-AGAGTAGA','CAGAGAGG-CTCTCTAT','CTCTCTAC-CTCTCTAT','PB11','CTCTCTAC-TAGATCGC','GGACTCCT-TAGATCGC','GGACTCCT-ACTGCATA','PB77','PB29','PB35','PB45','TCCTGAGC-CTAAGCCT','PB85','CAGAGAGG-CTAAGCCT','GGACTCCT-CTAAGCCT','PB57','PB21','TAAGGCGA-CTCTCTAT']
strains = ['CAGAGAGG-AAGGAGTA','CAGAGAGG-GTAAGGAG','CTCTCTAC-TATCCTCT','AAGAGGCA-CTCTCTAT','AAGAGGCA-TATCCTCT','AAGAGGCA-CTAAGCCT','GCTACGCT-ACTGCATA','PB39','PB32','PB40','PB88','PB16','PB79','PB87','AAGAGGCA-TAGATCGC','AAGAGGCA-GTAAGGAG','PB64','CAGAGAGG-TAGATCGC','CAGAGAGG-AGAGTAGA','CAGAGAGG-ACTGCATA','CTCTCTAC-CTAAGCCT','CTCTCTAC-ACTGCATA','PB37','PB_5','CAGAGAGG-TATCCTCT','CTCTCTAC-AGAGTAGA','CAGAGAGG-CTCTCTAT','CTCTCTAC-CTCTCTAT','PB11','CTCTCTAC-TAGATCGC','GGACTCCT-TAGATCGC','GGACTCCT-ACTGCATA','PB77','PB29','PB35','PB45','TCCTGAGC-CTAAGCCT','PB85','CAGAGAGG-CTAAGCCT','GGACTCCT-CTAAGCCT','PB57','PB21','TAAGGCGA-CTCTCTAT','TAAGGCGA-ACTGCATA','PB83','AGGCAGAA-CTCTCTAT','AGGCAGAA-TAGATCGC','PB53','GTAGAGGA-CTCTCTAT','TCCTGAGC-ACTGCATA','PB58','TCCTGAGC-AGAGTAGA','PB33','PB_4','PB84','TCCTGAGC-GTAAGGAG','PB41','PB25','GTAGAGGA-TATCCTCT','CGTACTAG-TATCCTCT','CGTACTAG-GTAAGGAG','PB50','AGGCAGAA-TATCCTCT','GTAGAGGA-TAGATCGC','TAAGGCGA-CTAAGCCT','PB82','PB90','PB27','AGGCAGAA-GTAAGGAG','CGTACTAG-CTAAGCCT','TAAGGCGA-TATCCTCT','PB18','CGTACTAG-AAGGAGTA','GTAGAGGA-AGAGTAGA','PB10','PB34','PB76','CGTACTAG-ACTGCATA','TCCTGAGC-TATCCTCT','CGTACTAG-AGAGTAGA','CGTACTAG-TAGATCGC','AGGCAGAA-AGAGTAGA','PB66','TAAGGCGA-AAGGAGTA','PB43','CGTACTAG-CTCTCTAT','PB_3','PB44','AGGCAGAA-AAGGAGTA','PB42']

"""
# Load in the SNP matrix
snp_matrix = pd.read_csv(full_matrix_path, index_col=0)
columns = snp_matrix.columns.tolist()
strains = snp_matrix.index.tolist()[1:]


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
"""
#print(singleton_counts)
singleton_counts = {'AAGAGGCA-AAGGAGTA': 483, 'PB36': 481, 'PB80': 646, 'PB23': 607, 'PB78': 547, 'TCCTGAGC-AAGGAGTA': 413, 'PB47': 454, 'PB48': 481, 'AAGAGGCA-ACTGCATA': 423, 'PB60': 549, 'PB_1': 554, 'PB91': 427, 'PB95': 523, 'PB19': 423, 'PB59': 584, 'GTAGAGGA-ACTGCATA': 465, 'PB81': 481, 'PB28': 545, 'PB68': 501, 'GTAGAGGA-GTAAGGAG': 472, 'PB73': 671, 'PB93': 24203, 'GGACTCCT-GTAAGGAG': 559, 'PB39': 359, 'TAAGGCGA-AGAGTAGA': 407, 'GTAGAGGA-CTAAGCCT': 333, 'PB32': 434, 'PB33': 381, 'PB79': 436, 'PB_9': 473, 'PB64': 442, 'PB66': 423, 'CAGAGAGG-ACTGCATA': 344, 'GTAGAGGA-CTCTCTAT': 349, 'CAGAGAGG-TATCCTCT': 365, 'CTCTCTAC-TATCCTCT': 274, 'AAGAGGCA-CTAAGCCT': 346, 'AGGCAGAA-AGAGTAGA': 396, 'AAGAGGCA-TATCCTCT': 276, 'CAGAGAGG-AAGGAGTA': 185, 'PB53': 263, 'CTCTCTAC-TAGATCGC': 363, 'GCTACGCT-ACTGCATA': 332, 'CAGAGAGG-TAGATCGC': 235, 'GGACTCCT-TAGATCGC': 420, 'CTCTCTAC-ACTGCATA': 223, 'CTCTCTAC-AGAGTAGA': 270, 'PB37': 241, 'PB16': 409, 'GTAGAGGA-TATCCTCT': 359, 'CAGAGAGG-AGAGTAGA': 174, 'AAGAGGCA-GTAAGGAG': 297, 'CAGAGAGG-GTAAGGAG': 187, 'PB77': 555, 'PB29': 256, 'CTCTCTAC-CTCTCTAT': 473, 'PB83': 187, 'TAAGGCGA-AAGGAGTA': 329, 'PB_5': 62, 'PB43': 341, 'TAAGGCGA-ACTGCATA': 54, 'PB40': 161, 'PB45': 312, 'CAGAGAGG-CTAAGCCT': 305, 'AGGCAGAA-CTCTCTAT': 132, 'PB34': 371, 'CTCTCTAC-CTAAGCCT': 57, 'PB85': 734, 'AAGAGGCA-TAGATCGC': 153, 'TCCTGAGC-CTAAGCCT': 314, 'AAGAGGCA-CTCTCTAT': 104, 'PB87': 87, 'CGTACTAG-AGAGTAGA': 219, 'PB25': 424, 'PB10': 307, 'TAAGGCGA-TATCCTCT': 195, 'PB57': 255, 'PB_4': 239, 'TAAGGCGA-CTCTCTAT': 44, 'GTAGAGGA-AGAGTAGA': 242, 'PB88': 109, 'GGACTCCT-ACTGCATA': 474, 'PB18': 307, 'CGTACTAG-AAGGAGTA': 198, 'AGGCAGAA-GTAAGGAG': 134, 'TCCTGAGC-TATCCTCT': 82, 'PB21': 216, 'CGTACTAG-ACTGCATA': 124, 'PB41': 458, 'CAGAGAGG-CTCTCTAT': 301, 'PB76': 247, 'AGGCAGAA-TAGATCGC': 107, 'PB_3': 268, 'TCCTGAGC-AGAGTAGA': 138, 'PB11': 218, 'PB27': 139, 'TAAGGCGA-CTAAGCCT': 62, 'PB44': 223, 'PB84': 115, 'AGGCAGAA-AAGGAGTA': 43, 'PB90': 140, 'CGTACTAG-CTAAGCCT': 62, 'GTAGAGGA-TAGATCGC': 35, 'PB35': 178, 'CGTACTAG-TATCCTCT': 19, 'TCCTGAGC-GTAAGGAG': 495, 'PB82': 51, 'CGTACTAG-TAGATCGC': 78, 'PB50': 37, 'GGACTCCT-CTAAGCCT': 118, 'CGTACTAG-GTAAGGAG': 34, 'AGGCAGAA-TATCCTCT': 31, 'CGTACTAG-CTCTCTAT': 56, 'TCCTGAGC-ACTGCATA': 43, 'PB58': 51, 'PB42': 41}


# Next, read in the file with the HGT blocks for each strain (disctionary-like format)

hgt_g_blocks_file_path55 = '/home/ada/Desktop/Shraiman_lab/srb_data/hgt_blocks_mixing_layer_green_all_strains_5_5.txt'
hgt_g_blocks_file_path53 = '/home/ada/Desktop/Shraiman_lab/srb_data/hgt_blocks_mixing_layer_green_all_strains_5_3.txt'
hgt_g_blocks_file_path105 = '/home/ada/Desktop/Shraiman_lab/srb_data/hgt_blocks_mixing_layer_green_all_strains_10_5.txt'
hgt_g_blocks_file_path57 = '/home/ada/Desktop/Shraiman_lab/srb_data/hgt_blocks_mixing_layer_green_all_strains_5_7.txt'

hgt_green_blocks_files = [hgt_g_blocks_file_path55, hgt_g_blocks_file_path53, hgt_g_blocks_file_path105,hgt_g_blocks_file_path57]

hgt_y_blocks_file_path55 = '/home/ada/Desktop/Shraiman_lab/srb_data/hgt_blocks_mixing_layer_yellow_all_strains_5_5.txt'
hgt_y_blocks_file_path53 = '/home/ada/Desktop/Shraiman_lab/srb_data/hgt_blocks_mixing_layer_yellow_all_strains_5_3.txt'
hgt_y_blocks_file_path105 = '/home/ada/Desktop/Shraiman_lab/srb_data/hgt_blocks_mixing_layer_yellow_all_strains_10_5.txt'
hgt_y_blocks_file_path57 = '/home/ada/Desktop/Shraiman_lab/srb_data/hgt_blocks_mixing_layer_yellow_all_strains_5_7.txt'

hgt_yellow_blocks_files = [hgt_y_blocks_file_path55, hgt_y_blocks_file_path53, hgt_y_blocks_file_path105, hgt_y_blocks_file_path57]

fig, axs = plt.subplots(2, 2)
for i in range(0, len(hgt_green_blocks_files)):

    hgt_blocks_green, hgt_blocks_yellow = {}, {}
    hgt_lengths_green = {}

    with open(hgt_green_blocks_files[i], 'r') as f:
        for line in f:
            line = line.split(':')
            strain, blocks = line[0], [int(i) for i in line[1].split(',')]
            blocks = [i for i in blocks if i != 0]   # clean zeros
            hgt_blocks_green[strain] = blocks
            hgt_lengths_green[strain] = np.sum(np.array(blocks))



    with open(hgt_yellow_blocks_files[i], 'r') as f:
        for line in f:
            line = line.split(':')
            strain, blocks = line[0], [int(i) for i in line[1].split(',')]
            blocks = [i for i in blocks if i != 0]   # clean zeros
            hgt_blocks_yellow[strain] = blocks


    # Plot the number of HGT blocks against the number of singleton SNPs for each mixing layer strain
    # Differentiate between green and yellow blocks

    
    # Combine green and yellow blocks
    x_g, y_g = [], []
    for strain in strains:
        total = 0
        if strain not in f_clade_strains and strain not in basal_AB_strains:
            total += len(hgt_blocks_green[strain])
            total += len(hgt_blocks_yellow[strain])

            if total < 70 and singleton_counts[strain] < 700:
                x_g.append(total)
                y_g.append(singleton_counts[strain])
            
    
    a, b = np.polyfit(y_g, x_g, 1)

    l = hgt_green_blocks_files[i].split('_')[-1].split('.')[0]
    penalty = hgt_green_blocks_files[i].split('_')[-2]
    label = 'penalty: ' + penalty + ' L: ' + l
    #plt.scatter(y_g, x_g, c='green')
    #plt.scatter(y_y, x_y, c='orange')
    y_g = np.array(y_g)
    #add line of best fit to plot
    plt.plot(y_g, a*y_g+b, label = label)

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




plt.suptitle('Number of HGT blocks vs number of singletons in the mixing layer strains (SRB): \n varying penalty variables + best fit line')
plt.xlabel('Singleton number')
plt.ylabel('HGT blocks number')
plt.legend()
plt.show()


# Check the pearson's correlation for singletons vs HGT block number
# # Use only one file 

hgt_blocks_green, hgt_blocks_yellow = {}, {}
hgt_lengths_green = {}

with open(hgt_g_blocks_file_path55, 'r') as f:
    for line in f:
        line = line.split(':')
        strain, blocks = line[0], [int(i) for i in line[1].split(',')]
        blocks = [i for i in blocks if i != 0]   # clean zeros
        hgt_blocks_green[strain] = blocks
        hgt_lengths_green[strain] = np.sum(np.array(blocks))



with open(hgt_y_blocks_file_path55, 'r') as f:
    for line in f:
        line = line.split(':')
        strain, blocks = line[0], [int(i) for i in line[1].split(',')]
        blocks = [i for i in blocks if i != 0]   # clean zeros
        hgt_blocks_yellow[strain] = blocks


# Plot the number of HGT blocks against the number of singleton SNPs for each mixing layer strain
# Differentiate between green and yellow blocks

# Combine green and yellow blocks
x_g, y_g = [], []
for strain in strains:
    total = 0
    if strain not in f_clade_strains and strain not in basal_AB_strains:
        total += len(hgt_blocks_green[strain])
        total += len(hgt_blocks_yellow[strain])

        if total < 70 and singleton_counts[strain] < 700:
            x_g.append(total)
            y_g.append(singleton_counts[strain])
        
y_g = np.array(y_g)
a, b = np.polyfit(y_g, x_g, 1)
plt.plot(y_g, a*y_g+b)
plt.scatter(y_g, x_g) 
plt.suptitle('Number of HGT blocks vs number of singletons in the mixing layer strains (SRB): \n varying penalty variables + best fit line')
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
