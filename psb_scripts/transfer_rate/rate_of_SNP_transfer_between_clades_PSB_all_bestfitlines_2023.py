
###############################################################
# What is the frequency of fixed SNPs transfer from one clade to another?
# Given all the transfer rate data for all clades, let's plot best fit lines
# and compare them for different clades (origins of transfer)
###############################################################

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit


def find_fixed_loci(clade, data_index, data_positions):

	# Create a mask that will only show the defined working clade values
	mask = [True if x in clade else False for x in data_index]
	fixed_loci = []

	# Check SNP frequency in the current clade.
	for i in data_positions:
		p = np.array(data[i].tolist())[mask]
		if sum(p==1) >= len(clade)/2:  				
			fixed_loci.append(i)

	return(fixed_loci)

def find_avg_counts_and_singletons(clade, other_clades):
	
    # Find all fixed topAB loci
    fixed = find_fixed_loci(clade, data_index, data_positions)

    # Calculate the average counts of the fixed topAB loci outside of the topAB clade
    # AND plot it against the avg number of singletons for each clade.
    # Whole tree
    mask = [True if x in fixed else False for x in data_positions]

    # Order: ladder, BC clade, E clade, basal clade
    avg_counts, avg_singletons = [], []

    for i in other_clades:
        counts, singletons = [], []
        for s in i:
            vector = np.array(data.loc[s][mask])
            counts.append(sum(vector==1))
            singletons.append(psb_singleton_counts[s])

        # Avg the counts and singletons

        counts, singletons = np.array(counts), np.array(singletons)
        avg_counts.append(np.mean(counts)/float(len(fixed)))  # Normalize the transferred counts to the total number of fixed loci
        avg_singletons.append(np.mean(singletons))

    return(avg_counts, avg_singletons)

## All plotting functions

def plot_bestfitline(avg_singletons, avg_counts, color, label):
    
    a, b = np.polyfit(avg_singletons, avg_counts, 1)
    avg_singletons = np.array(avg_singletons)
    plt.plot(avg_singletons, a*avg_singletons+b, color=color, label=label)

    return()

def func(x, a):
    return a * x


def plot_bestfitline_origin(avg_singletons, avg_counts, color, label):

    avg_singletons = [0] + list(avg_singletons)
    avg_counts = [0] + avg_counts

    popt, pcov = curve_fit(func, avg_singletons, avg_counts)
    plt.plot(avg_singletons, func(avg_singletons, popt),color=color, label=label)
    
    return()

def plot_bestfitline_exponential(avg_singletons, avg_counts, color, label):
     
    # Fit a weighted polynomial of degree 1 (a linear function) to the data
    p = np.polyfit(avg_singletons, np.log(avg_counts), 1,w=np.sqrt(avg_counts))

    # Convert the polynomial back into an exponential
    a = np.exp(p[1])
    b = p[0]
    x_fitted = np.linspace(np.min(avg_singletons), np.max(avg_singletons), 100)
    y_fitted = a * np.exp(b * x_fitted)

    plt.plot(x_fitted, y_fitted, color=color, label=label)
    
    return()


# Read in the SNP data
data = pd.read_csv('/home/ada/Desktop/Shraiman_lab/data/snp_data_2020.csv', index_col=0)

data_index = data.index.values.tolist()  # names of bacterial samples
data_positions = data.columns.tolist()

clade_9_strains = ['PB87','PB40','TCCTGAGC-TAGATCGC','TCCTGAGC-CTCTCTAT','PB73','AAGAGGCA-TATCCTCT','GGACTCCT-AGAGTAGA']
ladder = ['GGACTCCT-CTAAGCCT','PB80','PB63','AAGAGGCA-TAGATCGC','PB39','PB_8','PB24','PB64','PB55','AAGAGGCA-AGAGTAGA','GCTACGCT-GTAAGGAG','AAGAGGCA-ACTGCATA','AAGAGGCA-AAGGAGTA','GCTACGCT-ACTGCATA','PB31','AAGAGGCA-GTAAGGAG','AAGAGGCA-CTAAGCCT','AAGAGGCA-CTCTCTAT','PB32','PB16','PB47','PB48']
bc_strains = ['TAAGGCGA-TATCCTCT','PB59','PB45','PB77','GGACTCCT-TAGATCGC','GGACTCCT-ACTGCATA','PB78','PB61','PB13','PB85','PB37','PB_5','AGGCAGAA-ACTGCATA','AGGCAGAA-CTAAGCCT','PB25','PB67','PB76','PB59','PB84','PB52','PB44']
e_clade = ['PB11','CTCTCTAC-AAGGAGTA','GTAGAGGA-CTAAGCCT','CTCTCTAC-CTAAGCCT','CAGAGAGG-CTCTCTAT','CTCTCTAC-CTCTCTAT','CTCTCTAC-TAGATCGC','CTCTCTAC-ACTGCATA','CAGAGAGG-TAGATCGC','CAGAGAGG-CTAAGCCT','CAGAGAGG-GTAAGGAG','CTCTCTAC-TATCCTCT','CAGAGAGG-ACTGCATA','CAGAGAGG-AGAGTAGA','CAGAGAGG-TATCCTCT','CAGAGAGG-AAGGAGTA','CTCTCTAC-AGAGTAGA','CTCTCTAC-GTAAGGAG','PB34','PB26']
basal_clade = ['PB29','CGTACTAG-TATCCTCT','PB57','PB51','AGGCAGAA-GTAAGGAG','PB90','PB36','PB82','PB20','PB65','PB60','PB21','PB58','PB43','PB66','GTAGAGGA-AAGGAGTA','PB12','AGGCAGAA-AAGGAGTA','GTAGAGGA-CTCTCTAT','PB75','AGGCAGAA-AGAGTAGA','PB41','CGTACTAG-TAGATCGC','PB25','PB_2']
strains = ['PB87','PB40','TCCTGAGC-TAGATCGC','TCCTGAGC-CTCTCTAT','PB73','AAGAGGCA-TATCCTCT','GGACTCCT-AGAGTAGA','GGACTCCT-CTAAGCCT','GGACTCCT-CTAAGCCT','PB80','PB63','AAGAGGCA-TAGATCGC','PB39','PB_8','PB24','PB64','PB55','AAGAGGCA-AGAGTAGA','GCTACGCT-GTAAGGAG','AAGAGGCA-ACTGCATA','AAGAGGCA-AAGGAGTA','GCTACGCT-ACTGCATA','PB31','AAGAGGCA-GTAAGGAG','AAGAGGCA-CTAAGCCT','AAGAGGCA-CTCTCTAT','PB32','PB16','PB47','PB48','TAAGGCGA-TATCCTCT','PB59','PB45','PB77','GGACTCCT-TAGATCGC','GGACTCCT-ACTGCATA','PB78','PB61','PB13','PB85','PB37','PB_5','AGGCAGAA-ACTGCATA','AGGCAGAA-CTAAGCCT','PB25','PB67','PB76','PB59','PB84','PB52','PB44','PB11','CTCTCTAC-AAGGAGTA','GTAGAGGA-CTAAGCCT','CTCTCTAC-CTAAGCCT','CAGAGAGG-CTCTCTAT','CTCTCTAC-CTCTCTAT','CTCTCTAC-TAGATCGC','CTCTCTAC-ACTGCATA','CAGAGAGG-TAGATCGC','CAGAGAGG-CTAAGCCT','CAGAGAGG-GTAAGGAG','CTCTCTAC-TATCCTCT','CAGAGAGG-ACTGCATA','CAGAGAGG-AGAGTAGA','CAGAGAGG-TATCCTCT','CAGAGAGG-AAGGAGTA','CTCTCTAC-AGAGTAGA','CTCTCTAC-GTAAGGAG','PB34','PB26','PB29','CGTACTAG-TATCCTCT','PB57','PB51','AGGCAGAA-GTAAGGAG','PB90','PB36','PB82','PB20','PB65','PB60','PB21','PB58','PB43','PB66','GTAGAGGA-AAGGAGTA','PB12','AGGCAGAA-AAGGAGTA','GTAGAGGA-CTCTCTAT','PB75','AGGCAGAA-AGAGTAGA','PB41','CGTACTAG-TAGATCGC','PB25','PB_2']

psb_singleton_counts = {'PB63': 2826,'PB28': 1369, 'GGACTCCT-AGAGTAGA': 699, 'AAGAGGCA-AGAGTAGA': 966, 
                    'PB73': 767, 'PB76': 1250, 'AAGAGGCA-ACTGCATA': 1102, 'GGACTCCT-GTAAGGAG': 711, 
                    'GCTACGCT-GTAAGGAG': 931, 'TCCTGAGC-CTCTCTAT': 597, 'PB55': 855, 'PB_5': 653, 
                    'PB31': 997, 'AAGAGGCA-TATCCTCT': 809, 'PB24': 655, 'PB87': 398, 
                    'AAGAGGCA-AAGGAGTA': 813, 'PB84': 782, 'PB40': 366, 'AGGCAGAA-ACTGCATA': 951, 
                    'AGGCAGAA-CTAAGCCT': 886, 'PB67': 772, 'PB_8': 674, 'PB80': 680, 'PB39': 470, 
                    'PB47': 576, 'PB59': 607, 'PB37': 116, 'PB61': 289, 'TCCTGAGC-TAGATCGC': 504, 
                    'PB77': 150, 'PB69': 320, 'TAAGGCGA-AAGGAGTA': 587, 'PB13': 98, 'PB82': 336, 
                    'PB48': 432, 'PB52': 497, 'GGACTCCT-TAGATCGC': 96, 'PB11': 189, 
                    'GGACTCCT-CTAAGCCT': 551, 'PB44': 351, 'PB36': 397, 'AAGAGGCA-GTAAGGAG': 358, 
                    'PB16': 344, 'PB32': 271, 'PB20': 492, 'GCTACGCT-ACTGCATA': 337, 'PB85': 131, 
                    'PB65': 444, 'PB58': 297, 'PB64': 304, 'CTCTCTAC-GTAAGGAG': 266, 'PB60': 279, 
                    'TAAGGCGA-GTAAGGAG': 315, 'AGGCAGAA-GTAAGGAG': 101, 'PB21': 132, 
                    'TAAGGCGA-TATCCTCT': 194, 'AAGAGGCA-CTCTCTAT': 188, 'PB78': 129, 'PB_1': 86, 
                    'PB91': 120, 'AAGAGGCA-TAGATCGC': 180, 'CTCTCTAC-AAGGAGTA': 100, 
                    'GTAGAGGA-ACTGCATA': 96, 'CAGAGAGG-AAGGAGTA': 67, 'CTCTCTAC-AGAGTAGA': 74, 
                    'GTAGAGGA-GTAAGGAG': 96, 'PB43': 93, 'PB34': 179, 'AAGAGGCA-CTAAGCCT': 182, 
                    'PB92': 125, 'CTCTCTAC-TATCCTCT': 76, 'AGGCAGAA-AAGGAGTA': 65, 'PB90': 101, 
                    'PB12': 98, 'CAGAGAGG-TATCCTCT': 239, 'PB18': 87, 'PB26': 139, 
                    'GTAGAGGA-AAGGAGTA': 61, 'PB41': 37, 'AGGCAGAA-AGAGTAGA': 179, 'PB75': 104, 
                    'PB45': 117, 'CAGAGAGG-TAGATCGC': 85, 'TAAGGCGA-AGAGTAGA': 67, 
                    'CAGAGAGG-GTAAGGAG': 81, 'PB17': 55, 'CAGAGAGG-ACTGCATA': 96, 
                    'CTCTCTAC-CTCTCTAT': 62, 'PB51': 60, 'CGTACTAG-AGAGTAGA': 55, 'PB66': 97, 
                    'PB25': 59, 'PB83': 75, 'CGTACTAG-CTAAGCCT': 54, 'CTCTCTAC-TAGATCGC': 77, 
                    'TCCTGAGC-AAGGAGTA': 218, 'CAGAGAGG-AGAGTAGA': 107, 'PB_3': 66, 
                    'AGGCAGAA-CTCTCTAT': 80, 'CAGAGAGG-CTAAGCCT': 63, 'TAAGGCGA-TAGATCGC': 180, 
                    'PB81': 42, 'CGTACTAG-ACTGCATA': 36, 'CGTACTAG-AAGGAGTA': 61, 'PB89': 36, 
                    'TAAGGCGA-ACTGCATA': 17, 'PB74': 26, 'CTCTCTAC-ACTGCATA': 79, 'PB57': 40, 
                    'PB49': 39, 'GTAGAGGA-CTAAGCCT': 89, 'CGTACTAG-CTCTCTAT': 39, 
                    'GTAGAGGA-CTCTCTAT': 54, 'PB_9': 36, 'PB_4': 26, 'GTAGAGGA-AGAGTAGA': 36, 
                    'GTAGAGGA-TATCCTCT': 44, 'PB29': 29, 'PB10': 42, 'PB_2': 29, 'AGGCAGAA-TAGATCGC': 35,
                    'CGTACTAG-TATCCTCT': 32, 'PB53': 32, 'GGACTCCT-ACTGCATA': 410, 'PB19': 22, 
                    'CTCTCTAC-CTAAGCCT': 67, 'TCCTGAGC-AGAGTAGA': 30, 'CAGAGAGG-CTCTCTAT': 42, 
                    'PB35': 24, 'PB27': 18, 'GTAGAGGA-TAGATCGC': 26, 'CGTACTAG-TAGATCGC': 24, 
                    'TCCTGAGC-CTAAGCCT': 13, 'TAAGGCGA-CTAAGCCT': 30, 'TCCTGAGC-TATCCTCT': 19, 
                    'PB42': 24, 'TCCTGAGC-GTAAGGAG': 17, 'TAAGGCGA-CTCTCTAT': 18, 'TCCTGAGC-ACTGCATA': 23, 
                    'AGGCAGAA-TATCCTCT': 22, 'CGTACTAG-GTAAGGAG': 21}

###
# Data for 8 clade
###
other_clades = [ladder, bc_strains, e_clade, basal_clade]
avg_counts_8, avg_singletons_8 = find_avg_counts_and_singletons(clade_9_strains, other_clades)

###
# Data for ladder clade
###
other_clades = [clade_9_strains, bc_strains, e_clade, basal_clade]
avg_counts_ladder, avg_singletons_ladder = find_avg_counts_and_singletons(ladder, other_clades)

###
# Data for BC clade
###
other_clades = [clade_9_strains, ladder, e_clade, basal_clade]
avg_counts_bc, avg_singletons_bc = find_avg_counts_and_singletons(bc_strains, other_clades)

###
# Data for E clade
###
other_clades = [clade_9_strains, ladder, bc_strains, basal_clade]
avg_counts_e, avg_singletons_e = find_avg_counts_and_singletons(e_clade, other_clades)

###
# Plot the scatters on the same plot with different colors for different clade transfers.
###

plt.scatter(avg_singletons_8, avg_counts_8, color='r', label='8 clade')
plt.scatter(avg_singletons_ladder, avg_counts_ladder, color='blue', label='ladder')
plt.scatter(avg_singletons_bc, avg_counts_bc, color='g', label='BC clade')
plt.scatter(avg_singletons_e, avg_counts_e, color='m', label='E clade')

plt.title('What is the rate of transfer of fixed SNPs from one clade to other clades?\n all clades scatters')
plt.xlabel('average number of singletons per clade')
plt.ylabel('average number of SNPs transferred from different clades (normalized)')
plt.legend()
plt.show()

###
# Make a best fit line for the graph. For this don't fit through (0,0) origin
###
plot_bestfitline(avg_singletons_8, avg_counts_8, 'r', '8 clade')
plot_bestfitline(avg_singletons_ladder, avg_counts_ladder, 'blue', 'ladder')
plot_bestfitline(avg_singletons_bc, avg_counts_bc, 'g', 'BC clade')
plot_bestfitline(avg_singletons_e, avg_counts_e, 'm', 'E clade')

plt.title('What is the rate of transfer of fixed SNPs from one clade to other clades?\n all clades best fit lines')
plt.xlabel('average number of singletons per clade')
plt.ylabel('average number of SNPs transferred from different clades')
plt.legend()
plt.show()

###
# Make a best fit line for the graph. For this fit through (0,0) origin
###
plot_bestfitline_origin(avg_singletons_8, avg_counts_8, 'r', '8 clade')
plot_bestfitline_origin(avg_singletons_ladder, avg_counts_ladder, 'blue', 'ladder')
plot_bestfitline_origin(avg_singletons_bc, avg_counts_bc, 'g', 'BC clade')
plot_bestfitline_origin(avg_singletons_e, avg_counts_e, 'm', 'E clade')

plt.title('What is the rate of transfer of fixed SNPs from one clade to other clades?\n all clades best fit lines through origin')
plt.xlabel('average number of singletons per clade')
plt.ylabel('average number of SNPs transferred from different clades')
plt.legend()
plt.show()

###
# Make a best fit line for the graph. This time fit an exponential.
# We can convert an exponential function into a polynomial one
# The act of transforming a polynomial function into an exponential one has the effect
#   of increasing large values much more than it does small values, and thus it has the 
#   effect of increasing the distance to the fitted curve for large values more than it does 
#   for small values. This can be mitigated by adding a ‘weight’ proportional to y:
#   tell polyfit() to lend more importance to data points with a large y-value
#   from: https://rowannicholls.github.io/python/curve_fitting/exponential.html 
###
plot_bestfitline_exponential(avg_singletons_8, avg_counts_8, 'r', '8 clade')
plot_bestfitline_exponential(avg_singletons_ladder, avg_counts_ladder, 'blue', 'ladder')
plot_bestfitline_exponential(avg_singletons_bc, avg_counts_bc, 'g', 'BC clade')
plot_bestfitline_exponential(avg_singletons_e, avg_counts_e, 'm', 'E clade')

plt.title('What is the rate of transfer of fixed SNPs from one clade to other clades?\n all clades best fit lines weighted exponential')
plt.xlabel('average number of singletons per clade')
plt.ylabel('average number of SNPs transferred from different clades')
plt.legend()
plt.show()







