
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
            vector = np.array(data.loc[s][1:][mask])
            counts.append(sum(vector==1))
            singletons.append(srb_singleton_counts[s])

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
data = pd.read_csv('/home/ada/Desktop/Shraiman_lab/srb_data/srb_snp_data_2021_chr.csv', index_col=0)

data_index = data.index.values.tolist()  # names of bacterial samples
data_positions = data.columns.tolist()[1:]

basal_ab_clade = ['CGTACTAG-TATCCTCT','CGTACTAG-GTAAGGAG','PB50','AGGCAGAA-TATCCTCT','GTAGAGGA-TAGATCGC','TAAGGCGA-CTAAGCCT','PB82','PB90','PB27','AGGCAGAA-GTAAGGAG','CGTACTAG-CTAAGCCT','TAAGGCGA-TATCCTCT','PB18','CGTACTAG-AAGGAGTA','GTAGAGGA-AGAGTAGA','PB10','PB34','PB76']
f_clade = ['CAGAGAGG-AAGGAGTA','CAGAGAGG-GTAAGGAG','CTCTCTAC-TATCCTCT','AAGAGGCA-CTCTCTAT','AAGAGGCA-TATCCTCT','AAGAGGCA-CTAAGCCT','GCTACGCT-ACTGCATA','PB39','PB32','PB40','PB88','PB16','PB79','PB87','AAGAGGCA-TAGATCGC','AAGAGGCA-GTAAGGAG','PB64']
topAB_strains = ['TAAGGCGA-CTCTCTAT','TAAGGCGA-ACTGCATA','PB83','AGGCAGAA-CTCTCTAT','AGGCAGAA-TAGATCGC','PB53','GTAGAGGA-CTCTCTAT','TCCTGAGC-ACTGCATA','PB58','TCCTGAGC-AGAGTAGA','PB33','PB_4','PB84','TCCTGAGC-GTAAGGAG']
basal_clade = ['CGTACTAG-ACTGCATA','TCCTGAGC-TATCCTCT','CGTACTAG-AGAGTAGA','CGTACTAG-TAGATCGC','AGGCAGAA-AGAGTAGA','PB66','TAAGGCGA-AAGGAGTA','PB43','CGTACTAG-CTCTCTAT','PB_3','PB44','AGGCAGAA-AAGGAGTA','PB42']
mixing_layer = ['CAGAGAGG-TAGATCGC','CAGAGAGG-AGAGTAGA','CAGAGAGG-ACTGCATA','CTCTCTAC-CTAAGCCT','CTCTCTAC-ACTGCATA','PB37','PB_5','CAGAGAGG-TATCCTCT','CTCTCTAC-AGAGTAGA','CAGAGAGG-CTCTCTAT','CTCTCTAC-CTCTCTAT','PB11','CTCTCTAC-TAGATCGC','GGACTCCT-TAGATCGC','GGACTCCT-ACTGCATA','PB77','PB29','PB35','PB45','TCCTGAGC-CTAAGCCT','PB85','CAGAGAGG-CTAAGCCT','GGACTCCT-CTAAGCCT','PB57','PB21']
full_tree = ['CAGAGAGG-AAGGAGTA','CAGAGAGG-GTAAGGAG','CTCTCTAC-TATCCTCT','AAGAGGCA-CTCTCTAT','AAGAGGCA-TATCCTCT','AAGAGGCA-CTAAGCCT','GCTACGCT-ACTGCATA','PB39','PB32','PB40','PB88','PB16','PB79','PB87','AAGAGGCA-TAGATCGC','AAGAGGCA-GTAAGGAG','PB64','CAGAGAGG-TAGATCGC','CAGAGAGG-AGAGTAGA','CAGAGAGG-ACTGCATA','CTCTCTAC-CTAAGCCT','CTCTCTAC-ACTGCATA','PB37','PB_5','CAGAGAGG-TATCCTCT','CTCTCTAC-AGAGTAGA','CAGAGAGG-CTCTCTAT','CTCTCTAC-CTCTCTAT','PB11','CTCTCTAC-TAGATCGC','GGACTCCT-TAGATCGC','GGACTCCT-ACTGCATA','PB77','PB29','PB35','PB45','TCCTGAGC-CTAAGCCT','PB85','CAGAGAGG-CTAAGCCT','GGACTCCT-CTAAGCCT','PB57','PB21','TAAGGCGA-CTCTCTAT','TAAGGCGA-ACTGCATA','PB83','AGGCAGAA-CTCTCTAT','AGGCAGAA-TAGATCGC','PB53','GTAGAGGA-CTCTCTAT','TCCTGAGC-ACTGCATA','PB58','TCCTGAGC-AGAGTAGA','PB33','PB_4','PB84','TCCTGAGC-GTAAGGAG','PB41','PB25','GTAGAGGA-TATCCTCT','CGTACTAG-TATCCTCT','CGTACTAG-GTAAGGAG','PB50','AGGCAGAA-TATCCTCT','GTAGAGGA-TAGATCGC','TAAGGCGA-CTAAGCCT','PB82','PB90','PB27','AGGCAGAA-GTAAGGAG','CGTACTAG-CTAAGCCT','TAAGGCGA-TATCCTCT','PB18','CGTACTAG-AAGGAGTA','GTAGAGGA-AGAGTAGA','PB10','PB34','PB76','CGTACTAG-ACTGCATA','TCCTGAGC-TATCCTCT','CGTACTAG-AGAGTAGA','CGTACTAG-TAGATCGC','AGGCAGAA-AGAGTAGA','PB66','TAAGGCGA-AAGGAGTA','PB43','CGTACTAG-CTCTCTAT','PB_3','PB44','AGGCAGAA-AAGGAGTA','PB42']

srb_singleton_counts = {'AAGAGGCA-AAGGAGTA': 483, 'PB36': 481, 'PB80': 646, 'PB23': 607, 'PB78': 547, 
                    'TCCTGAGC-AAGGAGTA': 413, 'PB47': 454, 'PB48': 481, 'AAGAGGCA-ACTGCATA': 423, 
                    'PB60': 549, 'PB_1': 554, 'PB91': 427, 'PB95': 523, 'PB19': 423, 'PB59': 584, 
                    'GTAGAGGA-ACTGCATA': 465, 'PB81': 481, 'PB28': 545, 'PB68': 501, 
                    'GTAGAGGA-GTAAGGAG': 472, 'PB73': 671, 'PB93': 24203, 'GGACTCCT-GTAAGGAG': 559,
                    'PB39': 359, 'TAAGGCGA-AGAGTAGA': 407, 'GTAGAGGA-CTAAGCCT': 333, 'PB32': 434, 
                    'PB33': 381, 'PB79': 436, 'PB_9': 473, 'PB64': 442, 'PB66': 423, 
                    'CAGAGAGG-ACTGCATA': 344, 'GTAGAGGA-CTCTCTAT': 349, 'CAGAGAGG-TATCCTCT': 365,
                    'CTCTCTAC-TATCCTCT': 274, 'AAGAGGCA-CTAAGCCT': 346, 'AGGCAGAA-AGAGTAGA': 396,
                    'AAGAGGCA-TATCCTCT': 276, 'CAGAGAGG-AAGGAGTA': 185, 'PB53': 263, 
                    'CTCTCTAC-TAGATCGC': 363, 'GCTACGCT-ACTGCATA': 332, 'CAGAGAGG-TAGATCGC': 235,
                    'GGACTCCT-TAGATCGC': 420, 'CTCTCTAC-ACTGCATA': 223, 'CTCTCTAC-AGAGTAGA': 270, 
                    'PB37': 241, 'PB16': 409, 'GTAGAGGA-TATCCTCT': 359, 'CAGAGAGG-AGAGTAGA': 174,
                    'AAGAGGCA-GTAAGGAG': 297, 'CAGAGAGG-GTAAGGAG': 187, 'PB77': 555, 'PB29': 256, 
                    'CTCTCTAC-CTCTCTAT': 473, 'PB83': 187, 'TAAGGCGA-AAGGAGTA': 329, 'PB_5': 62, 
                    'PB43': 341, 'TAAGGCGA-ACTGCATA': 54, 'PB40': 161, 'PB45': 312, 
                    'CAGAGAGG-CTAAGCCT': 305, 'AGGCAGAA-CTCTCTAT': 132, 'PB34': 371, 
                    'CTCTCTAC-CTAAGCCT': 57, 'PB85': 734, 'AAGAGGCA-TAGATCGC': 153, 
                    'TCCTGAGC-CTAAGCCT': 314, 'AAGAGGCA-CTCTCTAT': 104, 'PB87': 87, 
                    'CGTACTAG-AGAGTAGA': 219, 'PB25': 424, 'PB10': 307, 'TAAGGCGA-TATCCTCT': 195, 
                    'PB57': 255, 'PB_4': 239, 'TAAGGCGA-CTCTCTAT': 44, 'GTAGAGGA-AGAGTAGA': 242, 
                    'PB88': 109, 'GGACTCCT-ACTGCATA': 474, 'PB18': 307, 'CGTACTAG-AAGGAGTA': 198, 
                    'AGGCAGAA-GTAAGGAG': 134, 'TCCTGAGC-TATCCTCT': 82, 'PB21': 216, 
                    'CGTACTAG-ACTGCATA': 124, 'PB41': 458, 'CAGAGAGG-CTCTCTAT': 301, 'PB76': 247, 
                    'AGGCAGAA-TAGATCGC': 107, 'PB_3': 268, 'TCCTGAGC-AGAGTAGA': 138, 'PB11': 218, 
                    'PB27': 139, 'TAAGGCGA-CTAAGCCT': 62, 'PB44': 223, 'PB84': 115, 
                    'AGGCAGAA-AAGGAGTA': 43, 'PB90': 140, 'CGTACTAG-CTAAGCCT': 62, 
                    'GTAGAGGA-TAGATCGC': 35, 'PB35': 178, 'CGTACTAG-TATCCTCT': 19, 
                    'TCCTGAGC-GTAAGGAG': 495, 'PB82': 51, 'CGTACTAG-TAGATCGC': 78, 'PB50': 37, 
                    'GGACTCCT-CTAAGCCT': 118, 'CGTACTAG-GTAAGGAG': 34, 'AGGCAGAA-TATCCTCT': 31, 
                    'CGTACTAG-CTCTCTAT': 56, 'TCCTGAGC-ACTGCATA': 43, 'PB58': 51, 'PB42': 41}


###
# Data for basal AB clade
###
other_clades = [mixing_layer, f_clade, topAB_strains, basal_clade]
avg_counts_basal_ab, avg_singletons_basal_ab = find_avg_counts_and_singletons(basal_ab_clade, other_clades)

###
# Data for F clade
###
other_clades = [mixing_layer, basal_ab_clade, topAB_strains, basal_clade]
avg_counts_f, avg_singletons_f = find_avg_counts_and_singletons(f_clade, other_clades)

###
# Data for top AB clade
###
other_clades = [mixing_layer, basal_ab_clade, f_clade, basal_clade]
avg_counts_topab, avg_singletons_topab = find_avg_counts_and_singletons(topAB_strains, other_clades)

###
# Data for mixing clade
###
other_clades = [f_clade, basal_ab_clade, topAB_strains, basal_clade]
avg_counts_mixing, avg_singletons_mixing = find_avg_counts_and_singletons(mixing_layer, other_clades)

###
# Plot the scatters on the same plot with different colors for different clade transfers.
###

plt.scatter(avg_singletons_f, avg_counts_f, color='blue', label='F clade')
plt.scatter(avg_singletons_basal_ab, avg_counts_basal_ab, color='r', label='basal AB clade')
plt.scatter(avg_singletons_topab, avg_counts_topab, color='g', label='top AB clade')
plt.scatter(avg_singletons_mixing, avg_counts_mixing, color='m', label='mixing layer')

plt.title('What is the rate of transfer of fixed SNPs from one clade to other clades?\n all clades scatters')
plt.xlabel('average number of singletons per clade')
plt.ylabel('average number of SNPs transferred (normalized)')
plt.legend()
plt.show()

###
# Make a best fit line for the graph. For this don't fit through (0,0) origin
###
plot_bestfitline(avg_singletons_basal_ab, avg_counts_basal_ab, 'r', 'basal AB clade')
plot_bestfitline(avg_singletons_f, avg_counts_f, 'blue', 'F clade')
plot_bestfitline(avg_singletons_topab, avg_counts_topab, 'g', 'top AB clade')
plot_bestfitline(avg_singletons_mixing, avg_counts_mixing, 'm', 'mixing layer')

plt.title('What is the rate of transfer of fixed SNPs from one clade to other clades?\n all clades best fit lines')
plt.xlabel('average number of singletons per clade')
plt.ylabel('average number of SNPs transferred (normalized)')
plt.legend()
plt.show()

###
# Make a best fit line for the graph. For this fit through (0,0) origin
###
plot_bestfitline_origin(avg_singletons_basal_ab, avg_counts_basal_ab, 'r', 'basal AB clade')
plot_bestfitline_origin(avg_singletons_f, avg_counts_f, 'blue', 'F clade')
plot_bestfitline_origin(avg_singletons_topab, avg_counts_topab, 'g', 'top AB clade')
plot_bestfitline_origin(avg_singletons_mixing, avg_counts_mixing, 'm', 'mixing layer')

plt.title('What is the rate of transfer of fixed SNPs from one clade to other clades?\n all clades best fit lines through origin')
plt.xlabel('average number of singletons per clade')
plt.ylabel('average number of SNPs transferred (normalized)')
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
plot_bestfitline_exponential(avg_singletons_basal_ab, avg_counts_basal_ab, 'r', 'basal AB clade')
plot_bestfitline_exponential(avg_singletons_f, avg_counts_f, 'blue', 'F clade')
plot_bestfitline_exponential(avg_singletons_topab, avg_counts_topab, 'g', 'top AB clade')
plot_bestfitline_exponential(avg_singletons_mixing, avg_counts_mixing, 'm', 'mixing layer')

plt.title('What is the rate of transfer of fixed SNPs from one clade to other clades?\n all clades best fit lines weighted exponential')
plt.xlabel('average number of singletons per clade')
plt.ylabel('average number of SNPs transferred (normalized)')
plt.legend()
plt.show()







