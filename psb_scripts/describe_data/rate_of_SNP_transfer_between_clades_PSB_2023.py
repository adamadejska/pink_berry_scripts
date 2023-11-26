###############################################################
# What is the frequency of fixed clade loci in different clades?
# Figure out the fixed loci in each clade
# and calculate frequencies of mutations at those positions
# for all other strains.
# PSB ANALYSIS
###############################################################

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def find_fixed_loci(clade, data_index, data_positions):

	# Create a mask that will only show the defined F clade values
	mask = [True if x in clade else False for x in data_index]
	fixed_loci = []

	# Check SNP frequency in the current clade.
	for i in data_positions:
		p = np.array(data[i].tolist())[mask]
		if sum(p==1) >= len(clade)/2:  				
			fixed_loci.append(i)

	return(fixed_loci)


# Read in the SNP data
# Read in the PSB snp data
data = pd.read_csv('/home/ada/Desktop/Shraiman_lab/data/snp_data_2020.csv', index_col=0)

# Index = names of bacterial samples
data_index = data.index.values.tolist()
data_positions = data.columns.tolist()

eight_clade = []
f_clade = []
e_clade = []
bc_clade = []
full_tree = []

working_clade = eight_clade

psb_singleton_counts = {'AAGAGGCA-AAGGAGTA': 483, 'PB36': 481, 'PB80': 646, 'PB23': 607, 'PB78': 547, 
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

# Find all fixed topAB loci
fixed = find_fixed_loci(working_clade, data_index, data_positions)

# Calculate the average counts of the fixed topAB loci outside of the topAB clade
# AND plot it against the avg number of singletons for each clade.
# Whole tree
mask = [True if x in fixed else False for x in data_positions]
# Order: F, mixing layer, basalAB, basal
avg_counts, avg_singletons = [], []
clades = [eight_clade, f_clade, bc_clade, e_clade]

for i in clades:
    counts, singletons = [], []
    for s in i:
        vector = np.array(data.loc[s][1:][mask])
        counts.append(sum(vector==1))
        singletons.append(psb_singleton_counts[s])

    # Avg the counts and singletons

    counts, singletons = np.array(counts), np.array(singletons)
    avg_counts.append(np.mean(counts))
    avg_singletons.append(np.mean(singletons))


print(avg_counts)
print(avg_singletons)
plt.scatter(avg_singletons,avg_counts)
plt.title('What is the rate of transfer of fixed SNPs from one clade to other clades?\n 8 clade')
plt.xlabel('average number of singletons per clade')
plt.ylabel('average number of SNPs transferred form 8 clade')
plt.show()
















