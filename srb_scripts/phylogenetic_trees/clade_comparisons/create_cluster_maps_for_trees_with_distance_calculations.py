####################################################################################
# Create files that will be used in the MATLAB script for trees creation from clusters
# Take big clusters and based exclusively on the positions present in each cluster
# create a tree of all samples
# Exclude strains with no mutations at all positions in the cluster.
# Produces a distance matrix.
#
# For the SRB data we need to calculate those distances for no fixed F clade, no fixed basal
# AB clade, and no fixed top AB clade. Each tree also has different clades based on 
# the tree. 
#
# This script focuses on tree without fixed all clades SNPs. 
####################################################################################

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.spatial.distance import pdist, squareform


def create_distance_matrix(current_strains, pos_index):
    # Create vector for each strain that contain positions only present in a particular cluster
    # Create strain_haplotypes for the positions in the cluster
    strain_data = {}
    for strain in current_strains:
        strain_positions = []
        ns_num = 0.0

        # Change Nans to 0.5 to treat them probabilistically in the distance calculations
        for i in pos_index:
            snp = data.loc[strain].values[i]
            value = 0.5 if np.isnan(snp) else snp
            if value == 0.5:
                ns_num += 1

            strain_positions.append(value)
            
        if ns_num/len(strain_positions) <= 0.3:   # Get rid of strains with a lot of missing data
            strain_data[strain] = strain_positions

    # Create a difference matrix where nans are considered 0.5 (treat it probabilistically)
    index = list(strain_data.keys())
    index.sort()

    x = []
    for strain in index:
        x.append(strain_data[strain])

    # Calculate distances
    distvec = pdist(x, 'correlation')
    m = squareform(distvec)
    
    return m



# Read in the SNP data
data = pd.read_csv('/home/ada/Desktop/Shraiman_lab/srb_data/srb_snp_data_2021.csv', index_col=0)

# Delete PB93 because it's a hyper mutator
data = data.drop('PB93')

f_strains = ['CAGAGAGG-AAGGAGTA','CAGAGAGG-GTAAGGAG','CTCTCTAC-TATCCTCT','AAGAGGCA-CTCTCTAT','AAGAGGCA-TATCCTCT','AAGAGGCA-CTAAGCCT','GCTACGCT-ACTGCATA','PB39','PB32','PB40','PB88','PB16','PB79','PB87','AAGAGGCA-TAGATCGC','AAGAGGCA-GTAAGGAG','PB64']
topAB_strains = ['TAAGGCGA-CTCTCTAT','TAAGGCGA-ACTGCATA','PB83','AGGCAGAA-CTCTCTAT','AGGCAGAA-TAGATCGC','PB53','GTAGAGGA-CTCTCTAT','TCCTGAGC-ACTGCATA','PB58','TCCTGAGC-AGAGTAGA','PB33','PB_4','PB84']
basalAB_strains = ['CGTACTAG-TATCCTCT','CGTACTAG-GTAAGGAG','PB50','AGGCAGAA-TATCCTCT','GTAGAGGA-TAGATCGC','TAAGGCGA-CTAAGCCT','PB82','PB90','PB27','AGGCAGAA-GTAAGGAG','CGTACTAG-CTAAGCCT','TAAGGCGA-TATCCTCT','PB18','CGTACTAG-AAGGAGTA','GTAGAGGA-AGAGTAGA','PB10','PB34','PB76']

all_strains_full_tree = ['AAGAGGCA-CTAAGCCT','AAGAGGCA-CTCTCTAT','AAGAGGCA-GTAAGGAG','AAGAGGCA-TAGATCGC','AAGAGGCA-TATCCTCT','AGGCAGAA-AAGGAGTA','AGGCAGAA-AGAGTAGA','AGGCAGAA-CTCTCTAT','AGGCAGAA-GTAAGGAG','AGGCAGAA-TAGATCGC','AGGCAGAA-TATCCTCT','CAGAGAGG-AAGGAGTA','CAGAGAGG-ACTGCATA','CAGAGAGG-AGAGTAGA','CAGAGAGG-CTAAGCCT','CAGAGAGG-CTCTCTAT','CAGAGAGG-GTAAGGAG','CAGAGAGG-TAGATCGC','CAGAGAGG-TATCCTCT','CGTACTAG-AAGGAGTA','CGTACTAG-ACTGCATA','CGTACTAG-AGAGTAGA','CGTACTAG-CTAAGCCT','CGTACTAG-CTCTCTAT','CGTACTAG-GTAAGGAG','CGTACTAG-TAGATCGC','CGTACTAG-TATCCTCT','CTCTCTAC-ACTGCATA','CTCTCTAC-AGAGTAGA','CTCTCTAC-CTAAGCCT','CTCTCTAC-CTCTCTAT','CTCTCTAC-TAGATCGC','CTCTCTAC-TATCCTCT','GCTACGCT-ACTGCATA','GGACTCCT-ACTGCATA','GGACTCCT-CTAAGCCT','GGACTCCT-TAGATCGC','GTAGAGGA-AGAGTAGA','GTAGAGGA-CTCTCTAT','GTAGAGGA-TAGATCGC','GTAGAGGA-TATCCTCT','PB10','PB11','PB16','PB18','PB21','PB25','PB27','PB29','PB32','PB33','PB34','PB35','PB37','PB39','PB40','PB41','PB42','PB43','PB44','PB45','PB50','PB53','PB57','PB58','PB64','PB66','PB76','PB77','PB79','PB82','PB83','PB84','PB85','PB87','PB88','PB90','PB_3','PB_4','PB_5','TAAGGCGA-AAGGAGTA','TAAGGCGA-ACTGCATA','TAAGGCGA-CTAAGCCT','TAAGGCGA-CTCTCTAT','TAAGGCGA-TATCCTCT','TCCTGAGC-ACTGCATA','TCCTGAGC-AGAGTAGA','TCCTGAGC-CTAAGCCT','TCCTGAGC-GTAAGGAG','TCCTGAGC-TATCCTCT']
all_strains_no_all_clades_tree = ['AAGAGGCA-CTAAGCCT','AAGAGGCA-CTCTCTAT','AAGAGGCA-GTAAGGAG','AAGAGGCA-TAGATCGC','AAGAGGCA-TATCCTCT','AGGCAGAA-AGAGTAGA','AGGCAGAA-CTCTCTAT','AGGCAGAA-TAGATCGC','CAGAGAGG-AAGGAGTA','CAGAGAGG-ACTGCATA','CAGAGAGG-AGAGTAGA','CAGAGAGG-CTAAGCCT','CAGAGAGG-CTCTCTAT','CAGAGAGG-GTAAGGAG','CAGAGAGG-TAGATCGC','CAGAGAGG-TATCCTCT','CGTACTAG-ACTGCATA','CGTACTAG-AGAGTAGA','CGTACTAG-CTAAGCCT','CGTACTAG-CTCTCTAT','CGTACTAG-TAGATCGC','CTCTCTAC-ACTGCATA','CTCTCTAC-AGAGTAGA','CTCTCTAC-CTAAGCCT','CTCTCTAC-CTCTCTAT','CTCTCTAC-TAGATCGC','CTCTCTAC-TATCCTCT','GCTACGCT-ACTGCATA','GGACTCCT-ACTGCATA','GGACTCCT-CTAAGCCT','GGACTCCT-TAGATCGC','GTAGAGGA-CTCTCTAT','PB11','PB16','PB21','PB25','PB29','PB33','PB35','PB37','PB39','PB40','PB41','PB42','PB43','PB44','PB45','PB53','PB57','PB58','PB64','PB66','PB77','PB83','PB84','PB85','PB87','PB88','PB_3','PB_4','PB_5','TAAGGCGA-AAGGAGTA','TAAGGCGA-ACTGCATA','TAAGGCGA-CTCTCTAT','TCCTGAGC-ACTGCATA','TCCTGAGC-AGAGTAGA','TCCTGAGC-CTAAGCCT','TCCTGAGC-GTAAGGAG','TCCTGAGC-TATCCTCT']

data_index = data.index.values.tolist()  # names of bacterial samples
data_positions = data.columns.tolist()

# Only consider positions that are not singletons. Include all other positions
positions = []
counter_used = 0
for i in data_positions:
    
    # Find positions that are singletons -> we don't want those.
    n = sum(np.array(data[i])==1)
    if n > 2:
        positions.append(i)
        counter_used += 1

# Find the indicies for the positions in the dataframe
pos_index = [data_positions.index(str(i)) for i in positions]


full_distances = create_distance_matrix(all_strains_full_tree, pos_index)

g = sns.clustermap(full_distances, yticklabels=True, xticklabels=True, metric='correlation', method='single')
g.fig.suptitle('SRB full SNP matrix similarity cluster map')
plt.show()

# Only consider positions that are not singletons and that are not fixed in our current test clade.
positions = []
counter_used = 0
# Create a mask that will only show the defined test clade values
mask_f = [True if x in f_strains else False for x in data_index]
mask_bAB = [True if x in basalAB_strains else False for x in data_index]
mask_tAB = [True if x in topAB_strains else False for x in data_index]
for i in data_positions:
    
    # Find positions that are singletons -> we don't want those.
    n = sum(np.array(data[i])==1)
    if n > 2:
        # Calculate the frequency and find positions that are fixed in the current clade -> we don't want those.
        frequency_f = sum(np.array(data[i].tolist())[mask_f]==1) / float(len(f_strains))
        frequency_bAB = sum(np.array(data[i].tolist())[mask_bAB]==1) / float(len(basalAB_strains))
        frequency_tAB = sum(np.array(data[i].tolist())[mask_tAB]==1) / float(len(topAB_strains))

        if frequency_f < 0.5 and frequency_bAB < 0.5 and frequency_tAB < 0.5:
            positions.append(i)
            counter_used += 1


# Find the indicies for the positions in the dataframe
pos_index_no_fixed = [data_positions.index(str(i)) for i in positions]


no_fixed_distances = create_distance_matrix(all_strains_no_all_clades_tree, pos_index_no_fixed)

g = sns.clustermap(no_fixed_distances, yticklabels=True, xticklabels=True, metric='correlation', method='single')
g.fig.suptitle('SRB no fixed F / topAB / basalAB clades SNPs matrix similarity cluster map')
plt.show()
