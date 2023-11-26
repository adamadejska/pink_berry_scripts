####################################################################################
# Create files that will be used in the MATLAB script for trees creation from clusters
# Take big clusters and based exclusively on the positions present in each cluster
# create a tree of all samples
# Exclude strains with no mutations at all positions in the cluster.
# Produces a distance matrix.
#
# For the PSB data, there is only one possible versions that can be made:
# no 8 clade loci.
####################################################################################

import matplotlib.pyplot as plt
import math
import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform
import seaborn as sns

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

    # Instead of calculating the distances manually, we will use pdist and its different methods.
    # Each strain's SNPs will be considered separate dimentions
    # Create a list of lists x where each list inside is the information about separate strain. 
    index = list(strain_data.keys())
    index.sort()

    x = []
    for strain in index:
        x.append(strain_data[strain])

    # Calculate distances
    distvec = pdist(x, 'cosine')
    m = squareform(distvec)
    print(m)
    
    return m



# Read in the SNP data
data = pd.read_csv('/home/ada/Desktop/Shraiman_lab/data/snp_data_2020.csv', index_col=0)
all_strains_full_tree = ['AAGAGGCA-AAGGAGTA','AAGAGGCA-CTAAGCCT','AAGAGGCA-CTCTCTAT','AAGAGGCA-GTAAGGAG','AAGAGGCA-TAGATCGC','AGGCAGAA-AGAGTAGA','CAGAGAGG-AAGGAGTA','CAGAGAGG-ACTGCATA','CAGAGAGG-AGAGTAGA','CAGAGAGG-CTAAGCCT','CAGAGAGG-CTCTCTAT','CAGAGAGG-GTAAGGAG','CAGAGAGG-TAGATCGC','CAGAGAGG-TATCCTCT','CTCTCTAC-AAGGAGTA','CTCTCTAC-ACTGCATA','CTCTCTAC-AGAGTAGA','CTCTCTAC-CTAAGCCT','CTCTCTAC-CTCTCTAT','CTCTCTAC-GTAAGGAG','CTCTCTAC-TAGATCGC','CTCTCTAC-TATCCTCT','GCTACGCT-ACTGCATA','GGACTCCT-ACTGCATA','GGACTCCT-AGAGTAGA','GGACTCCT-CTAAGCCT','GGACTCCT-GTAAGGAG','GGACTCCT-TAGATCGC','GTAGAGGA-CTAAGCCT','PB11','PB13','PB16','PB24','PB26','PB32','PB34','PB37','PB39','PB40','PB44','PB45','PB47','PB48','PB52','PB55','PB59','PB61','PB64','PB67','PB69','PB73','PB77','PB78','PB80','PB85','PB87','PB_5','PB_8','TAAGGCGA-TATCCTCT','TCCTGAGC-CTCTCTAT','TCCTGAGC-TAGATCGC']
all_strains_no_8_clade_tree = ['AAGAGGCA-AAGGAGTA','AAGAGGCA-CTAAGCCT','AAGAGGCA-CTCTCTAT','AAGAGGCA-GTAAGGAG','AAGAGGCA-TAGATCGC','AGGCAGAA-AGAGTAGA','CAGAGAGG-AAGGAGTA','CAGAGAGG-ACTGCATA','CAGAGAGG-AGAGTAGA','CAGAGAGG-CTAAGCCT','CAGAGAGG-CTCTCTAT','CAGAGAGG-GTAAGGAG','CAGAGAGG-TAGATCGC','CAGAGAGG-TATCCTCT','CTCTCTAC-AAGGAGTA','CTCTCTAC-ACTGCATA','CTCTCTAC-AGAGTAGA','CTCTCTAC-CTAAGCCT','CTCTCTAC-CTCTCTAT','CTCTCTAC-GTAAGGAG','CTCTCTAC-TAGATCGC','CTCTCTAC-TATCCTCT','GCTACGCT-ACTGCATA','GGACTCCT-ACTGCATA','GGACTCCT-AGAGTAGA','GGACTCCT-CTAAGCCT','GGACTCCT-GTAAGGAG','GGACTCCT-TAGATCGC','GTAGAGGA-CTAAGCCT','PB11','PB13','PB16','PB24','PB26','PB32','PB34','PB37','PB39','PB40','PB44','PB45','PB47','PB48','PB52','PB59','PB61','PB64','PB67','PB69','PB73','PB77','PB78','PB80','PB85','PB87','PB_5','PB_8','TAAGGCGA-TATCCTCT','TCCTGAGC-CTCTCTAT','TCCTGAGC-TAGATCGC']

clade_8_strains = ['PB87','PB40','TCCTGAGC-TAGATCGC','TCCTGAGC-CTCTCTAT','PB73','AAGAGGCA-TATCCTCT','GGACTCCT-AGAGTAGA','GGACTCCT-CTAAGCCT']

data_index = data.index.values.tolist()  # names of bacterial samples
data_positions = data.columns.tolist()

# Only consider positions that are not singletons. Include all other positions
positions = []
counter_used = 0
for i in data_positions:
    
    # Find positions that are singletons -> we don't want those.
    n = sum(np.array(data[i])==1)
    if n > 1:
        positions.append(i)
        counter_used += 1

# Find the indicies for the positions in the dataframe
pos_index = [data_positions.index(str(i)) for i in positions]


full_distances = create_distance_matrix(all_strains_full_tree, pos_index)

g = sns.clustermap(full_distances, yticklabels=True, xticklabels=True, metric='cosine', method='median')
g.fig.suptitle('PSB full SNP matrix similarity cluster map')
plt.show()

# Only consider positions that are not singletons and that are not fixed in our current test clade.
positions = []
counter_used = 0
for i in data_positions:
    
    # Find positions that are singletons -> we don't want those.
    n = sum(np.array(data[i])==1)
    if n > 1:
        total = np.nansum(np.array(data.loc[clade_8_strains,i]))
        if total < len(clade_8_strains)/2:
            positions.append(i)
            counter_used += 1


# Find the indicies for the positions in the dataframe
pos_index_no_fixed = [data_positions.index(str(i)) for i in positions]

no_fixed_distances = create_distance_matrix(all_strains_no_8_clade_tree, pos_index)

g = sns.clustermap(no_fixed_distances, yticklabels=True, xticklabels=True, metric='cosine', method='median')
g.fig.suptitle('PSB no fixed 8 clade SNPs matrix similarity cluster map')
plt.show()