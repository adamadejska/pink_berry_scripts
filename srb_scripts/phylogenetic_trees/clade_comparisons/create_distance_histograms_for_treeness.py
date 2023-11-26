####################################################################################
# A script that checks if we have improved the 'treeness' of the SRB tree when we
# get rid of loci that are fixed in F clade, basal AB clade, top AB clade.
#
# We calculate the distances between leaves of the same clade (clade leaf vs clade leaf).
# We create a histogram based on those distances. Then we calculate the distances between
# clade leafs and leafs outside the clade (clade leaf vs leaf outside). We create another
# histogram that shows those distances and compare those distributions. 
# 
# We expect them to be far apart if the clades are well defined. 
#
####################################################################################

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

def create_distance_matrix(current_strains, pos_index):
    # Create vector for each strain that contain positions only present in a particular cluster
    # Create strain_haplotypes for the positions in the cluster
    strain_data = {}
    for strain in current_strains:
        strain_positions = []

        # Change Nans to 0.5 to treat them probabilistically in the distance calculations
        for i in pos_index:
            snp = data.loc[strain].values[i]
            value = 0.5 if np.isnan(snp) else snp
            strain_positions.append(value)
            
        strain_data[strain] = strain_positions

    # Create a difference matrix where nans are considered 0.5 (treat it probabilistically)
    index = list(strain_data.keys())
    index.sort()
    distances = []
    for s1 in range(0, len(index)):
        strain1 = index[s1]
        a = np.array(strain_data[strain1])
        for s2 in range(s1+1, len(index)):
            strain2 = index[s2]
            b = np.array(strain_data[strain2])

            difference = sum(np.square(a-b))
            difference = difference / float(len(a))
            distances.append(difference)
    
    return distances


def calculate_distances_full_tree(clade_leaves, all_leaves, pos_index):
   
    # Calculate distances between a clade leaves and all other leaves in the clade
    strain_data = {}
    for strain in all_leaves:
        strain_positions = []

        # Change Nans to 0.5 to treat them probabilistically in the distance calculations
        for i in pos_index:
            snp = data.loc[strain].values[i]
            value = 0.5 if np.isnan(snp) else snp
            strain_positions.append(value)
            
        strain_data[strain] = strain_positions

    # Create a difference matrix where nans are considered 0.5 (treat it probabilistically)
    index = list(strain_data.keys())
    index.sort()
    distances_out = []
    for s1 in range(0, len(clade_leaves)):
        strain1 = clade_leaves[s1]
        a = np.array(strain_data[strain1])
        for s2 in range(0, len(index)):
            strain2 = index[s2]
            if  strain2 not in clade_leaves:
                b = np.array(strain_data[strain2])

                difference = sum(np.square(a-b))
                difference = difference / float(len(a))
                distances_out.append(difference)
                #print(strain1 + ' ' + strain2 + ' ' + str(difference))

    return distances_out

# Read in the SNP data
data = pd.read_csv('/home/ada/Desktop/Shraiman_lab/srb_data/srb_snp_data_2021.csv', index_col=0)

# Delete PB93 because it's a hyper mutator
data = data.drop('PB93')

f_strains = ['CAGAGAGG-AAGGAGTA','CAGAGAGG-GTAAGGAG','CTCTCTAC-TATCCTCT','AAGAGGCA-CTCTCTAT','AAGAGGCA-TATCCTCT','AAGAGGCA-CTAAGCCT','GCTACGCT-ACTGCATA','PB39','PB32','PB40','PB88','PB16','PB79','PB87','AAGAGGCA-TAGATCGC','AAGAGGCA-GTAAGGAG','PB64']
topAB_strains = ['TAAGGCGA-CTCTCTAT','TAAGGCGA-ACTGCATA','PB83','AGGCAGAA-CTCTCTAT','AGGCAGAA-TAGATCGC','PB53','GTAGAGGA-CTCTCTAT','TCCTGAGC-ACTGCATA','PB58','TCCTGAGC-AGAGTAGA','PB33','PB_4','PB84']
basalAB_strains = ['CGTACTAG-TATCCTCT','CGTACTAG-GTAAGGAG','PB50','AGGCAGAA-TATCCTCT','GTAGAGGA-TAGATCGC','TAAGGCGA-CTAAGCCT','PB82','PB90','PB27','AGGCAGAA-GTAAGGAG','CGTACTAG-CTAAGCCT','TAAGGCGA-TATCCTCT','PB18','CGTACTAG-AAGGAGTA','GTAGAGGA-AGAGTAGA','PB10','PB34','PB76']

clade1 = ['GGACTCCT-TAGATCGC', 'GGACTCCT-ACTGCATA', 'PB11', 'CTCTCTAC-TAGATCGC', 'PB77', 'PB35', 'PB29', 'GGACTCCT-CTAAGCCT', 'CAGAGAGG-CTCTCTAT', 'CTCTCTAC-CTCTCTAT', 'TCCTGAGC-CTAAGCCT', 'PB85', 'CAGAGAGG-CTAAGCCT']
clade2 = ['PB40', 'PB88', 'PB16', 'PB64', 'AAGAGGCA-TAGATCGC', 'AAGAGGCA-GTAAGGAG', 'PB87']
clade3 = ['CAGAGAGG-AAGGAGTA', 'CAGAGAGG-GTAAGGAG', 'CTCTCTAC-TATCCTCT', 'AAGAGGCA-TATCCTCT', 'AAGAGGCA-CTCTCTAT', 'AAGAGGCA-CTAAGCCT', 'GCTACGCT-ACTGCATA', 'PB39', 'CAGAGAGG-TAGATCGC', 'CAGAGAGG-AGAGTAGA', 'CAGAGAGG-ACTGCATA']
clade4 = ['CTCTCTAC-CTAAGCCT', 'CTCTCTAC-ACTGCATA', 'PB37', 'PB_5', 'CAGAGAGG-TATCCTCT', 'PB45', 'CTCTCTAC-AGAGTAGA']
clade5 = ['AGGCAGAA-AGAGTAGA', 'PB66', 'PB43', 'CGTACTAG-ACTGCATA', 'TCCTGAGC-TATCCTCT', 'CGTACTAG-TAGATCGC', 'TAAGGCGA-AAGGAGTA', 'CGTACTAG-AGAGTAGA', 'CGTACTAG-CTCTCTAT', 'PB_3', 'PB44', 'PB41', 'TCCTGAGC-GTAAGGAG', 'PB25']
clade6 = ['TAAGGCGA-ACTGCATA', 'TAAGGCGA-CTCTCTAT', 'AGGCAGAA-CTCTCTAT', 'AGGCAGAA-TAGATCGC', 'PB83', 'PB33', 'GTAGAGGA-CTCTCTAT', 'PB53']

whole_tree_no_hgt = ['GGACTCCT-TAGATCGC', 'GGACTCCT-ACTGCATA', 'PB11', 'CTCTCTAC-TAGATCGC', 'PB77', 'PB35', 'PB29', 'GGACTCCT-CTAAGCCT', 'CAGAGAGG-CTCTCTAT', 'CTCTCTAC-CTCTCTAT', 'TCCTGAGC-CTAAGCCT', 'PB85', 'CAGAGAGG-CTAAGCCT', 'PB40', 'PB88', 'PB16', 'PB64', 'AAGAGGCA-TAGATCGC', 'AAGAGGCA-GTAAGGAG', 'PB87', 'CAGAGAGG-AAGGAGTA', 'CAGAGAGG-GTAAGGAG', 'CTCTCTAC-TATCCTCT', 'AAGAGGCA-TATCCTCT', 'AAGAGGCA-CTCTCTAT', 'AAGAGGCA-CTAAGCCT', 'GCTACGCT-ACTGCATA', 'PB39', 'CAGAGAGG-TAGATCGC', 'CAGAGAGG-AGAGTAGA', 'CAGAGAGG-ACTGCATA','CTCTCTAC-CTAAGCCT', 'CTCTCTAC-ACTGCATA', 'PB37', 'PB_5', 'CAGAGAGG-TATCCTCT', 'PB45', 'CTCTCTAC-AGAGTAGA', 'PB57', 'PB21','AGGCAGAA-AGAGTAGA', 'PB66', 'PB43', 'CGTACTAG-ACTGCATA', 'TCCTGAGC-TATCCTCT', 'CGTACTAG-TAGATCGC', 'TAAGGCGA-AAGGAGTA', 'CGTACTAG-AGAGTAGA', 'CGTACTAG-CTCTCTAT', 'PB_3', 'PB44', 'PB41', 'TCCTGAGC-GTAAGGAG', 'PB25', 'PB_4', 'TAAGGCGA-ACTGCATA', 'TAAGGCGA-CTCTCTAT', 'AGGCAGAA-CTCTCTAT', 'AGGCAGAA-TAGATCGC', 'PB83', 'PB33', 'GTAGAGGA-CTCTCTAT', 'PB53','TCCTGAGC-AGAGTAGA', 'PB84', 'PB58', 'TCCTGAGC-ACTGCATA', 'PB42', 'CGTACTAG-CTAAGCCT']

data_index = data.index.values.tolist()  # names of bacterial samples
data_positions = data.columns.tolist()


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

# Create lists of distances for each clade
no_fixed_distances_1 = create_distance_matrix(clade1, pos_index_no_fixed)
no_fixed_distances_2 = create_distance_matrix(clade2, pos_index_no_fixed)
no_fixed_distances_3 = create_distance_matrix(clade3, pos_index_no_fixed)
no_fixed_distances_4 = create_distance_matrix(clade4, pos_index_no_fixed)
no_fixed_distances_5 = create_distance_matrix(clade5, pos_index_no_fixed)
no_fixed_distances_6 = create_distance_matrix(clade6, pos_index_no_fixed)

no_fixed_outside_distances_1 = calculate_distances_full_tree(clade1, whole_tree_no_hgt, pos_index_no_fixed)
no_fixed_outside_distances_2 = calculate_distances_full_tree(clade2, whole_tree_no_hgt, pos_index_no_fixed)
no_fixed_outside_distances_3 = calculate_distances_full_tree(clade2, whole_tree_no_hgt, pos_index_no_fixed)
no_fixed_outside_distances_4 = calculate_distances_full_tree(clade2, whole_tree_no_hgt, pos_index_no_fixed)
no_fixed_outside_distances_5 = calculate_distances_full_tree(clade2, whole_tree_no_hgt, pos_index_no_fixed)
no_fixed_outside_distances_6 = calculate_distances_full_tree(clade2, whole_tree_no_hgt, pos_index_no_fixed)

clade_distances = [no_fixed_distances_1, no_fixed_distances_2, no_fixed_distances_3, no_fixed_distances_4, no_fixed_distances_5, no_fixed_distances_6]
outside_distances = [no_fixed_outside_distances_1, no_fixed_outside_distances_2, no_fixed_outside_distances_3, no_fixed_outside_distances_4, no_fixed_outside_distances_5, no_fixed_outside_distances_6]
title = ['clade 1', 'clade 2', 'clade 3', 'clade 4', 'clade 5', 'clade 6']

for i in range(0, len(clade_distances)):
    bins = np.array(range(0, 36))/100
    plt.hist(outside_distances[i], alpha=0.5, bins=bins, label='outside leaves distances')
    plt.hist(clade_distances[i], alpha=0.5, bins=bins, label='within clade leaf distances')
    plt.title('Pairwise distances (clade leaf vs clade leaf) and (clade leaf vs outside leaf)\n' + title[i])
    plt.xlabel('pairwise distance')
    plt.ylabel('number of pairwise distances')
    plt.legend()
    plt.show()


