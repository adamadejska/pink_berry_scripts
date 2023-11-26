####################################################################################
# Create files that will be used in the MATLAB script for trees creation from clusters
# Take big clusters and based exclusively on the positions present in each cluster
# create a tree of all samples
# Exclude strains with no mutations at all positions in the cluster.
# Produces two files: strain names and similarity matrix.
#
# This script makes trees based on data matrix without the fixed loci (the ones that
# we consider to be involved in horizontal transfer). 
#
# For the SRB data, the are multiple possible versions that can be made:
# no F clade loci, no topAB clade loci, and no basalAB clade loci.
####################################################################################

import numpy as np
import pandas as pd


# Read in the SNP data
data = pd.read_csv('/home/ada/Desktop/Shraiman_lab/srb_data/srb_snp_data_2021.csv', index_col=0)

# Delete PB93 because it's a hyper mutator
data = data.drop('PB93')

f_strains = ['CAGAGAGG-AAGGAGTA','CAGAGAGG-GTAAGGAG','CTCTCTAC-TATCCTCT','AAGAGGCA-CTCTCTAT','AAGAGGCA-TATCCTCT','AAGAGGCA-CTAAGCCT','GCTACGCT-ACTGCATA','PB39','PB32','PB40','PB88','PB16','PB79','PB87','AAGAGGCA-TAGATCGC','AAGAGGCA-GTAAGGAG','PB64']
topAB_strains = ['TAAGGCGA-CTCTCTAT','TAAGGCGA-ACTGCATA','PB83','AGGCAGAA-CTCTCTAT','AGGCAGAA-TAGATCGC','PB53','GTAGAGGA-CTCTCTAT','TCCTGAGC-ACTGCATA','PB58','TCCTGAGC-AGAGTAGA','PB33','PB_4','PB84']
basalAB_strains = ['CGTACTAG-TATCCTCT','CGTACTAG-GTAAGGAG','PB50','AGGCAGAA-TATCCTCT','GTAGAGGA-TAGATCGC','TAAGGCGA-CTAAGCCT','PB82','PB90','PB27','AGGCAGAA-GTAAGGAG','CGTACTAG-CTAAGCCT','TAAGGCGA-TATCCTCT','PB18','CGTACTAG-AAGGAGTA','GTAGAGGA-AGAGTAGA','PB10','PB34','PB76']

data_index = data.index.values.tolist()  # names of bacterial samples
data_positions = data.columns.tolist()

# Create a mask that will only show the defined test clade values
mask_f = [True if x in f_strains else False for x in data_index]
mask_bAB = [True if x in basalAB_strains else False for x in data_index]
mask_tAB = [True if x in topAB_strains else False for x in data_index]

# Only consider positions that are not singletons and that are not fixed in our current test clade.
positions = []
counter_fixed = 0
counter_used = 0
for i in data_positions:
    
    # Calculate the frequency and find positions that are fixed in the current clade -> we don't want those.
    frequency_f = sum(np.array(data[i].tolist())[mask_f]==1) / float(len(f_strains))
    frequency_bAB = sum(np.array(data[i].tolist())[mask_bAB]==1) / float(len(basalAB_strains))
    frequency_tAB = sum(np.array(data[i].tolist())[mask_tAB]==1) / float(len(topAB_strains))

    n = sum(np.array(data[i])==1)
    if n > 2:
        if frequency_f < 0.5 and frequency_bAB < 0.5 and frequency_tAB < 0.5:
            positions.append(i)
            counter_used += 1
        else:
            counter_fixed += 1

print('fixed: ' + str(counter_fixed))
print('used: ' + str(counter_used))
print(len(data_positions))

# Find the indicies for the positions in the dataframe
pos_index = [data_positions.index(str(i)) for i in positions]

# Create vector for each strain that contain positions only present in a particular cluster
# Create strain_haplotypes for the positions in the cluster
strain_data = {}
for strain in data_index:
    strain_positions = []
    ns_num = 0.0
    snp_num = 1.0

    # Change Nans to 0.5 to treat them probabilistically in the distance calculations
    for i in pos_index:
        snp = data.loc[strain].values[i]
        value = 0.5 if np.isnan(snp) else snp
        if value == 0.5:
            ns_num += 1
        elif value == 1: # Check the number of actual SNPs in each strain
            snp_num += 1

        strain_positions.append(value)

    

        
    # Add number of missing data to the name
    #name = strain + '_' + str(ns_num)

    #print(strain + ' ' + str(ns_num/len(strain_positions)))
    if ns_num/len(strain_positions) <= 0.3 and snp_num/len(strain_positions) >= 0.005:   
        # Get rid of strains with a lot of missing data and not a lot of SNPs
        strain_data[strain] = strain_positions
        print(strain + ' ' + str(snp_num))


# Add an out group to the strain data
strain_data['outgroup'] = [0]*len(positions)

# Create a difference matrix where nans are considered 0.5 (treat it probabilistically)
index = list(strain_data.keys())
index.sort()
sim_df = pd.DataFrame(index=index)
for s1 in index:
    a = np.array(strain_data[s1])
    col_a = []
    for s2 in index:
        b = np.array(strain_data[s2])

        difference = sum(np.square(a-b))
        difference = difference / float(len(a))
        col_a.append(difference)
    sim_df[s1] = col_a

matlab_data = []
for i in range(0, len(index)):
    for j in range(i+1, len(index)):
        temp = (sim_df.loc[index[j],index[i]])
        matlab_data.append(temp)

# Save the names and values to two csv files that will be inported to matlab
# to create trees.

matlab_data = [str(i) for i in matlab_data]
out1 = open('/home/ada/Desktop/PinkBerry_scripts_paper/data/srb/tree_matrixes/SRB_names_no1s_no_all_clades.csv','w')
out2 = open('/home/ada/Desktop/PinkBerry_scripts_paper/data/srb/tree_matrixes/SRB_values_no1s_no_all_clades.csv', 'w')
out1.write(','.join(index))
out2.write(','.join(matlab_data))
out1.close()
out2.close()
