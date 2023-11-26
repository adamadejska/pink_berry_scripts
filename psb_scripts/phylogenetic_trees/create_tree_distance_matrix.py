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
# For the PSB data, there is only one possible versions that can be made:
# no 8 clade loci.
####################################################################################

import numpy as np
import pandas as pd


# Read in the SNP data
data = pd.read_csv('/home/ada/Desktop/Shraiman_lab/data/snp_data_2020.csv', index_col=0)
clade_8_strains = ['PB87','PB40','TCCTGAGC-TAGATCGC','TCCTGAGC-CTCTCTAT','PB73','AAGAGGCA-TATCCTCT','GGACTCCT-AGAGTAGA','GGACTCCT-CTAAGCCT']

data_index = data.index.values.tolist()  # names of bacterial samples
data_positions = data.columns.tolist()

# Only consider positions that are not singletons and that are not fixed in our current test clade.
positions = []
counter_used = 0
for i in data_positions:
    
    # Find positions that are singletons -> we don't want those.
    n = sum(np.array(data[i])==1)
    if n > 1:
        positions.append(i)
        counter_used += 1

print('used: ' + str(counter_used))
print(len(data_positions))

# Find the indicies for the positions in the dataframe
pos_index = [data_positions.index(str(i)) for i in positions]

strains_kept = []
# Check how many strains have very low number of SNPs
for strain in data_index:
    strain_positions = []
    ns_num = 0.0

    # Change Nans to 0.5 to treat them probabilistically in the distance calculations
    for i in pos_index:
        snp = data.loc[strain].values[i]
        value = 0.5 if np.isnan(snp) else snp
        if value == 0.5:
            ns_num += 1

        strain_positions.append(value)
    
    if float(sum(np.array(strain_positions)==1))/len(strain_positions) >= 0.01:
        strains_kept.append(strain)
    
    
    #print(strain + ' ' + str(sum(np.array(strain_positions)==1)) + 
    #      ' / ' + str(len(strain_positions)) + ' ' + 
    #      str(float(sum(np.array(strain_positions)==1))/len(strain_positions)))

# Create vector for each strain that contain positions only present in a particular cluster
# Create strain_haplotypes for the positions in the cluster
#strains_kept = data_index
strain_data = {}
for strain in strains_kept:
    strain_positions = []
    ns_num = 0.0

    # Change Nans to 0.5 to treat them probabilistically in the distance calculations
    for i in pos_index:
        snp = data.loc[strain].values[i]
        value = 0.5 if np.isnan(snp) else snp
        if value == 0.5:
            ns_num += 1

        strain_positions.append(value)
        
    # Add number of missing data to the name
    #name = strain + '_' + str(ns_num)
    strain_data[strain] = strain_positions
    
    if ns_num/len(strain_positions) >= 0.2:
        print(strain + ' ' + str(ns_num/len(strain_positions)))
    """
    if ns_num/len(strain_positions) <= 0.3:   # Get rid of strains with a lot of missing data
        #strain_data[name] = strain_positions
        strain_data[strain] = strain_positions
        #print(strain)
    """
    

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
out1 = open('/home/ada/Desktop/Shraiman_lab/data/tree_data/PSB_names_no1s.csv','w')
out2 = open('/home/ada/Desktop/Shraiman_lab/data/tree_data/PSB_values_no1s.csv', 'w')
out1.write(','.join(index))
out2.write(','.join(matlab_data))
out1.close()
out2.close()

