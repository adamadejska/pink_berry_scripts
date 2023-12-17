###########################################################################
# Using the information from the haplotyping algorithm, we can figure
# out which loci are wild type. In this script we will figure out
# which loci are wild type in all strains. We will use this information
# for building trees
###########################################################################

import ast
import numpy as np
import pandas as pd

# Read in the wt haplotype file
wt_hap_file = '/home/ada/Desktop/PinkBerry_scripts_paper/psb_scripts/haplotyping/haplotypes_full_data.txt'

# Create a numpy matrix to store the data divided by the strain
m = []
with open(wt_hap_file, 'r') as f:
    for line in f:
        if not line.startswith('#'):    # we found a line with wt regions listed
            #line = line[1:-1]
            list_of_regions = list(ast.literal_eval(line))
            wt_positions = np.array([])
            for i in list_of_regions:
                wt_positions = np.append(wt_positions, range(int(i[0]), int(i[1])))   # Make a list of all possible loci within the WT regions
            
            m.append(wt_positions)

print('figuring out all possible common loci')
# Find common loci between all strains
common = m[0]
for i in m:
    common = np.intersect1d(common, i)

print('found all common possible loci')
# Read in the SNP data
data = pd.read_csv('/home/ada/Desktop/Shraiman_lab/data/snp_data_2020.csv', index_col=0)
clade_8_strains = ['PB87','PB40','TCCTGAGC-TAGATCGC','TCCTGAGC-CTCTCTAT','PB73','AAGAGGCA-TATCCTCT','GGACTCCT-AGAGTAGA','GGACTCCT-CTAAGCCT']

data_index = data.index.values.tolist()  # names of bacterial samples
data_positions = np.array([int(i) for i in data.columns.tolist()])

common_actual_loci = np.intersect1d(data_positions, common)

print(len(common_actual_loci))

# Only consider positions that are not singletons and that are not fixed in our current test clade.
positions = [str(i).split('.')[0] for i in common_actual_loci.tolist()]
data_positions = [str(i).split('.')[0] for i in data_positions.tolist()]

# Find the indicies for the positions in the dataframe
pos_index = [data_positions.index(i) for i in positions]


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
    
# Create vector for each strain that contain positions only present in a particular cluster
# Create strain_haplotypes for the positions in the cluster
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
out1 = open('/home/ada/Desktop/PinkBerry_scripts_paper/data/psb/tree_matrices/PSB_names_wt_regions.csv','w')
out2 = open('/home/ada/Desktop/PinkBerry_scripts_paper/data/psb/tree_matrices/PSB_values_wt_region.csv', 'w')
out1.write(','.join(index))
out2.write(','.join(matlab_data))
out1.close()
out2.close()
