###############################################################
# Ada Madejska, 2023
# Create a base tree for PSB. Given our analysis of SNP frequency 
# across the genome, create a tree based on only the loci that
# were in low frequency regions in that previous analysis.
# We want to establish some kind of base using well conserved regions.  
###############################################################

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Read in the PSB snp data
data = pd.read_csv('/home/ada/Desktop/Shraiman_lab/srb_data/srb_snp_data_2021.csv', index_col=0)

# Delete PB93 because it's a hyper mutator
data = data.drop('PB93')

# Index = names of bacterial samples
data_index = data.index.values.tolist()
data_positions = data.columns.tolist()

# Count the number of SNPs across all strains in each 1 kb window.
window = 1000
densities = {}
windows = []
for p in data_positions:
    p_int = p.split('.')[0]
    int_div = int(p_int) // window
    if int_div not in densities.keys():
        densities[int_div] = 1
        windows.append(int_div)
    else:
        densities[int_div] += 1

# Plot the densities across the genome
windows = np.array(windows) / window

# Filter the windows based on the number of mutations they contain
# We want only the windows that have less than 10 SNPs per 1 kb window (for PSB)

low_density_regions = []
for window, density in densities.items():
    if density <= 20:
        low_density_regions.append(window)

print('finding low density loci')
# Find the loci that are in the low density regions. We will use them for distance matrix construction.
positions = []
for p in data_positions:
    p_int = p.split('.')[0]
    int_div = int(p_int) // window
    if int_div in low_density_regions:
        positions.append(p)

print('creating matrix')
####
# Create a distance matrix to be used in MATlab to create a tree
####

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
    strain_data[strain] = strain_positions
    
    if ns_num/len(strain_positions) <= 0.3:   # Get rid of strains with a lot of missing data
        strain_data[strain] = strain_positions


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
out1 = open('/home/ada/Desktop/PinkBerry_scripts_paper/data/srb/SRB_names_baseline.csv','w')
out2 = open('/home/ada/Desktop/PinkBerry_scripts_paper/data/srb/SRB_values_baseline.csv', 'w')
out1.write(','.join(index))
out2.write(','.join(matlab_data))
out1.close()
out2.close()