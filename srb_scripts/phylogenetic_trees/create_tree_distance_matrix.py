####################################################################################
# Create files that will be used in the MATLAB script for trees creation from clusters
# Take big clusters and based exclusively on the positions present in each cluster
# create a tree of all samples
# Exclude strains with no mutations at all positions in the cluster.
# Produces two files: strain names and similarity matrix.
####################################################################################

import numpy as np
import pandas as pd
import re


# Read in the SNP data
data = pd.read_csv('/home/ada/Desktop/Shraiman_lab/srb_data/srb_snp_data_2021.csv', index_col=0)

data_index = data.index.values.tolist()  # names of bacterial samples
data_positions = data.columns.tolist()

# Only consider positions that are not singletons or doubletons.
positions = []
for i in data_positions:
	n = sum(np.array(data[i])==1)
	if n > 2:
		positions.append(i)

print(len(positions))

# Find the indicies for the positions in the dataframe
pos_index = [data_positions.index(str(i)) for i in positions]

# Create vector for each strain that contain positions only present in a particular cluster
# Create strain_haplotypes for the positions in the cluster
strain_data = {}
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
	
	# Add number of missing data to the name
	#name = strain + '_' + str(ns_num)

	#print(strain + ' ' + str(ns_num/len(strain_positions)))
	if ns_num/len(strain_positions) <= 0.3:   # Get rid of strains with a lot of missing data
		#strain_data[name] = strain_positions
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
out1 = open('/home/ada/Desktop/PinkBerry_scripts_paper/data/srb/SRB_names_no1s_2023.csv','w')
out2 = open('/home/ada/Desktop/PinkBerry_scripts_paper/data/srb/SRB_values_no1s_2023.csv', 'w')
out1.write(','.join(index))
out2.write(','.join(matlab_data))
out1.close()
out2.close()

