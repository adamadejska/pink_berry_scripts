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
import numpy as np
import pandas as pd
from scipy.stats import ks_2samp

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

    # Add an out group to the strain data
    strain_data['outgroup'] = [0]*len(positions)

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


def plot_boxplots(df):
    vals, names, xs = [],[],[]
    for i, col in enumerate(df.columns):
        vals.append(df[col].values)
        names.append(col)
        xs.append(np.random.normal(i + 1, 0.04, df[col].values.shape[0])) 
        # adds jitter to the data points - can be adjusted

    plt.boxplot(vals, labels=names)
    palette = ['r', 'g']

    for x, val, c in zip(xs, vals, palette):
        plt.scatter(x, val, alpha=0.4, color=c)


# Read in the SNP data
data = pd.read_csv('/home/ada/Desktop/Shraiman_lab/data/snp_data_2020.csv', index_col=0)
clade_8_strains = ['PB87','PB40','TCCTGAGC-TAGATCGC','TCCTGAGC-CTCTCTAT','PB73','AAGAGGCA-TATCCTCT','GGACTCCT-AGAGTAGA','GGACTCCT-CTAAGCCT']
clade_1 = ['GGACTCCT-AGAGTAGA', 'GGACTCCT-GTAAGGAG', 'GGACTCCT-CTAAGCCT', 'TCCTGAGC-TAGATCGC', 'TCCTGAGC-CTCTCTAT', 'PB73', 'PB67', 'PB59', 'PB52', 'PB44', 'AGGCAGAA-AGAGTAGA']
clade_2 = ['PB78', 'PB61', 'GGACTCCT-TAGATCGC', 'GGACTCCT-ACTGCATA', 'PB45', 'PB77', 'PB69', 'PB37', 'PB_5', 'PB13', 'PB85', 'TAAGGCGA-TATCCTCT']
clade_3 = ['AAGAGGCA-CTAAGCCT', 'AAGAGGCA-CTCTCTAT', 'PB32', 'PB16', 'PB47', 'AAGAGGCA-GTAAGGAG', 'GCTACGCT-ACTGCATA', 'PB87', 'PB40', 'PB80', 'AAGAGGCA-TAGATCGC', 'PB39', 'PB_8', 'PB48', 'PB64']
clade_4 = ['CTCTCTAC-AAGGAGTA','GTAGAGGA-CTAAGCCT','CTCTCTAC-CTAAGCCT','CTCTCTAC-CTCTCTAT','CTCTCTAC-TAGATCGC','CTCTCTAC-ACTGCATA','CAGAGAGG-TAGATCGC','CAGAGAGG-CTAAGCCT','CAGAGAGG-GTAAGGAG','CTCTCTAC-TATCCTCT','CAGAGAGG-ACTGCATA','CAGAGAGG-AGAGTAGA','CAGAGAGG-TATCCTCT','CAGAGAGG-AAGGAGTA','CTCTCTAC-AGAGTAGA','CTCTCTAC-GTAAGGAG']


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
print(len(pos_index))

full_distances_1 = create_distance_matrix(clade_1, pos_index)
full_distances_2 = create_distance_matrix(clade_2, pos_index)
full_distances_3 = create_distance_matrix(clade_3, pos_index)
full_distances_4 = create_distance_matrix(clade_4, pos_index)

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
print(len(pos_index_no_fixed))

# Create lists of distances for each clade
no_fixed_distances_1 = create_distance_matrix(clade_1, pos_index_no_fixed)
no_fixed_distances_2 = create_distance_matrix(clade_2, pos_index_no_fixed)
no_fixed_distances_3 = create_distance_matrix(clade_3, pos_index_no_fixed)
no_fixed_distances_4 = create_distance_matrix(clade_4, pos_index_no_fixed)

# Plot
df = pd.DataFrame({'Full SNP data clade 1': full_distances_1, 'No fixed \n8 clade SNPs clade 1': no_fixed_distances_1})
plot_boxplots(df)
plt.title('Pairwise distances used in tree building \n clade 1')
#plt.show()

df = pd.DataFrame({'Full SNP data clade 2': full_distances_2, 'No fixed \n8 clade SNPs clade 2': no_fixed_distances_2})
plot_boxplots(df)
plt.title('Pairwise distances used in tree building \n clade 2')
#plt.show()

df = pd.DataFrame({'Full SNP data clade 3': full_distances_3, 'No fixed \n8 clade SNPs clade 3': no_fixed_distances_3})
plot_boxplots(df)
plt.title('Pairwise distances used in tree building \n clade 3')
#plt.show()

df = pd.DataFrame({'Full SNP data clade 4': full_distances_4, 'No fixed \n8 clade SNPs clade 4': no_fixed_distances_4})
plot_boxplots(df)
plt.title('Pairwise distances used in tree building \n clade 4')
#plt.show()


# Statistical test
# Under the null hypothesis the two distributions are identical.
print(ks_2samp(full_distances_1, no_fixed_distances_1, alternative='less'))
print(ks_2samp(full_distances_2, no_fixed_distances_2, alternative='less'))
print(ks_2samp(full_distances_3, no_fixed_distances_3, alternative='less'))
print(ks_2samp(full_distances_4, no_fixed_distances_4, alternative='less'))
