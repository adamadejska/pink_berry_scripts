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
            strain_positions.append(value)
            
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


data_index = data.index.values.tolist()  # names of bacterial samples
data_positions = data.columns.tolist()

# Only consider positions that are not singletons. Include all other positions
positions = []
counter_used = 0
for i in data_positions:
    
    # Find positions that are singletons ot doubletons -> we don't want those.
    n = sum(np.array(data[i])==1)
    if n > 2:
        positions.append(i)
        counter_used += 1

# Find the indicies for the positions in the dataframe
pos_index = [data_positions.index(str(i)) for i in positions]
print(len(pos_index))
full_distances_1 = create_distance_matrix(clade1, pos_index)
full_distances_2 = create_distance_matrix(clade2, pos_index)
full_distances_3 = create_distance_matrix(clade3, pos_index)
full_distances_4 = create_distance_matrix(clade4, pos_index)
full_distances_5 = create_distance_matrix(clade5, pos_index)
full_distances_6 = create_distance_matrix(clade6, pos_index)

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
print(len(pos_index_no_fixed))

# Plot
df = pd.DataFrame({'Full SNP data clade 1': full_distances_1, 'No fixed \nF / top AB / basal AB clade SNPs\n clade 1': no_fixed_distances_1})
plot_boxplots(df)
plt.title('Pairwise distances used in tree building \n clade 1')
plt.show()

df = pd.DataFrame({'Full SNP data clade 2': full_distances_2, 'No fixed \nF / top AB / basal AB clade SNPs\n clade 2': no_fixed_distances_2})
plot_boxplots(df)
plt.title('Pairwise distances used in tree building \n clade 2')
plt.show()

df = pd.DataFrame({'Full SNP data clade 3': full_distances_3, 'No fixed \nF / top AB / basal AB clade SNPs\n clade 3': no_fixed_distances_3})
plot_boxplots(df)
plt.title('Pairwise distances used in tree building \n clade 3')
plt.show()

df = pd.DataFrame({'Full SNP data clade 4': full_distances_4, 'No fixed \nF / top AB / basal AB clade SNPs\n clade 4': no_fixed_distances_4})
plot_boxplots(df)
plt.title('Pairwise distances used in tree building \n clade 4')
plt.show()

df = pd.DataFrame({'Full SNP data clade 5': full_distances_5, 'No fixed \nF / top AB / basal AB clade SNPs\n clade 5': no_fixed_distances_5})
plot_boxplots(df)
plt.title('Pairwise distances used in tree building \n clade 5')
plt.show()

df = pd.DataFrame({'Full SNP data clade 6': full_distances_6, 'No fixed \nF / top AB / basal AB clade SNPs\n clade 6': no_fixed_distances_6})
plot_boxplots(df)
plt.title('Pairwise distances used in tree building \n clade 6')
plt.show()

# Statistical test
# Under the null hypothesis the two distributions are identical.
print(ks_2samp(full_distances_1, no_fixed_distances_1, alternative='less'))
print(ks_2samp(full_distances_2, no_fixed_distances_2, alternative='less'))
print(ks_2samp(full_distances_3, no_fixed_distances_3, alternative='less'))
print(ks_2samp(full_distances_4, no_fixed_distances_4, alternative='less'))
print(ks_2samp(full_distances_5, no_fixed_distances_5, alternative='less'))
print(ks_2samp(full_distances_6, no_fixed_distances_6, alternative='less'))

