###########################################################################
# Estimate the recombination rate. 
###########################################################################

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Info on haplotypes and dataset
snp_matrix_path = '/home/ada/Desktop/Shraiman_lab/srb_data/snp_matrix_fixed_F_basalAB_clade_threshold50_NaNs30wholegenome_fixed.csv'
basal_AB_strains = ['CGTACTAG-TATCCTCT','CGTACTAG-GTAAGGAG','PB50','AGGCAGAA-TATCCTCT','GTAGAGGA-TAGATCGC','TAAGGCGA-CTAAGCCT','PB82','PB90','PB27','AGGCAGAA-GTAAGGAG','CGTACTAG-CTAAGCCT','TAAGGCGA-TATCCTCT','PB18','CGTACTAG-AAGGAGTA','GTAGAGGA-AGAGTAGA','PB10','PB34','PB76']
f_clade_strains = ['chr', 'CAGAGAGG-AAGGAGTA','CAGAGAGG-GTAAGGAG','CTCTCTAC-TATCCTCT','AAGAGGCA-CTCTCTAT','AAGAGGCA-TATCCTCT','AAGAGGCA-CTAAGCCT','GCTACGCT-ACTGCATA','PB39','PB32','PB40','PB88','PB16','PB79','PB87','AAGAGGCA-TAGATCGC','AAGAGGCA-GTAAGGAG','PB64']

# Read in the data matrix
snp_matrix = pd.read_csv(snp_matrix_path, index_col=0)
columns = snp_matrix.columns.tolist()
strains = snp_matrix.index.tolist()[1:]

# Make a matrix containing only F clade loci. We'll use it to check which columns belong to the F matrix (aka which columns have
# SNPs that are majority in the F clade.)
f_matrix = snp_matrix.loc[f_clade_strains, :]

# We will store distances between loci based on their colors (either the change of color or all distances in the matrix)
# (conditional and unconditional arrays)
transition_distances = []
all_distances = []
current_chr = snp_matrix.iloc[0, 0]   # [row_position, column_position]
f_n = len(f_clade_strains)   # number of F clade strains.

big_dist = []
# Iterrate through the matrix, column by column.
for i in range(1, len(columns)-1):

    if current_chr == snp_matrix.iloc[0, i+1]:   # Check if the next column belong to the same chr.
        col1 = np.array(f_matrix.iloc[:, i])
        col2 = np.array(f_matrix.iloc[:, i+1])

        if sum(col1==1) >= f_n/2 and sum(col2==1) >= f_n/2:    # Check if the current columns belong to the F clade fixed matrix.

            for s in range(1, len(strains)):                # Go through all strains (all rows), not just the F clade strains
                current_snp = snp_matrix.iloc[s, i]
                next_snp = snp_matrix.iloc[s, i+1]

                # Check if we transition from one color to another.
                if (current_snp == -1 and next_snp != -1) or (current_snp != -1 and next_snp == -1):  
                    dist = int(float(columns[i+1])) - int(float(columns[i]))
                    transition_distances.append(dist)
                
                
                # Add all distances to the unconditional array.
                dist = int(float(columns[i+1])) - int(float(columns[i]))
                all_distances.append(dist)

                if dist > 10000:
                    big_dist.append(columns[i])


    current_chr = snp_matrix.iloc[0, i+1]

#print(list(set(big_dist)))
# Bin the distances into bins of size 100 bp.
# Calculate the cumulant for the green -> white and white -> green distances of F clade
t_bin_height, bins, _ = plt.hist(transition_distances, bins=range(0,21000, 100), edgecolor='black')
all_bin_height, bins, _ = plt.hist(all_distances, bins=range(0,21000, 100), edgecolor='black')
t_bin_height = np.array(t_bin_height)
all_bin_height = np.array(all_bin_height)

print(list(t_bin_height))
print(list(all_bin_height))
print(list(bins))

plt.clf()
plt.hist(all_distances, bins=range(0,21000, 100), edgecolor='black', label='unconditional')
plt.hist(transition_distances, bins=range(0,21000, 100), edgecolor='black', label='transitions')
plt.legend()
plt.xlabel('distance (bp) (bin_size=100bp)')
plt.ylabel('counts of loci pairs (log)')
plt.yscale('log')
plt.show()
all_bin_height = all_bin_height + 1   # add pseudocounts to avoid division by zero

rtl = -1 * np.log(1 - (t_bin_height / all_bin_height))

plt.plot(bins[:-1], rtl)
plt.ylabel('rtl = -log(1 - (n/N))')
plt.xlabel('distance (bp) (bin_size=100bp)')
plt.title('Estimation of the recombination rate for the SRB SNP matrix')
plt.show()