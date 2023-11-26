##############################################################
# Make a tree matrix based on fixed SRB F clade loci (frequency 20+ / 24)
# Rows: 1st-> chr number, next -> strains in the tree order, columns: fixed loci (sortedby contig position)
##############################################################


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def make_matrix(positions, strain_order, data):
# Make a dataframe that has data on the particular loci and strains specified.

	og_strain_order =  data.index.values.tolist()  # names of bacterial samples
	index = ['chr'] + strain_order

	# Create a dataframe to store info about snps in cluster strains
	df = pd.DataFrame(index=index)
	for i in positions:
		col = []
		col.append(data.loc['chr',i])
		for j in strain_order:
			vector = data[str(i)]
			vector[vector == 0] = -1
			vector = np.nan_to_num(vector)
			val = vector[og_strain_order.index(j)]
			col.append(val)
		df[i] = col

	return(df)


# Read in the SNP data
data = pd.read_csv('/home/ada/Desktop/Shraiman_lab/srb_data/srb_snp_data_2021_chr.csv', index_col=0)

data_index = data.index.values.tolist()  # names of bacterial samples
data_positions = data.columns.tolist()
"""
# List of strains that are in the srb f clade
srb_f_clade = ['CAGAGAGG-AAGGAGTA','CAGAGAGG-GTAAGGAG','CTCTCTAC-TATCCTCT','AAGAGGCA-CTCTCTAT','AAGAGGCA-TATCCTCT','AAGAGGCA-CTAAGCCT','GCTACGCT-ACTGCATA','PB39','PB32','PB40','PB88','PB16','PB79','PB87','AAGAGGCA-TAGATCGC','AAGAGGCA-GTAAGGAG','PB64']

# Create a mask that will only show the defined F clade values
mask = [True if x in srb_f_clade else False for x in data_index]

fixed_loci = []

# Check which loci are fixed in the SRB F clade. (more than 50%)
for i in data_positions:
	p = np.array(data[i].tolist())[mask]
	if sum(p==1)/(len(srb_f_clade)) >= 0.5:
		fixed_loci.append(i)

print(len(fixed_loci))

tree_order = ['CAGAGAGG-AAGGAGTA','CAGAGAGG-GTAAGGAG','CTCTCTAC-TATCCTCT','AAGAGGCA-CTCTCTAT','AAGAGGCA-TATCCTCT','AAGAGGCA-CTAAGCCT','GCTACGCT-ACTGCATA','PB39','PB32','PB40','PB88','PB16','PB79','PB87','AAGAGGCA-TAGATCGC','AAGAGGCA-GTAAGGAG','PB64','CAGAGAGG-TAGATCGC','CAGAGAGG-AGAGTAGA','CAGAGAGG-ACTGCATA','CTCTCTAC-CTAAGCCT','CTCTCTAC-ACTGCATA','PB37','PB_5','CAGAGAGG-TATCCTCT','CTCTCTAC-AGAGTAGA','CAGAGAGG-CTCTCTAT','CTCTCTAC-CTCTCTAT','PB93','PB11','CTCTCTAC-TAGATCGC','GGACTCCT-TAGATCGC','GGACTCCT-ACTGCATA','PB77','PB29','PB35','PB45','TCCTGAGC-CTAAGCCT','PB85','CAGAGAGG-CTAAGCCT','GGACTCCT-CTAAGCCT','PB57','PB21','TAAGGCGA-CTCTCTAT','TAAGGCGA-ACTGCATA','PB83','AGGCAGAA-CTCTCTAT','AGGCAGAA-TAGATCGC','PB53','GTAGAGGA-CTCTCTAT','TCCTGAGC-ACTGCATA','PB58','TCCTGAGC-AGAGTAGA','PB33','PB_4','PB84','TCCTGAGC-GTAAGGAG','PB41','PB25','GTAGAGGA-TATCCTCT','CGTACTAG-TATCCTCT','CGTACTAG-GTAAGGAG','PB50','AGGCAGAA-TATCCTCT','GTAGAGGA-TAGATCGC','TAAGGCGA-CTAAGCCT','PB82','PB90','PB27','AGGCAGAA-GTAAGGAG','CGTACTAG-CTAAGCCT','TAAGGCGA-TATCCTCT','PB18','CGTACTAG-AAGGAGTA','GTAGAGGA-AGAGTAGA','PB10','PB34','PB76','CGTACTAG-ACTGCATA','TCCTGAGC-TATCCTCT','CGTACTAG-AGAGTAGA','CGTACTAG-TAGATCGC','AGGCAGAA-AGAGTAGA','PB66','TAAGGCGA-AAGGAGTA','PB43','CGTACTAG-CTCTCTAT','PB_3','PB44','AGGCAGAA-AAGGAGTA','PB42']

# Get rid of strains that have more than 30% missing data in the designated positions. 
pos_mask = [True if x in fixed_loci else False for x in data_positions]



snp_df = make_matrix(fixed_loci, tree_order, data)

file = '/home/ada/Desktop/Shraiman_lab/srb_data/snp_matrix_fixed_f_clade_threshold50_NaNs30wholegenome_fixed.csv'
snp_df.to_csv(file)
"""
##############################################################
# Make the same SNP matrix for fixed loci in the AB clade.

srb_basal_ab_clade = ['CGTACTAG-TATCCTCT','CGTACTAG-GTAAGGAG','PB50','AGGCAGAA-TATCCTCT','GTAGAGGA-TAGATCGC','TAAGGCGA-CTAAGCCT','PB82','PB90','PB27','AGGCAGAA-GTAAGGAG','CGTACTAG-CTAAGCCT','TAAGGCGA-TATCCTCT','PB18','CGTACTAG-AAGGAGTA','GTAGAGGA-AGAGTAGA','PB10','PB34','PB76']
srb_top_ab_clade = ['TAAGGCGA-ACTGCATA','PB83','AGGCAGAA-CTCTCTAT','AGGCAGAA-TAGATCGC','PB53','GTAGAGGA-CTCTCTAT','TCCTGAGC-ACTGCATA','PB58','TCCTGAGC-AGAGTAGA','PB33','PB_4','PB84','TCCTGAGC-GTAAGGAG','PB41','PB25','GTAGAGGA-TATCCTCT']

# Create a mask that will only show the defined F clade values
mask = [True if x in srb_top_ab_clade else False for x in data_index]

fixed_loci = []

# Check which loci are fixed in the SRB basal/top AB clade. (more than 50%)
for i in data_positions:
	p = np.array(data[i].tolist())[mask]
	if sum(p==1) >= len(srb_top_ab_clade)/2.0:
		fixed_loci.append(i)

print(len(fixed_loci))

tree_order = ['CAGAGAGG-AAGGAGTA','CAGAGAGG-GTAAGGAG','CTCTCTAC-TATCCTCT','AAGAGGCA-CTCTCTAT','AAGAGGCA-TATCCTCT','AAGAGGCA-CTAAGCCT','GCTACGCT-ACTGCATA','PB39','PB32','PB40','PB88','PB16','PB79','PB87','AAGAGGCA-TAGATCGC','AAGAGGCA-GTAAGGAG','PB64','CAGAGAGG-TAGATCGC','CAGAGAGG-AGAGTAGA','CAGAGAGG-ACTGCATA','CTCTCTAC-CTAAGCCT','CTCTCTAC-ACTGCATA','PB37','PB_5','CAGAGAGG-TATCCTCT','CTCTCTAC-AGAGTAGA','CAGAGAGG-CTCTCTAT','CTCTCTAC-CTCTCTAT','PB93','PB11','CTCTCTAC-TAGATCGC','GGACTCCT-TAGATCGC','GGACTCCT-ACTGCATA','PB77','PB29','PB35','PB45','TCCTGAGC-CTAAGCCT','PB85','CAGAGAGG-CTAAGCCT','GGACTCCT-CTAAGCCT','PB57','PB21','TAAGGCGA-CTCTCTAT','TAAGGCGA-ACTGCATA','PB83','AGGCAGAA-CTCTCTAT','AGGCAGAA-TAGATCGC','PB53','GTAGAGGA-CTCTCTAT','TCCTGAGC-ACTGCATA','PB58','TCCTGAGC-AGAGTAGA','PB33','PB_4','PB84','TCCTGAGC-GTAAGGAG','PB41','PB25','GTAGAGGA-TATCCTCT','CGTACTAG-TATCCTCT','CGTACTAG-GTAAGGAG','PB50','AGGCAGAA-TATCCTCT','GTAGAGGA-TAGATCGC','TAAGGCGA-CTAAGCCT','PB82','PB90','PB27','AGGCAGAA-GTAAGGAG','CGTACTAG-CTAAGCCT','TAAGGCGA-TATCCTCT','PB18','CGTACTAG-AAGGAGTA','GTAGAGGA-AGAGTAGA','PB10','PB34','PB76','CGTACTAG-ACTGCATA','TCCTGAGC-TATCCTCT','CGTACTAG-AGAGTAGA','CGTACTAG-TAGATCGC','AGGCAGAA-AGAGTAGA','PB66','TAAGGCGA-AAGGAGTA','PB43','CGTACTAG-CTCTCTAT','PB_3','PB44','AGGCAGAA-AAGGAGTA','PB42']

# Get rid of strains that have more than 30% missing data in the designated positions. 
pos_mask = [True if x in fixed_loci else False for x in data_positions]

new_tree_order = []
for s in tree_order:
	vector = np.array(data.loc[s])
	vector[vector == 0] = -1
	vector = np.nan_to_num(vector)
	if sum(vector==0)/float(len(vector)) <= 0.3:
		new_tree_order.append(s)


snp_df = make_matrix(fixed_loci, tree_order, data)

file = '/home/ada/Desktop/Shraiman_lab/srb_data/snp_matrix_fixed_TOP_AB_clade_threshold50_NaNs30wholegenome.csv'
snp_df.to_csv(file)



