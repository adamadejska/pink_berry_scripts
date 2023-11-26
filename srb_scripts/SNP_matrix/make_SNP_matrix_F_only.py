##############################################################
# Make a SNP matrix for loci that are fixed (50%+) in F clade and basal_AB clade.
# Rows: 1st-> chr number, next -> strains in the tree order, columns: fixed loci (sortedby contig position)
##############################################################


import numpy as np
import pandas as pd


# Read in the SNP data
data = pd.read_csv('/home/ada/Desktop/Shraiman_lab/srb_data/srb_snp_data_2021_chr.csv', index_col=0)

data_index = data.index.values.tolist()  # names of bacterial samples
data_positions = data.columns.tolist()

# List of strains that are in the srb f clade
f_clade_strains = ['CAGAGAGG-AAGGAGTA','CAGAGAGG-GTAAGGAG','CTCTCTAC-TATCCTCT','AAGAGGCA-CTCTCTAT','AAGAGGCA-TATCCTCT','AAGAGGCA-CTAAGCCT','GCTACGCT-ACTGCATA','PB39','PB32','PB40','PB88','PB16','PB79','PB87','AAGAGGCA-TAGATCGC','AAGAGGCA-GTAAGGAG','PB64']

# Create a mask that will only show the defined F clade values
mask_F = [True if x in f_clade_strains else False for x in data_index]

fixed_loci = []
clade = []
chr = []

# Check which loci are fixed in the SRB F clade. (more than 50%)
for i in range(0,len(data_positions)):

    p = np.array(data[data_positions[i]].tolist())[mask_F]
    if sum(p==1)/(len(f_clade_strains)) >= 0.5:
        fixed_loci.append(data_positions[i])
        clade.append('F')
        chr.append(data.iloc[0,i])
		

tree_order = ['CAGAGAGG-AAGGAGTA','CAGAGAGG-GTAAGGAG','CTCTCTAC-TATCCTCT','AAGAGGCA-CTCTCTAT','AAGAGGCA-TATCCTCT','AAGAGGCA-CTAAGCCT','GCTACGCT-ACTGCATA','PB39','PB32','PB40','PB88','PB16','PB79','PB87','AAGAGGCA-TAGATCGC','AAGAGGCA-GTAAGGAG','PB64','CAGAGAGG-TAGATCGC','CAGAGAGG-AGAGTAGA','CAGAGAGG-ACTGCATA','CTCTCTAC-CTAAGCCT','CTCTCTAC-ACTGCATA','PB37','PB_5','CAGAGAGG-TATCCTCT','CTCTCTAC-AGAGTAGA','CAGAGAGG-CTCTCTAT','CTCTCTAC-CTCTCTAT','PB93','PB11','CTCTCTAC-TAGATCGC','GGACTCCT-TAGATCGC','GGACTCCT-ACTGCATA','PB77','PB29','PB35','PB45','TCCTGAGC-CTAAGCCT','PB85','CAGAGAGG-CTAAGCCT','GGACTCCT-CTAAGCCT','PB57','PB21','TAAGGCGA-CTCTCTAT','TAAGGCGA-ACTGCATA','PB83','AGGCAGAA-CTCTCTAT','AGGCAGAA-TAGATCGC','PB53','GTAGAGGA-CTCTCTAT','TCCTGAGC-ACTGCATA','PB58','TCCTGAGC-AGAGTAGA','PB33','PB_4','PB84','TCCTGAGC-GTAAGGAG','PB41','PB25','GTAGAGGA-TATCCTCT','CGTACTAG-TATCCTCT','CGTACTAG-GTAAGGAG','PB50','AGGCAGAA-TATCCTCT','GTAGAGGA-TAGATCGC','TAAGGCGA-CTAAGCCT','PB82','PB90','PB27','AGGCAGAA-GTAAGGAG','CGTACTAG-CTAAGCCT','TAAGGCGA-TATCCTCT','PB18','CGTACTAG-AAGGAGTA','GTAGAGGA-AGAGTAGA','PB10','PB34','PB76','CGTACTAG-ACTGCATA','TCCTGAGC-TATCCTCT','CGTACTAG-AGAGTAGA','CGTACTAG-TAGATCGC','AGGCAGAA-AGAGTAGA','PB66','TAAGGCGA-AAGGAGTA','PB43','CGTACTAG-CTCTCTAT','PB_3','PB44','AGGCAGAA-AAGGAGTA','PB42']
index = ['chr'] + tree_order
og_strain_order =  data.index.values.tolist()  # names of bacterial samples

# Create a dataframe to store info about snps in cluster strains
snp_df = pd.DataFrame(index=index)
for i in range(0,len(fixed_loci)):
	col = []
	col.append(chr[i])
	for j in tree_order:
		vector = data[str(fixed_loci[i])]
		vector[vector == 0] = -1
		vector = np.nan_to_num(vector)

		# Distinguish between F clade SNPs and basal_AB SNPs
		# F clade SNP = 1, basal_AB SNP = 2
		if clade[i] == 'F':	
			vector[vector == 1] = 1

		val = vector[og_strain_order.index(j)]
		col.append(val)
	snp_df[fixed_loci[i]] = col

file = '/home/ada/Desktop/Shraiman_lab/srb_data/snp_matrix_fixed_F_ONLY_threshold50_NaNs30wholegenome_fixed.csv'
snp_df.to_csv(file)