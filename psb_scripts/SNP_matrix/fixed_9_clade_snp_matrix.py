#####################################################################
# Given the SNP data for the PSB, make a SNP matrix based on the tree.
# Every column is a locus in the genome (sorted) and every row is a different
# bacteria (ordered by the tree order).
# Focus only on SNPs that are considered fixed (over 50% in a clade) in 
# the 9 clade. 
# Ignore strains with a lot of missing data (over 40%).  
#####################################################################


import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Open the original file with all snps
snp_matrix_path = '/home/ada/Desktop/Shraiman_lab/data/snp_data_2020.csv'
snp_matrix = pd.read_csv(snp_matrix_path, index_col=0)
columns = snp_matrix.columns.tolist()
strains = snp_matrix.index.tolist()

clade_9_strains = ['PB87','PB40','GGACTCCT-AGAGTAGA','GGACTCCT-GTAAGGAG','GGACTCCT-CTAAGCCT','TCCTGAGC-TAGATCGC','TCCTGAGC-CTCTCTAT','PB73','PB80']
strains_order = ['PB87','PB40','GGACTCCT-AGAGTAGA','GGACTCCT-GTAAGGAG','GGACTCCT-CTAAGCCT','TCCTGAGC-TAGATCGC','TCCTGAGC-CTCTCTAT','PB73','PB80','AAGAGGCA-TATCCTCT','PB63','AAGAGGCA-TAGATCGC','PB39','PB_8','GCTACGCT-GTAAGGAG','PB24','PB55','AAGAGGCA-ACTGCATA','AAGAGGCA-AGAGTAGA','AAGAGGCA-GTAAGGAG','AAGAGGCA-AAGGAGTA','GCTACGCT-ACTGCATA','PB31','AAGAGGCA-CTAAGCCT','AAGAGGCA-CTCTCTAT','PB32','PB16','PB47','PB48','PB64','AGGCAGAA-ACTGCATA','AGGCAGAA-CTAAGCCT','TAAGGCGA-TATCCTCT','PB78','PB61','GGACTCCT-ACTGCATA','GGACTCCT-TAGATCGC','PB45','PB77','PB69','PB_5','PB37','PB13','PB85','PB28','PB67','PB76','PB59','PB52','PB84','TAAGGCGA-AAGGAGTA','PB26','PB34','PB11','CTCTCTAC-AAGGAGTA','GTAGAGGA-CTAAGCCT','CTCTCTAC-CTAAGCCT','CTCTCTAC-CTCTCTAT','CTCTCTAC-TAGATCGC','CTCTCTAC-ACTGCATA','CAGAGAGG-TAGATCGC','CAGAGAGG-CTAAGCCT','CAGAGAGG-GTAAGGAG','CTCTCTAC-TATCCTCT','CAGAGAGG-ACTGCATA','CAGAGAGG-AGAGTAGA','CAGAGAGG-TATCCTCT','CAGAGAGG-AAGGAGTA','CTCTCTAC-AGAGTAGA','CTCTCTAC-GTAAGGAG','PB44','PB36','PB58','PB82','PB20','PB60','PB_3','PB51','PB_2','PB25','CGTACTAG-CTCTCTAT','PB41','CGTACTAG-TAGATCGC','PB57','PB_9','PB43','AGGCAGAA-CTCTCTAT','PB12','PB66','TAAGGCGA-AGAGTAGA','PB_1','PB92','PB90','PB21','PB18','PB35','TCCTGAGC-ACTGCATA','PB42','TCCTGAGC-CTAAGCCT','TCCTGAGC-AGAGTAGA','GTAGAGGA-CTCTCTAT','TCCTGAGC-TATCCTCT','CGTACTAG-ACTGCATA','CGTACTAG-GTAAGGAG','PB75','TCCTGAGC-AAGGAGTA','GTAGAGGA-GTAAGGAG']

fixed_snps = []
for i in columns:
    total = np.nansum(np.array(snp_matrix.loc[clade_9_strains,i]))
    if total > len(clade_9_strains)/2:
        fixed_snps.append(i)

fixed_matrix = snp_matrix.loc[strains_order, fixed_snps]
fixed_matrix = fixed_matrix.reindex(index = strains_order)
fixed_matrix = fixed_matrix.replace(0, -1)
fixed_matrix = fixed_matrix.fillna(0)

strains = fixed_matrix.index.tolist()

# Check the number of NaNs in each strain, get rid of any strains that have a lot of missing data: 40%+
nansum_total = []
for s in strains:
    nansum = float(np.sum(np.array(fixed_matrix.loc[s,:])==0))
    if nansum/len(np.array(fixed_matrix.loc[s,:])) > 0.4:
        fixed_matrix = fixed_matrix.drop(s)

strains = ['PB87','PB40','GGACTCCT-AGAGTAGA','GGACTCCT-GTAAGGAG','GGACTCCT-CTAAGCCT','TCCTGAGC-TAGATCGC','TCCTGAGC-CTCTCTAT','PB73','PB80','AAGAGGCA-TATCCTCT','PB63','AAGAGGCA-TAGATCGC','PB39','PB_8','GCTACGCT-GTAAGGAG','PB24','PB55','AAGAGGCA-GTAAGGAG','AAGAGGCA-AAGGAGTA','GCTACGCT-ACTGCATA','PB31','AAGAGGCA-CTAAGCCT','AAGAGGCA-CTCTCTAT','PB32','PB16','PB47','PB48','PB64','AGGCAGAA-ACTGCATA','AGGCAGAA-CTAAGCCT','TAAGGCGA-TATCCTCT','PB78','PB61','GGACTCCT-ACTGCATA','GGACTCCT-TAGATCGC','PB45','PB77','PB69','PB_5','PB37','PB13','PB85','PB28','PB67','PB76','PB59','PB52','PB84','TAAGGCGA-AAGGAGTA','PB26','PB34','PB11','CTCTCTAC-AAGGAGTA','GTAGAGGA-CTAAGCCT','CTCTCTAC-CTAAGCCT','CTCTCTAC-CTCTCTAT','CTCTCTAC-TAGATCGC','CTCTCTAC-ACTGCATA','CAGAGAGG-TAGATCGC','CAGAGAGG-CTAAGCCT','CAGAGAGG-GTAAGGAG','CTCTCTAC-TATCCTCT','CAGAGAGG-ACTGCATA','CAGAGAGG-AGAGTAGA','CAGAGAGG-TATCCTCT','CAGAGAGG-AAGGAGTA','CTCTCTAC-AGAGTAGA','CTCTCTAC-GTAAGGAG','PB44','PB_3','CGTACTAG-CTCTCTAT','PB_9','AGGCAGAA-CTCTCTAT','TAAGGCGA-AGAGTAGA','PB_1','PB92','PB18','PB35','TCCTGAGC-ACTGCATA','PB42','TCCTGAGC-CTAAGCCT','TCCTGAGC-AGAGTAGA','TCCTGAGC-TATCCTCT','CGTACTAG-ACTGCATA','CGTACTAG-GTAAGGAG','TCCTGAGC-AAGGAGTA','GTAGAGGA-GTAAGGAG']


# Check the number of fixed 9 clade SNPs, if the number is very low, we would like to ignore those strains (they are very close to the root)
wt_total = []
for s in strains:
    wtsum = float(np.sum(np.array(fixed_matrix.loc[s,:])==-1))
    nansum = float(np.sum(np.array(fixed_matrix.loc[s,:])==0))
    fraction = (wtsum+nansum)/len(np.array(fixed_matrix.loc[s,:]))
    if fraction >= 0.99:
        fixed_matrix = fixed_matrix.drop(s)

print(len(fixed_matrix.columns.tolist()))
fixed_matrix.to_csv('/home/ada/Desktop/PinkBerry_scripts_paper/data/psb/fixed_9_clade_snp_matrix.csv', index=True)
