##############################################
# Create a bar graph that shows how many of the loci that are fixed in the 
# 8 clade are also present in other strains across the tree.
##############################################



import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def find_fixed_loci(clade, data_index, data_positions):

	# Create a mask that will only show the defined working clade values
	mask = [True if x in clade else False for x in data_index]
	fixed_loci = []

	# Check SNP frequency in the current clade.
	for i in data_positions:
		p = np.array(data[i].tolist())[mask]
		if sum(p==1) >= len(clade)/2:  				
			fixed_loci.append(i)

	return(fixed_loci)


# Read in the SNP data
data = pd.read_csv('/home/ada/Desktop/Shraiman_lab/data/snp_data_2020.csv', index_col=0)

data_index = data.index.values.tolist()  # names of bacterial samples
data_positions = data.columns.tolist()

clade_9_strains = ['PB87','PB40','TCCTGAGC-TAGATCGC','TCCTGAGC-CTCTCTAT','PB73','AAGAGGCA-TATCCTCT','GGACTCCT-AGAGTAGA','GGACTCCT-CTAAGCCT']
ladder = ['GGACTCCT-CTAAGCCT','PB80','PB63','AAGAGGCA-TAGATCGC','PB39','PB_8','PB24','PB64','PB55','AAGAGGCA-AGAGTAGA','GCTACGCT-GTAAGGAG','AAGAGGCA-ACTGCATA','AAGAGGCA-AAGGAGTA','GCTACGCT-ACTGCATA','PB31','AAGAGGCA-GTAAGGAG','AAGAGGCA-CTAAGCCT','AAGAGGCA-CTCTCTAT','PB32','PB16','PB47','PB48']
bc_strains = ['TAAGGCGA-TATCCTCT','PB59','PB45','PB77','GGACTCCT-TAGATCGC','GGACTCCT-ACTGCATA','PB78','PB61','PB13','PB85','PB37','PB_5','AGGCAGAA-ACTGCATA','AGGCAGAA-CTAAGCCT','PB25','PB67','PB76','PB59','PB84','PB52','PB44']
e_clade = ['PB11','CTCTCTAC-AAGGAGTA','GTAGAGGA-CTAAGCCT','CTCTCTAC-CTAAGCCT','CAGAGAGG-CTCTCTAT','CTCTCTAC-CTCTCTAT','CTCTCTAC-TAGATCGC','CTCTCTAC-ACTGCATA','CAGAGAGG-TAGATCGC','CAGAGAGG-CTAAGCCT','CAGAGAGG-GTAAGGAG','CTCTCTAC-TATCCTCT','CAGAGAGG-ACTGCATA','CAGAGAGG-AGAGTAGA','CAGAGAGG-TATCCTCT','CAGAGAGG-AAGGAGTA','CTCTCTAC-AGAGTAGA','CTCTCTAC-GTAAGGAG','PB34','PB26']
basal_clade = ['PB29','CGTACTAG-TATCCTCT','PB57','PB51','AGGCAGAA-GTAAGGAG','PB90','PB36','PB82','PB20','PB65','PB60','PB21','PB58','PB43','PB66','GTAGAGGA-AAGGAGTA','PB12','AGGCAGAA-AAGGAGTA','GTAGAGGA-CTCTCTAT','PB75','AGGCAGAA-AGAGTAGA','PB41','CGTACTAG-TAGATCGC','PB25','PB_2']
strains = ['PB87','PB40','TCCTGAGC-TAGATCGC','TCCTGAGC-CTCTCTAT','PB73','AAGAGGCA-TATCCTCT','GGACTCCT-AGAGTAGA','GGACTCCT-CTAAGCCT','GGACTCCT-CTAAGCCT','PB80','PB63','AAGAGGCA-TAGATCGC','PB39','PB_8','PB24','PB64','PB55','AAGAGGCA-AGAGTAGA','GCTACGCT-GTAAGGAG','AAGAGGCA-ACTGCATA','AAGAGGCA-AAGGAGTA','GCTACGCT-ACTGCATA','PB31','AAGAGGCA-GTAAGGAG','AAGAGGCA-CTAAGCCT','AAGAGGCA-CTCTCTAT','PB32','PB16','PB47','PB48','TAAGGCGA-TATCCTCT','PB59','PB45','PB77','GGACTCCT-TAGATCGC','GGACTCCT-ACTGCATA','PB78','PB61','PB13','PB85','PB37','PB_5','AGGCAGAA-ACTGCATA','AGGCAGAA-CTAAGCCT','PB25','PB67','PB76','PB59','PB84','PB52','PB44','PB11','CTCTCTAC-AAGGAGTA','GTAGAGGA-CTAAGCCT','CTCTCTAC-CTAAGCCT','CAGAGAGG-CTCTCTAT','CTCTCTAC-CTCTCTAT','CTCTCTAC-TAGATCGC','CTCTCTAC-ACTGCATA','CAGAGAGG-TAGATCGC','CAGAGAGG-CTAAGCCT','CAGAGAGG-GTAAGGAG','CTCTCTAC-TATCCTCT','CAGAGAGG-ACTGCATA','CAGAGAGG-AGAGTAGA','CAGAGAGG-TATCCTCT','CAGAGAGG-AAGGAGTA','CTCTCTAC-AGAGTAGA','CTCTCTAC-GTAAGGAG','PB34','PB26','PB29','CGTACTAG-TATCCTCT','PB57','PB51','AGGCAGAA-GTAAGGAG','PB90','PB36','PB82','PB20','PB65','PB60','PB21','PB58','PB43','PB66','GTAGAGGA-AAGGAGTA','PB12','AGGCAGAA-AAGGAGTA','GTAGAGGA-CTCTCTAT','PB75','AGGCAGAA-AGAGTAGA','PB41','CGTACTAG-TAGATCGC','PB25','PB_2']

# Find all fixed 9 clade loci
fixed = find_fixed_loci(ladder, data_index, data_positions)


mask = [True if x in fixed else False for x in data_positions]

counts = []
for s in strains:
    vector = np.array(data.loc[s][mask])
    counts.append(sum(vector==1))

# Make a bar graph based on tree order
plt.bar(strains, counts)
plt.title('How persistent are the ladder fixed loci  in the whole PSB tree?')
plt.xlabel('Strain names (tree order)')
plt.ylabel('counts of fixed ladder loci')
plt.xticks(rotation=90)
plt.show()