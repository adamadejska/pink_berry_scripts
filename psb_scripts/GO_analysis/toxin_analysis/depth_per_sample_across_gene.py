#########################################################
# Given a location of a gene on the chromosome, check how
# many alternative alleles are there at each position.
# We have been looking at strict biallelic loci but we want 
# to now explore the rest of the data for those specific genes.
#########################################################


from pysam import VariantFile
import sys
import pandas as pd
from pysam import VariantFile
import numpy as np


# Select which region of the chromosome we will look at.
#start, end = 3362633, 3366313  # wapA_2
#start, end = 5049505, 5058900  # wapA_3
#start, end = 3362633, 3375078  # wapA_2 and following ORFs
#start, end = 5037116, 5058900  # wapA_3 and following ORFs
#start, end = 6430294, 6435702  # wapA_6
#start, end = 6427067, 6435702  # wapA_6 and following ORFs
#start, end = 3350557, 3354828  # rhsC_3
#start, end = 3350557, 3359484  # rhsC_3 and following ORFs
#start, end = 4219239, 4234403  # rhsC_4
start, end = 4213680, 4234403   # rhsC_4 and following ORFs

# Make a simple gene map without figuring out the kinds of mutations
# Read in the CSV file of the variance data.
data = pd.read_csv('/home/ada/Desktop/Shraiman_lab/data/snp_data_2020.csv', index_col=0)

# Index = names of bacterial samples
index = data.index.values.tolist()
positions = data.columns.tolist()

# calculate the coverage for each locus in the range
np.set_printoptions(threshold=sys.maxsize)
bcf_in = VariantFile("/home/ada/Desktop/Shraiman_lab/data/psb-scaff03.minimap2-freebayes.bcf")  # auto-detect input format
bcf_iter = bcf_in.fetch()

# samples that have most loci covered
clade_9_strains = ['PB87','PB40','TCCTGAGC-TAGATCGC','TCCTGAGC-CTCTCTAT','PB73','AAGAGGCA-TATCCTCT','GGACTCCT-AGAGTAGA','GGACTCCT-CTAAGCCT']
ladder = ['GGACTCCT-CTAAGCCT','PB80','PB63','AAGAGGCA-TAGATCGC','PB39','PB_8','PB24','PB64','PB55','AAGAGGCA-AGAGTAGA','GCTACGCT-GTAAGGAG','AAGAGGCA-ACTGCATA','AAGAGGCA-AAGGAGTA','GCTACGCT-ACTGCATA','PB31','AAGAGGCA-GTAAGGAG','AAGAGGCA-CTAAGCCT','AAGAGGCA-CTCTCTAT','PB32','PB16','PB47','PB48']
bc_strains = ['AGGCAGAA-ACTGCATA','AGGCAGAA-CTAAGCCT','TAAGGCGA-TATCCTCT','PB78','PB61','GGACTCCT-TAGATCGC','PB45','PB77','PB69','PB_5','PB37','PB13','PB85','PB28','PB67','PB59','PB52','PB26','PB34','PB11']
e_clade = ['PB11','CTCTCTAC-AAGGAGTA','GTAGAGGA-CTAAGCCT','CTCTCTAC-CTAAGCCT','CAGAGAGG-CTCTCTAT','CTCTCTAC-CTCTCTAT','CTCTCTAC-TAGATCGC','CTCTCTAC-ACTGCATA','CAGAGAGG-TAGATCGC','CAGAGAGG-CTAAGCCT','CAGAGAGG-GTAAGGAG','CTCTCTAC-TATCCTCT','CAGAGAGG-ACTGCATA','CAGAGAGG-AGAGTAGA','CAGAGAGG-TATCCTCT','CAGAGAGG-AAGGAGTA','CTCTCTAC-AGAGTAGA','CTCTCTAC-GTAAGGAG','PB34','PB26']
basal_clade = ['PB29','CGTACTAG-TATCCTCT','PB57','PB51','AGGCAGAA-GTAAGGAG','PB90','PB36','PB82','PB20','PB65','PB60','PB21','PB58','PB43','PB66','GTAGAGGA-AAGGAGTA','PB12','AGGCAGAA-AAGGAGTA','GTAGAGGA-CTCTCTAT','PB75','AGGCAGAA-AGAGTAGA','PB41','CGTACTAG-TAGATCGC','PB25','PB_2']
strains = ['PB87','PB40','TCCTGAGC-TAGATCGC','TCCTGAGC-CTCTCTAT','PB73','AAGAGGCA-TATCCTCT','GGACTCCT-AGAGTAGA','GGACTCCT-CTAAGCCT','GGACTCCT-CTAAGCCT','PB80','PB63','AAGAGGCA-TAGATCGC','PB39','PB_8','PB24','PB64','PB55','AAGAGGCA-AGAGTAGA','GCTACGCT-GTAAGGAG','AAGAGGCA-ACTGCATA','AAGAGGCA-AAGGAGTA','GCTACGCT-ACTGCATA','PB31','AAGAGGCA-GTAAGGAG','AAGAGGCA-CTAAGCCT','AAGAGGCA-CTCTCTAT','PB32','PB16','PB47','PB48','TAAGGCGA-TATCCTCT','PB59','PB45','PB77','GGACTCCT-TAGATCGC','PB78','PB61','PB13','PB85','PB37','PB_5','AGGCAGAA-ACTGCATA','AGGCAGAA-CTAAGCCT','PB25','PB67','PB76','PB59','PB84','PB52','PB44','PB11','CTCTCTAC-AAGGAGTA','GTAGAGGA-CTAAGCCT','CTCTCTAC-CTAAGCCT','CAGAGAGG-CTCTCTAT','CTCTCTAC-CTCTCTAT','CTCTCTAC-TAGATCGC','CTCTCTAC-ACTGCATA','CAGAGAGG-TAGATCGC','CAGAGAGG-CTAAGCCT','CAGAGAGG-GTAAGGAG','CTCTCTAC-TATCCTCT','CAGAGAGG-ACTGCATA','CAGAGAGG-AGAGTAGA','CAGAGAGG-TATCCTCT','CAGAGAGG-AAGGAGTA','CTCTCTAC-AGAGTAGA','CTCTCTAC-GTAAGGAG','PB34','PB26','PB29','CGTACTAG-TATCCTCT','PB57','PB51','AGGCAGAA-GTAAGGAG','PB90','PB36','PB82','PB20','PB65','PB60','PB21','PB58','PB43','PB66','GTAGAGGA-AAGGAGTA','PB12','AGGCAGAA-AAGGAGTA','GTAGAGGA-CTCTCTAT','PB75','AGGCAGAA-AGAGTAGA','PB41','CGTACTAG-TAGATCGC','PB25','PB_2']

# Get rid of 'GGACTCCT-ACTGCATA' from the bc clade (outlier in depth)
#good_samples =  e_clade + ladder  + basal_clade + clade_9_strains + bc_strains
good_samples =  clade_9_strains

positions = []
coverage = []
alternative_alleles = []

for pos in bcf_iter:
    print(pos.pos)
    
    # gather coverage as alleles, mark non-calls, complex alleles and coverage <=3 as 'N'
    tmp_coverage = ['0' if pos.samples[s]['DP'] is None else str(pos.samples[s]['DP']) for s in good_samples]
    tmp_coverage = np.array(tmp_coverage)

    if pos.pos >= start and pos.pos <= end:
        coverage.append(tmp_coverage)
        positions.append(str(pos.pos))

    elif pos.pos > end:
        break

coverage = np.array(coverage)
coverage = coverage.transpose()
#out_file = open('/home/ada/Desktop/PinkBerry_scripts_paper/data/psb/rhsC_4_allele_depths_8_clade.txt', 'w')
out_file = open('/home/ada/Desktop/PinkBerry_scripts_paper/data/psb/rhsC_4_and_further_allele_depths_8_clade.txt', 'w')
#out_file = open('/home/ada/Desktop/PinkBerry_scripts_paper/data/psb/wapA_2_and_further_allele_depths_E_clade.txt', 'w')
out_file.write(','.join(positions) + '\n')

for i in range(0, len(coverage)):
    out_file.write(','.join(coverage[i]) + '\n')

out_file.close()

