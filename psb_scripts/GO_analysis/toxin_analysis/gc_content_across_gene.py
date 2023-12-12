#########################################################
# Given a location of a gene on the chromosome, check 
# the GC content at those positions.
# We have been looking at strict biallelic loci but we want 
# to now explore the rest of the data for those specific genes.
#########################################################


from collections import Counter 
from pysam import VariantFile
import sys
import pandas as pd
from pysam import VariantFile
import numpy as np


# Select which region of the chromosome we will look at.
start, end = 3362633, 3366313  # wapA_2
#start, end = 5049505, 5058900  # wapA_3

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
bc_strains = ['TAAGGCGA-TATCCTCT','PB59','PB45','PB77','GGACTCCT-TAGATCGC','GGACTCCT-ACTGCATA','PB78','PB61','PB13','PB85','PB37','PB_5','AGGCAGAA-ACTGCATA','AGGCAGAA-CTAAGCCT','PB25','PB67','PB76','PB59','PB84','PB52','PB44']
e_clade = ['PB11','CTCTCTAC-AAGGAGTA','GTAGAGGA-CTAAGCCT','CTCTCTAC-CTAAGCCT','CAGAGAGG-CTCTCTAT','CTCTCTAC-CTCTCTAT','CTCTCTAC-TAGATCGC','CTCTCTAC-ACTGCATA','CAGAGAGG-TAGATCGC','CAGAGAGG-CTAAGCCT','CAGAGAGG-GTAAGGAG','CTCTCTAC-TATCCTCT','CAGAGAGG-ACTGCATA','CAGAGAGG-AGAGTAGA','CAGAGAGG-TATCCTCT','CAGAGAGG-AAGGAGTA','CTCTCTAC-AGAGTAGA','CTCTCTAC-GTAAGGAG','PB34','PB26']
basal_clade = ['PB29','CGTACTAG-TATCCTCT','PB57','PB51','AGGCAGAA-GTAAGGAG','PB90','PB36','PB82','PB20','PB65','PB60','PB21','PB58','PB43','PB66','GTAGAGGA-AAGGAGTA','PB12','AGGCAGAA-AAGGAGTA','GTAGAGGA-CTCTCTAT','PB75','AGGCAGAA-AGAGTAGA','PB41','CGTACTAG-TAGATCGC','PB25','PB_2']
strains = ['PB87','PB40','TCCTGAGC-TAGATCGC','TCCTGAGC-CTCTCTAT','PB73','AAGAGGCA-TATCCTCT','GGACTCCT-AGAGTAGA','GGACTCCT-CTAAGCCT','GGACTCCT-CTAAGCCT','PB80','PB63','AAGAGGCA-TAGATCGC','PB39','PB_8','PB24','PB64','PB55','AAGAGGCA-AGAGTAGA','GCTACGCT-GTAAGGAG','AAGAGGCA-ACTGCATA','AAGAGGCA-AAGGAGTA','GCTACGCT-ACTGCATA','PB31','AAGAGGCA-GTAAGGAG','AAGAGGCA-CTAAGCCT','AAGAGGCA-CTCTCTAT','PB32','PB16','PB47','PB48','TAAGGCGA-TATCCTCT','PB59','PB45','PB77','GGACTCCT-TAGATCGC','GGACTCCT-ACTGCATA','PB78','PB61','PB13','PB85','PB37','PB_5','AGGCAGAA-ACTGCATA','AGGCAGAA-CTAAGCCT','PB25','PB67','PB76','PB59','PB84','PB52','PB44','PB11','CTCTCTAC-AAGGAGTA','GTAGAGGA-CTAAGCCT','CTCTCTAC-CTAAGCCT','CAGAGAGG-CTCTCTAT','CTCTCTAC-CTCTCTAT','CTCTCTAC-TAGATCGC','CTCTCTAC-ACTGCATA','CAGAGAGG-TAGATCGC','CAGAGAGG-CTAAGCCT','CAGAGAGG-GTAAGGAG','CTCTCTAC-TATCCTCT','CAGAGAGG-ACTGCATA','CAGAGAGG-AGAGTAGA','CAGAGAGG-TATCCTCT','CAGAGAGG-AAGGAGTA','CTCTCTAC-AGAGTAGA','CTCTCTAC-GTAAGGAG','PB34','PB26','PB29','CGTACTAG-TATCCTCT','PB57','PB51','AGGCAGAA-GTAAGGAG','PB90','PB36','PB82','PB20','PB65','PB60','PB21','PB58','PB43','PB66','GTAGAGGA-AAGGAGTA','PB12','AGGCAGAA-AAGGAGTA','GTAGAGGA-CTCTCTAT','PB75','AGGCAGAA-AGAGTAGA','PB41','CGTACTAG-TAGATCGC','PB25','PB_2']

#good_samples = e_clade + ladder + bc_strains + basal_clade
good_samples = e_clade

positions = []
coverage = []
counter = 0
alternative_alleles = []

for pos in bcf_iter:
    counter += 1
    print(pos.pos)
    
    if pos.pos >= start and pos.pos <= end:

        # gather coverage as alleles, mark non-calls, complex alleles and coverage <=3 as 'N'
        tmp_coverage = [np.nan if pos.samples[s]['DP'] is None else pos.samples[s]['DP'] for s in good_samples]
        tmp_coverage = np.array(tmp_coverage)

        tmp_ref = ['N' if pos.ref is None or tmp_coverage[si] < 3
                    else pos.ref for si,s in enumerate(good_samples)]
        
        tmp_ref_count = [np.nan if pos.samples[s]['RO'] is None else pos.samples[s]['RO'] for s in good_samples]
        tmp_ref_count = np.array(tmp_ref_count)

        tmp_allele_count = [np.nan if pos.samples[s]['AO'] is None else pos.samples[s]['AO'] for s in good_samples]
        #tmp_allele_count = np.array(tmp_allele_count)

        tmp_alleles = ['N' if pos.samples[s].alleles[0] is None or tmp_coverage[si] < 3
                    else pos.samples[s].alleles for si,s in enumerate(good_samples)]
        tmp_alleles = np.array(tmp_alleles)

        coverage.append(np.nanmean(tmp_coverage))
        positions.append(str(pos.pos))
        pos_data = zip(tmp_ref, tmp_ref_count, tmp_allele_count, tmp_alleles)
        print(list(pos_data))
        #print(tmp_alleles)
        alternative_alleles.append(tmp_alleles)
        
    elif pos.pos > end:
        break

# Create a plot where we calculate the frequency of different alleles at each position
# Start with Ns
n_frequency, a_frequency, t_frequency, c_frequency, g_frequency, other_frequency = [], [], [], [], [], []
for i in range(0, len(alternative_alleles)):
    locus = alternative_alleles[i]
    frequency = Counter(locus)
    #print(frequency)
    try:
        n_freq = frequency['N']/len(good_samples)
        n_frequency.append(str(n_freq))
    except KeyError:
        n_frequency.append('0')

    try:
        a_freq = frequency['A']/len(good_samples) 
        a_frequency.append(str(a_freq))
    except KeyError:
        a_frequency.append('0')

    try:
        t_freq = frequency['T']/len(good_samples) 
        t_frequency.append(str(t_freq))
    except KeyError:
        t_frequency.append('0')

    try:   
        c_freq = frequency['C']/len(good_samples)
        c_frequency.append(str(c_freq))
    except KeyError:
        c_frequency.append('0')

    try:
        g_freq = frequency['G']/len(good_samples) 
        g_frequency.append(str(g_freq))
    except KeyError:
        g_frequency.append('0')

    # Check if there are any changes that are not single nucleotides (insertions / deletions)
    other_counter = 0
    for k in frequency.keys():
        if k not in ['A', 'C', 'T', 'G', 'N']:
            other_counter += 1
    
    other_frequency.append(str(other_counter/len(good_samples)))
"""
out_file = open('/home/ada/Desktop/PinkBerry_scripts_paper/data/psb/wapA_3_allele_frequency_all_clades.txt', 'w')
out_file.write(','.join(positions) + '\n')
out_file.write(','.join(n_frequency) + '\n')
out_file.write(','.join(a_frequency) + '\n')
out_file.write(','.join(t_frequency) + '\n')
out_file.write(','.join(c_frequency) + '\n')
out_file.write(','.join(g_frequency) + '\n')
out_file.write(','.join(other_frequency))
out_file.close()
"""
