#########################################################
# Given a location of a gene on the chromosome, check how
# many alternative alleles are there at each position.
# We have been looking at strict biallelic loci but we want 
# to now explore the rest of the data for those specific genes.
#########################################################


from collections import Counter 
from pysam import VariantFile
import sys
import matplotlib.pyplot as plt
import pandas as pd
from pysam import VariantFile
import numpy as np


# Select which region of the chromosome we will look at.
start, end = 3362633, 3366313

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

good_samples = clade_9_strains + ladder + bc_strains + basal_clade

positions = []
coverage = []
counter = 0
alternative_alleles = []

for pos in bcf_iter:
    counter += 1
    print(pos.pos)
    
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
                   else pos.samples[s].alleles[0] for si,s in enumerate(good_samples)]
    tmp_alleles = np.array(tmp_alleles)
    
    if pos.pos >= start and pos.pos <= end:
        coverage.append(np.nanmean(tmp_coverage))
        positions.append(str(pos.pos))
        pos_data = zip(tmp_ref, tmp_ref_count, tmp_allele_count, tmp_alleles)
        #print(list(pos_data))

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

out_file = open('/home/ada/Desktop/PinkBerry_scripts_paper/data/psb/wapA_2_allele_frequency_all_clades.txt', 'w')
out_file.write(','.join(positions) + '\n')
out_file.write(','.join(n_frequency) + '\n')
out_file.write(','.join(a_frequency) + '\n')
out_file.write(','.join(t_frequency) + '\n')
out_file.write(','.join(c_frequency) + '\n')
out_file.write(','.join(g_frequency) + '\n')
out_file.write(','.join(other_frequency))
out_file.close()

"""
def divide_into_three_parts(freq_list, l):
    """
    #Divide a given list into three parts (for three graphs)
"""
    freq_list_one = np.array(freq_list[:int(l/3)])
    freq_list_two = np.array(freq_list[int(l/3):int(l/3)*2])
    freq_list_three = np.array(freq_list[int(l/3)*2:])

    return(freq_list_one, freq_list_two, freq_list_three)


l = len(positions)
n_frequency_top, n_frequency_middle, n_frequency_bottom = divide_into_three_parts(n_frequency, l)
a_frequency_top, a_frequency_middle, a_frequency_bottom = divide_into_three_parts(a_frequency, l)
t_frequency_top, t_frequency_middle, t_frequency_bottom = divide_into_three_parts(t_frequency, l)
c_frequency_top, c_frequency_middle, c_frequency_bottom = divide_into_three_parts(c_frequency, l)
g_frequency_top, g_frequency_middle, g_frequency_bottom = divide_into_three_parts(g_frequency, l)
other_frequency_top, other_frequency_middle, other_frequency_bottom = divide_into_three_parts(other_frequency, l)

positions_top, positions_middle, positions_bottom = divide_into_three_parts(positions, l)


# Plot
fig, (ax1, ax2, ax3) = plt.subplots(3, 1)
ax1.bar(positions_top, n_frequency_top, color='gray', label='N')
ax1.bar(positions_top, a_frequency_top, bottom=n_frequency_top, color='red', label='A')
ax1.bar(positions_top, t_frequency_top, bottom=n_frequency_top+a_frequency_top, color='green', label='T')
ax1.bar(positions_top, c_frequency_top, bottom=n_frequency_top+a_frequency_top+t_frequency_top, color='orange', label='C')
ax1.bar(positions_top, g_frequency_top, bottom=n_frequency_top+a_frequency_top+t_frequency_top+c_frequency_top, color='magenta', label='G')
ax1.bar(positions_top, other_frequency_top, bottom=n_frequency_top+a_frequency_top+t_frequency_top+c_frequency_top+g_frequency_top, color='black', label='other')

ax2.bar(positions_middle, n_frequency_middle, color='gray', label='N')
ax2.bar(positions_middle, a_frequency_middle, bottom=n_frequency_middle, color='red', label='A')
ax2.bar(positions_middle, t_frequency_middle, bottom=n_frequency_middle+a_frequency_middle, color='green', label='T')
ax2.bar(positions_middle, c_frequency_middle, bottom=n_frequency_middle+a_frequency_middle+t_frequency_middle, color='orange', label='C')
ax2.bar(positions_middle, g_frequency_middle, bottom=n_frequency_middle+a_frequency_middle+t_frequency_middle+c_frequency_middle, color='magenta', label='G')
ax1.bar(positions_middle, other_frequency_middle, bottom=n_frequency_middle+a_frequency_middle+t_frequency_middle+c_frequency_middle+g_frequency_middle, color='black', label='other')

ax3.bar(positions_bottom, n_frequency_bottom, color='gray', label='N')
ax3.bar(positions_bottom, a_frequency_bottom, bottom=n_frequency_bottom, color='red', label='A')
ax3.bar(positions_bottom, t_frequency_bottom, bottom=n_frequency_bottom+a_frequency_bottom, color='green', label='T')
ax3.bar(positions_bottom, c_frequency_bottom, bottom=n_frequency_bottom+a_frequency_bottom+t_frequency_bottom, color='orange', label='C')
ax3.bar(positions_bottom, g_frequency_bottom, bottom=n_frequency_bottom+a_frequency_bottom+t_frequency_bottom+c_frequency_bottom, color='magenta', label='G')
ax1.bar(positions_bottom, other_frequency_bottom, bottom=n_frequency_bottom+a_frequency_bottom+t_frequency_bottom+c_frequency_bottom+g_frequency_bottom, color='black', label='other')

#plt.ylabel('frequency (percentage)')
fig.text(0.06, 0.5, 'frequency (percentage)', ha='center', va='center', rotation='vertical')
plt.xlabel('chromosomal position')
fig.suptitle('WapA_2 allele distribution across the gene for the E clade')
ax1.legend()
plt.show()
"""