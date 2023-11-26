##################################################################
# Continue the analysis of sample PB93 which is a hypermutator. 
# Is this sample well sequenced? What's the coverage for it?
# How does it compare to the average coverage of the whole dataset?
##################################################################

from pysam import VariantFile
import numpy as np
import sys
import pandas as pd
import matplotlib.pyplot as plt


np.set_printoptions(threshold=sys.maxsize)

# Read in the SNP matrix to get the relevant sample names. 
data = pd.read_csv('/home/ada/Desktop/Shraiman_lab/srb_data/srb_snp_data_2021_chr.csv', index_col=0)

# Index = names of bacterial samples
data_index = data.index.values.tolist()[1:]

# Read in the original BCF file.
bcf_in = VariantFile('/home/ada/Desktop/Shraiman_lab/srb_data/srb.minimap2-freebayes.bcf')  # auto-detect input format
bcf_iter = bcf_in.fetch()
all_samples = list(bcf_in.header.samples)

# Samples that have most loci covered. (Ignore the worst 50 strains)
good_samples = data_index
#sample = 'PB93'

coverage = []
positions = []
rec = 0
counter = 0
rec_max=20000
chr = '0'

# Read through the BCF file to extract the coverage data.
for pos in bcf_iter:
    counter += 1
    print(counter)
    # gather coverage as alleles, mark non-calls, complex alleles and coverage <=3 as 'N'
    tmp_coverage = [-1 if pos.samples[s]['DP'] is None else pos.samples[s]['DP'] for s in good_samples]
    # Get a major allele for each bacteria
    tmp_alleles = ['N' if pos.samples[s].alleles[0] is None or len(pos.samples[s].alleles[0])>1 or tmp_coverage[si]<3
                   else pos.samples[s].alleles[0] for si,s in enumerate(good_samples)]
    
    # Skip if only one proper call or not a proper biallelic.
    unique_alleles = set(tmp_alleles)
    if (len(unique_alleles)==2 and 'N' not in unique_alleles) or (len(unique_alleles)==3 and 'N' in unique_alleles):
        if pos.chrom.split('_')[1] == '0':
            coverage.append(tmp_coverage)
            positions.append(pos.pos)
        else:
            break
    else:
        continue

# Coverage is a matrix of positions x strains. Find which index is the PB93 to extract its coverage values
index = good_samples.index('PB93')
pb93_coverage = [coverage[i][index] for i in range(0,len(coverage))]
#print(pb93_coverage)
for i in good_samples:
    index = good_samples.index(i)
    line = [coverage[i][index] for i in range(0,len(coverage))]
    plt.plot(positions, line, color='gray', alpha=0.5)
plt.plot(positions, pb93_coverage, color='red')
plt.xlabel('contig 0 position')
plt.ylabel('coverage')
plt.title('Coverage of all SRB samples for contig 0 vs coverage of PB93')
plt.show()


# Calculate average coverage of each position across strains
avg_coverage = []
for pos in range(0, len(pb93_coverage)):
    tmp_cov_list = []
    for strain in range(0, len(good_samples)):
        tmp_cov_list.append(coverage[pos][strain])

    avg_coverage.append(np.mean(np.array(tmp_cov_list)))

plt.plot(positions, avg_coverage, color='gray')
plt.plot(positions, pb93_coverage, color='red')
plt.xlabel('contig 0 position')
plt.ylabel('coverage')
plt.title('Average coverage of all SRB samples for contig 0 vs coverage of PB93')
plt.show()


# Calculate average of average coverage across a window size.
window = 100.0
all_window_coverage = {}
for i in range(0, len(positions)):
    current_window = positions[i]//window
    if current_window not in all_window_coverage.keys():
        all_window_coverage[current_window] = [avg_coverage[i]]
    else:
        all_window_coverage[current_window].append(avg_coverage[i])

# For every key in the dictionary, calculate the average.
avg_window_cov = []
all_window_list = []
for k,v in all_window_coverage.items():
    avg_window_cov.append(np.mean(np.array(v)))
    all_window_list.append(k)

# Do the same for PB93 coverage.

pb93_window_coverage = {}
for i in range(0, len(positions)):
    current_window = positions[i]//window
    if current_window not in pb93_window_coverage.keys():
        pb93_window_coverage[current_window] = [pb93_coverage[i]]
    else:
        pb93_window_coverage[current_window].append(pb93_coverage[i])

pb93_avg_window_cov = []
pb93_window_list = []
for k,v in pb93_window_coverage.items():
    pb93_avg_window_cov.append(np.mean(np.array(v)))
    pb93_window_list.append(k)

plt.plot(all_window_list, avg_window_cov, color='gray')
plt.plot(pb93_window_list, pb93_avg_window_cov, color='red')
plt.xlabel('contig 0 position (window size=100pb)')
plt.ylabel('coverage')
plt.title('Average coverage of all SRB samples for contig 0 vs coverage of PB93\nwindow size = 100 bp')
plt.show()