######################################################
# Using the original bcf file, create a dataset of original
# positions and changes for the SNP matrix positions. 
######################################################

from collections import Counter
from pysam import VariantFile
import numpy as np
import sys
import pandas as pd
import operator
import matplotlib.pyplot as plt

np.set_printoptions(threshold=sys.maxsize)

bcf_in = VariantFile("/home/ada/Desktop/Shraiman_lab/srb_data/srb.minimap2-freebayes.bcf")  # auto-detect input format
bcf_iter = bcf_in.fetch()
all_samples = list(bcf_in.header.samples)

# samples that have most loci covered
good_samples = ['AAGAGGCA-AAGGAGTA', 'PB36', 'PB80', 'PB23', 'PB78', 'TCCTGAGC-AAGGAGTA', 'PB47', 'PB48', 'AAGAGGCA-ACTGCATA', 'PB60', 'PB_1', 'PB91', 'PB95', 'PB19', 'PB59', 'GTAGAGGA-ACTGCATA', 'PB81', 'PB28', 'PB68', 'GTAGAGGA-GTAAGGAG', 'PB73', 'PB93', 'GGACTCCT-GTAAGGAG', 'PB39', 'TAAGGCGA-AGAGTAGA', 'GTAGAGGA-CTAAGCCT', 'PB32', 'PB33', 'PB79', 'PB_9', 'PB64', 'PB66', 'CAGAGAGG-ACTGCATA', 'GTAGAGGA-CTCTCTAT', 'CAGAGAGG-TATCCTCT', 'CTCTCTAC-TATCCTCT', 'AAGAGGCA-CTAAGCCT', 'AGGCAGAA-AGAGTAGA', 'AAGAGGCA-TATCCTCT', 'CAGAGAGG-AAGGAGTA', 'PB53', 'CTCTCTAC-TAGATCGC', 'GCTACGCT-ACTGCATA', 'CAGAGAGG-TAGATCGC', 'GGACTCCT-TAGATCGC', 'CTCTCTAC-ACTGCATA', 'CTCTCTAC-AGAGTAGA', 'PB37', 'PB16', 'GTAGAGGA-TATCCTCT', 'CAGAGAGG-AGAGTAGA', 'AAGAGGCA-GTAAGGAG', 'CAGAGAGG-GTAAGGAG', 'PB77', 'PB29', 'CTCTCTAC-CTCTCTAT', 'PB83', 'TAAGGCGA-AAGGAGTA', 'PB_5', 'PB43', 'TAAGGCGA-ACTGCATA', 'PB40', 'PB45', 'CAGAGAGG-CTAAGCCT', 'AGGCAGAA-CTCTCTAT', 'PB34', 'CTCTCTAC-CTAAGCCT', 'PB85', 'AAGAGGCA-TAGATCGC', 'TCCTGAGC-CTAAGCCT', 'AAGAGGCA-CTCTCTAT', 'PB87', 'CGTACTAG-AGAGTAGA', 'PB25', 'PB10', 'TAAGGCGA-TATCCTCT', 'PB57', 'PB_4', 'TAAGGCGA-CTCTCTAT', 'GTAGAGGA-AGAGTAGA', 'PB88', 'GGACTCCT-ACTGCATA', 'PB18', 'CGTACTAG-AAGGAGTA', 'AGGCAGAA-GTAAGGAG', 'TCCTGAGC-TATCCTCT', 'PB21', 'CGTACTAG-ACTGCATA', 'PB41', 'CAGAGAGG-CTCTCTAT', 'PB76', 'AGGCAGAA-TAGATCGC', 'PB_3', 'TCCTGAGC-AGAGTAGA', 'PB11', 'PB27', 'TAAGGCGA-CTAAGCCT', 'PB44', 'PB84', 'AGGCAGAA-AAGGAGTA', 'PB90', 'CGTACTAG-CTAAGCCT', 'GTAGAGGA-TAGATCGC', 'PB35', 'CGTACTAG-TATCCTCT', 'TCCTGAGC-GTAAGGAG', 'PB82', 'CGTACTAG-TAGATCGC', 'PB50', 'GGACTCCT-CTAAGCCT', 'CGTACTAG-GTAAGGAG', 'AGGCAGAA-TATCCTCT', 'CGTACTAG-CTCTCTAT', 'TCCTGAGC-ACTGCATA', 'PB58', 'PB42']

positions = []
coverage = []
called = []
alleles = []
counter = 0
chromosome = []
# max 100000 SNPs
rec_max=100000

f = open('/home/ada/Desktop/Shraiman_lab/srb_data/wt_vs_alt_srb.csv', 'w')
f.write('pos,chr,wt,alt\n')
for pos in bcf_iter:
    counter += 1
    print(counter)
    # gather coverage as alleles, mark non-calls, complex alleles and coverage <=3 as 'N'
    tmp_coverage = [-1 if pos.samples[s]['DP'] is None else pos.samples[s]['DP'] for s in good_samples]
    # Get a major allele for each bacteria
    tmp_alleles = ['N' if pos.samples[s].alleles[0] is None or len(pos.samples[s].alleles[0])>1 or tmp_coverage[si]<3
                   else pos.samples[s].alleles[0] for si,s in enumerate(good_samples)]

    # skip if only one proper call or not a proper biallelic
    unique_alleles = set(tmp_alleles)
    if (len(unique_alleles)==2 and 'N' not in unique_alleles) or (len(unique_alleles)==3 and 'N' in unique_alleles):       
        f.write(str(pos.pos) + ',' + str(pos.chrom.split('_')[1]) + ','+ pos.alleles[0] + ',' + pos.alleles[1] + '\n')
    else:
        continue

f.close()