
from pysam import VariantFile
import numpy as np
import sys
import pandas as pd
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
meta = []
rec = 0
chromosome = []
# max 100000 SNPs
max_count=850
counter = 0
for pos in bcf_iter:
    counter += 1
    if counter > max_count:
        break
    # gather coverage as alleles, mark non-calls, complex alleles and coverage <=3 as 'N'
    tmp_coverage = [-1 if pos.samples[s]['DP'] is None else pos.samples[s]['DP'] for s in good_samples]
    tmp_alleles = ['N' if pos.samples[s].alleles[0] is None else pos.samples[s].alleles[0] for si,s in enumerate(good_samples)]
    if 'N' in tmp_alleles:
        print(pos.pos)
        print([pos.samples[s].alleles[0] for si,s in enumerate(good_samples)])
        print([s for si,s in enumerate(good_samples)])
