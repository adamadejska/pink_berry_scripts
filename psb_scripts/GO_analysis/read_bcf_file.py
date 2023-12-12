# Read bcf file

from pysam import VariantFile
import sys
import matplotlib.pyplot as plt
import pandas as pd
from pysam import VariantFile
import numpy as np

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
good_samples = index

positions = []
coverage = []
counter = 0
start, end = 3362633, 3375078

for pos in bcf_iter:
    counter += 1
    # gather coverage as alleles, mark non-calls, complex alleles and coverage <=3 as 'N'
    tmp_coverage = [np.nan if pos.samples[s]['DP'] is None else pos.samples[s]['DP'] for s in good_samples]
    tmp_coverage = np.array(tmp_coverage)

    tmp_ref_count = [np.nan if pos.samples[s]['RO'] is None else pos.samples[s]['RO'] for s in good_samples]
    tmp_ref_count = np.array(tmp_ref_count)

    tmp_allele_count = [np.nan if pos.samples[s]['AO'] is None else pos.samples[s]['AO'] for s in good_samples]
    tmp_allele_count = np.array(tmp_allele_count)

    tmp_alleles = ['N' if pos.samples[s].alleles[0] is None or tmp_coverage[si] < 3
                   else pos.samples[s].alleles[0] for si,s in enumerate(good_samples)]

    if pos.pos >= start and pos.pos <= end:
        coverage.append(np.nanmean(tmp_coverage))
        positions.append(pos.pos)
        pos_data = zip(tmp_ref_count, tmp_allele_count, tmp_alleles)
        print(list(pos_data))
    elif pos.pos > end:
        break