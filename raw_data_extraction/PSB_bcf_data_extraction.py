#############################################################
# Extract SNPs from the original BCF file.
# Look only for SNPs with one major and one minor allele (biallelic).
# Don't ignore missing data. If missing data is present (less than 4 reads)
# include that info in the SNP matrix.
#############################################################

from pysam import VariantFile
import numpy as np
import sys
import pandas as pd

np.set_printoptions(threshold=sys.maxsize)

# Read in the original BCF file.
bcf_in = VariantFile("psb-scaff03.minimap2-freebayes.bcf")  # auto-detect input format
bcf_iter = bcf_in.fetch()
all_samples = list(bcf_in.header.samples)

# good_samples = [all_samples[i] for i in worst_to_best[50:]]
# Samples that have most loci covered. (Ignore the worst 50 strains)
good_samples = ['PB63', 'PB28', 'GGACTCCT-AGAGTAGA', 'AAGAGGCA-AGAGTAGA', 'PB73', 'PB76', 'AAGAGGCA-ACTGCATA', 'GGACTCCT-GTAAGGAG', 'GCTACGCT-GTAAGGAG', 'TCCTGAGC-CTCTCTAT', 'PB55', 'PB_5', 'PB31', 'AAGAGGCA-TATCCTCT', 'PB24', 'PB87', 'AAGAGGCA-AAGGAGTA', 'PB84', 'PB40', 'AGGCAGAA-ACTGCATA', 'AGGCAGAA-CTAAGCCT', 'PB67', 'PB_8', 'PB80', 'PB39', 'PB47', 'PB59', 'PB37', 'PB61', 'TCCTGAGC-TAGATCGC', 'PB77', 'PB69', 'TAAGGCGA-AAGGAGTA', 'PB13', 'PB82', 'PB48', 'PB52', 'GGACTCCT-TAGATCGC', 'PB11', 'GGACTCCT-CTAAGCCT', 'PB44', 'PB36', 'AAGAGGCA-GTAAGGAG', 'PB16', 'PB32', 'PB20', 'GCTACGCT-ACTGCATA', 'PB85', 'PB65', 'PB58', 'PB64', 'CTCTCTAC-GTAAGGAG', 'PB60', 'TAAGGCGA-GTAAGGAG', 'AGGCAGAA-GTAAGGAG', 'PB21', 'TAAGGCGA-TATCCTCT', 'AAGAGGCA-CTCTCTAT', 'PB78', 'PB_1', 'PB91', 'AAGAGGCA-TAGATCGC', 'CTCTCTAC-AAGGAGTA', 'GTAGAGGA-ACTGCATA', 'CAGAGAGG-AAGGAGTA', 'CTCTCTAC-AGAGTAGA', 'GTAGAGGA-GTAAGGAG', 'PB43', 'PB34', 'AAGAGGCA-CTAAGCCT', 'PB92', 'CTCTCTAC-TATCCTCT', 'AGGCAGAA-AAGGAGTA', 'PB90', 'PB12', 'CAGAGAGG-TATCCTCT', 'PB18', 'PB26', 'GTAGAGGA-AAGGAGTA', 'PB41', 'AGGCAGAA-AGAGTAGA', 'PB75', 'PB45', 'CAGAGAGG-TAGATCGC', 'TAAGGCGA-AGAGTAGA','CAGAGAGG-GTAAGGAG','PB17','CAGAGAGG-ACTGCATA','CTCTCTAC-CTCTCTAT','PB51','CGTACTAG-AGAGTAGA','PB66','PB25','PB83','CGTACTAG-CTAAGCCT','CTCTCTAC-TAGATCGC','TCCTGAGC-AAGGAGTA','CAGAGAGG-AGAGTAGA','PB_3','AGGCAGAA-CTCTCTAT','CAGAGAGG-CTAAGCCT','TAAGGCGA-TAGATCGC','PB81','CGTACTAG-ACTGCATA','CGTACTAG-AAGGAGTA','PB89','TAAGGCGA-ACTGCATA','PB74','CTCTCTAC-ACTGCATA','PB57','PB49','GTAGAGGA-CTAAGCCT','CGTACTAG-CTCTCTAT','GTAGAGGA-CTCTCTAT','PB_9','PB_4','GTAGAGGA-AGAGTAGA','GTAGAGGA-TATCCTCT','PB29','PB10','PB_2','AGGCAGAA-TAGATCGC','CGTACTAG-TATCCTCT','PB53','GGACTCCT-ACTGCATA','PB19','CTCTCTAC-CTAAGCCT','TCCTGAGC-AGAGTAGA','CAGAGAGG-CTCTCTAT','PB35','PB27','GTAGAGGA-TAGATCGC','CGTACTAG-TAGATCGC','TCCTGAGC-CTAAGCCT','TAAGGCGA-CTAAGCCT','TCCTGAGC-TATCCTCT','PB42','TCCTGAGC-GTAAGGAG','TAAGGCGA-CTCTCTAT','TCCTGAGC-ACTGCATA','AGGCAGAA-TATCCTCT','CGTACTAG-GTAAGGAG']


positions = []
coverage = []
called = []
alleles = []
meta = []
rec = 0
counter = 0
# max 100000 SNPs
rec_max=100000
for pos in bcf_iter:
    counter += 1
    print(counter)
    # gather coverage as alleles, mark non-calls, complex alleles and coverage <=3 as 'N'
    tmp_coverage = [-1 if pos.samples[s]['DP'] is None else pos.samples[s]['DP'] for s in good_samples]
    # Get a major allele for each bacteria
    tmp_alleles = ['N' if pos.samples[s].alleles[0] is None or len(pos.samples[s].alleles[0])>1 or tmp_coverage[si]<3
                   else pos.samples[s].alleles[0] for si,s in enumerate(good_samples)]
    tmp_meta = [[tmp_coverage[si], pos.samples[s]['RO'], pos.samples[s]['AO']] for si, s in enumerate(good_samples)]

    # Skip if only one proper call or not a proper biallelic.
    unique_alleles = set(tmp_alleles)
    if (len(unique_alleles)==2 and 'N' not in unique_alleles) or (len(unique_alleles)==3 and 'N' in unique_alleles):
        alleles.append(tmp_alleles)
        positions.append(pos.pos)
        coverage.append(tmp_coverage)
        called.append([pos.samples[s].alleles[0] is not None for s in good_samples]) # Not a blank space
        meta.append(tmp_meta)
        rec+=1
        if rec>rec_max:
            break
    else:
        continue

# Convert to array.
called = np.array(called).T
coverage = np.array(coverage).T
alleles = np.array(alleles).T


# Make a masked array of SNPs to keep NaNs.
snps = np.ma.array(alleles, fill_value='N')
snps.mask = (snps=='-')|(snps=='N')  # We want Ns and gaps to stay labeled 

# Get allele counts and major alleles
alpha = np.array(['A', 'C', 'G', 'T'])
af = np.array([np.sum(snps==a, axis=0) for a in 'ACGT'])
maj = alpha[af.argmax(axis=0)]

print(maj)


# SNPs are not major alleles > SNPs = 1, maj = 0
maj_snps = snps!=maj
maj_snps = maj_snps.astype(int)

# Create dataframes to make saving easier.
df = pd.DataFrame(maj_snps, index=good_samples, columns=positions)
mdf = pd.DataFrame(meta, index=positions, columns=good_samples)
mdf = mdf.T


# Save the matrices to CSV
df.to_csv('snp_data_2020.csv')
mdf.to_csv('snp_meta_data_2020.csv')
