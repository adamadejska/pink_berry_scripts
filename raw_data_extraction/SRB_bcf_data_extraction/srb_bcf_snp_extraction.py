#######################################################
# This script extracts data about biallelic SNPs from 
# SRB binary vcf SNP file. Since SRB does not have a 
# closed chromosome, there are multiple contigs that
# we need to keep track on for linkage disequilibrium
# analysis. 
# This script outputs a csv dataframe with biallelic SNPs
# for each strain as well as the info on chromosome #
#######################################################


from pysam import VariantFile
import numpy as np
import sys
import pandas as pd
import matplotlib.pyplot as plt

np.set_printoptions(threshold=sys.maxsize)

bcf_in = VariantFile("/home/ada/Desktop/Shraiman_lab/srb_data/srb.minimap2-freebayes.bcf")  # auto-detect input format
bcf_iter = bcf_in.fetch()
all_samples = list(bcf_in.header.samples)

#worst_to_best =[110,67,101,35,10,191,121,137,105,70,8,135,41,178,109,5,179,32,160,79,134,52,161,163,106,95,83,156,127,98,0,6,190,140,36,177,144,9,34,68,29,21,87,147,112,63,171,90,187,153,12,62,168,7,170,107,91,43,38,117,173,69,17,61,1,55,82,166,169,30,176,123,48,181,172,159,46,49,104,56,115,116,142,158,180,15,75,4,94,26,74,188,99,88,2,66,182,103,33,72,136,77,19,71,185,45,64,22,76,143,102,108,148,27,111,18,122,146,125,96,131,184,11,164,47,89,157,154,50,186,81,23,151,165,80,113,141,60,174,97,59,162,40,14,139,53,167,73,145,25,138,28,54,100,78,3,126,130,149,152,51,93,189,24,124,31,119,118,120,128,37,175,86,84,150,16,114,129,132,65,58,133,20,57,39,44,183,13,155,85,42,92]
# Discard any samples with less than 80% bases called.
#good_samples = [all_samples[i] for i in worst_to_best[76:]]
#print(good_samples)

# samples that have most loci covered
#good_samples = ['AAGAGGCA-AAGGAGTA', 'PB36', 'PB80', 'PB23', 'PB78', 'TCCTGAGC-AAGGAGTA', 'PB47', 'PB48', 'AAGAGGCA-ACTGCATA', 'PB60', 'PB_1', 'PB91', 'PB95', 'PB19', 'PB59', 'GTAGAGGA-ACTGCATA', 'PB81', 'PB28', 'PB68', 'GTAGAGGA-GTAAGGAG', 'PB73', 'PB93', 'GGACTCCT-GTAAGGAG', 'PB39', 'TAAGGCGA-AGAGTAGA', 'GTAGAGGA-CTAAGCCT', 'PB32', 'PB33', 'PB79', 'PB_9', 'PB64', 'PB66', 'CAGAGAGG-ACTGCATA', 'GTAGAGGA-CTCTCTAT', 'CAGAGAGG-TATCCTCT', 'CTCTCTAC-TATCCTCT', 'AAGAGGCA-CTAAGCCT', 'AGGCAGAA-AGAGTAGA', 'AAGAGGCA-TATCCTCT', 'CAGAGAGG-AAGGAGTA', 'PB53', 'CTCTCTAC-TAGATCGC', 'GCTACGCT-ACTGCATA', 'CAGAGAGG-TAGATCGC', 'GGACTCCT-TAGATCGC', 'CTCTCTAC-ACTGCATA', 'CTCTCTAC-AGAGTAGA', 'PB37', 'PB16', 'GTAGAGGA-TATCCTCT', 'CAGAGAGG-AGAGTAGA', 'AAGAGGCA-GTAAGGAG', 'CAGAGAGG-GTAAGGAG', 'PB77', 'PB29', 'CTCTCTAC-CTCTCTAT', 'PB83', 'TAAGGCGA-AAGGAGTA', 'PB_5', 'PB43', 'TAAGGCGA-ACTGCATA', 'PB40', 'PB45', 'CAGAGAGG-CTAAGCCT', 'AGGCAGAA-CTCTCTAT', 'PB34', 'CTCTCTAC-CTAAGCCT', 'PB85', 'AAGAGGCA-TAGATCGC', 'TCCTGAGC-CTAAGCCT', 'AAGAGGCA-CTCTCTAT', 'PB87', 'CGTACTAG-AGAGTAGA', 'PB25', 'PB10', 'TAAGGCGA-TATCCTCT', 'PB57', 'PB_4', 'TAAGGCGA-CTCTCTAT', 'GTAGAGGA-AGAGTAGA', 'PB88', 'GGACTCCT-ACTGCATA', 'PB18', 'CGTACTAG-AAGGAGTA', 'AGGCAGAA-GTAAGGAG', 'TCCTGAGC-TATCCTCT', 'PB21', 'CGTACTAG-ACTGCATA', 'PB41', 'CAGAGAGG-CTCTCTAT', 'PB76', 'AGGCAGAA-TAGATCGC', 'PB_3', 'TCCTGAGC-AGAGTAGA', 'PB11', 'PB27', 'TAAGGCGA-CTAAGCCT', 'PB44', 'PB84', 'AGGCAGAA-AAGGAGTA', 'PB90', 'CGTACTAG-CTAAGCCT', 'GTAGAGGA-TAGATCGC', 'PB35', 'CGTACTAG-TATCCTCT', 'TCCTGAGC-GTAAGGAG', 'PB82', 'CGTACTAG-TAGATCGC', 'PB50', 'GGACTCCT-CTAAGCCT', 'CGTACTAG-GTAAGGAG', 'AGGCAGAA-TATCCTCT', 'CGTACTAG-CTCTCTAT', 'TCCTGAGC-ACTGCATA', 'PB58', 'PB42']
good_samples = all_samples

positions = []
coverage = []
called = []
alleles = []
meta = []
rec = 0
counter = 0
chromosome = []
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

    # skip if only one proper call or not a proper biallelic
    unique_alleles = set(tmp_alleles)
    if (len(unique_alleles)==2 and 'N' not in unique_alleles) or (len(unique_alleles)==3 and 'N' in unique_alleles):
        alleles.append(tmp_alleles)
        #print(tmp_alleles)
        positions.append(str(pos.pos))     # Positions are specific to the contig so we also need to keep track of the contig number
        chromosome.append(pos.chrom.split('_')[1])
        coverage.append(tmp_coverage)
        called.append([pos.samples[s].alleles[0] is not None for s in good_samples]) # Not a blank space (doesn't consider the coverage requirement so coverage will be higher than what we consider)
        meta.append(tmp_meta)
        rec+=1
        if rec>rec_max:
            break
    else:
        continue

# convert to array (data per strain)
called = np.array(called).T
coverage = np.array(coverage).T
alleles = np.array(alleles).T


# make a masked array of SNPs
snps = np.ma.array(alleles, fill_value='N')
snps.mask = (snps=='-')|(snps=='N')  # We want Ns and gaps to stay labeled 
# get allele counts and major alleles
alpha = np.array(['A', 'C', 'G', 'T'])
af = np.array([np.sum(snps==a, axis=0) for a in 'ACGT'])
maj = alpha[af.argmax(axis=0)]

# SNPs are not major alleles > SNPs = 1, maj = 0
maj_snps = snps!=maj
maj_snps = maj_snps.astype(int)

chr_df = pd.DataFrame(columns=positions)
chr_df.loc['chr'] = chromosome
df = pd.DataFrame(maj_snps, index=good_samples, columns=positions)
df_final = pd.concat([chr_df,df], ignore_index=False)
mdf = pd.DataFrame(meta, index=positions, columns=good_samples)
mdf = mdf.T

#df_final.to_csv('/home/ada/Desktop/Shraiman_lab/srb_data/srb_snp_data_2021_chr.csv')
#mdf.to_csv('../../srb_data/srb_snp_meta_data_2021_concat.csv')



fraction_called = called.mean(axis=1)
fraction_called.sort()
print(fraction_called)
worst_to_best = fraction_called.argsort()

#print(worst_to_best)

