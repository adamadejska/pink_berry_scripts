################################################################################
# Since SRB genome is divided into contigs, make fasta sequences for each contig.
# Go through the bcf file and create 'chromosomes' based on the contig ID and the reference allele
################################################################################

from pysam import VariantFile
import numpy as np
import sys
import pandas as pd
import matplotlib.pyplot as plt


np.set_printoptions(threshold=sys.maxsize)

bcf_in = VariantFile("/home/ada/Desktop/Shraiman_lab/srb_data/srb.minimap2-freebayes.bcf")  # auto-detect input format
bcf_iter = bcf_in.fetch()

contig_to_seq = {}

# Go through each position in the bcf file. For each contig, create a sequence from the reference allele.
for pos in bcf_iter:
    if pos.chrom not in contig_to_seq.keys():
        contig_to_seq[pos.chrom] = pos.ref
    else:
        contig_to_seq[pos.chrom] += pos.ref


# Save the created chromosome sequences into a fasta-like file
out = open('/home/ada/Desktop/Shraiman_lab/srb_data/srb_contigs.fasta', 'w')
for chrom in contig_to_seq.keys():
    out.write('>' + chrom + '\n')
    out.write(contig_to_seq[chrom] + '\n')
    
out.close()