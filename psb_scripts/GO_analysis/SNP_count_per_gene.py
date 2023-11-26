####################################################################
# Use the background set of GO terms from all genes in the gff file
# to check if experimental set are overrepresented in certain terms. 
#
# How many SNPs are in certain genes of interest? 
# Given a gene code, find how many SNPs and where they are. 
# PSB
####################################################################

from collections import Counter
import re
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np


gff_file_path = '/home/ada/Desktop/Shraiman_lab/data/psb_scaff03_noG.gff'
go_file_path = '/home/ada/Desktop/Shraiman_lab/GO_analysis/PSB_Uniprot_GO_all_genes.tsv'
data = pd.read_csv('/home/ada/Desktop/Shraiman_lab/data/snp_data_2020.csv', index_col=0)

# Index = names of bacterial samples
data_index = data.index.values.tolist()
data_positions = data.columns.tolist()

# Example GO file entry
# P03004	P03004	DNAA_ECOLI	Chromosomal replication initiator protein DnaA	dnaA b3702 JW3679	Escherichia coli (strain K12)	DNA replication [GO:0006260]; DNA replication initiation [GO:0006270]; DNA unwinding involved in DNA replication [GO:0006268]; negative regulation of DNA-templated DNA replication initiation [GO:0032297]; positive regulation of DNA-templated DNA replication initiation [GO:0032298]; regulation of DNA replication [GO:0006275]; regulation of DNA-templated DNA replication initiation [GO:0030174]	cytoplasmic side of plasma membrane [GO:0009898]; cytosol [GO:0005829]; DnaA-DiaA complex [GO:1990102]; DnaA-Dps complex [GO:1990084]; DnaA-HU complex [GO:1990103]; DnaA-L2 complex [GO:1990082]; DnaA-oriC complex [GO:1990101]; plasma membrane [GO:0005886]; replication inhibiting complex [GO:1990078]	ATP binding [GO:0005524]; DNA binding [GO:0003677]; DNA replication origin binding [GO:0003688]; identical protein binding [GO:0042802]; sequence-specific DNA binding [GO:0043565]	cytoplasmic side of plasma membrane [GO:0009898]; cytosol [GO:0005829]; DnaA-DiaA complex [GO:1990102]; DnaA-Dps complex [GO:1990084]; DnaA-HU complex [GO:1990103]; DnaA-L2 complex [GO:1990082]; DnaA-oriC complex [GO:1990101]; plasma membrane [GO:0005886]; replication inhibiting complex [GO:1990078]; ATP binding [GO:0005524]; DNA binding [GO:0003677]; DNA replication origin binding [GO:0003688]; identical protein binding [GO:0042802]; sequence-specific DNA binding [GO:0043565]; DNA replication [GO:0006260]; DNA replication initiation [GO:0006270]; DNA unwinding involved in DNA replication [GO:0006268]; negative regulation of DNA-templated DNA replication initiation [GO:0032297]; positive regulation of DNA-templated DNA replication initiation [GO:0032298]; regulation of DNA replication [GO:0006275]; regulation of DNA-templated DNA replication initiation [GO:0030174]

# Example gff entry
# psb-scaff03	Prodigal:2.6	CDS	1660	2769	.	+	0	ID=P1s_00002;Name=dnaN;dbxref=COG:COG0592;gene=dnaN;inference=ab initio prediction:Prodigal:2.6,similar to AA sequence:UniProtKB:Q9I7C4;locus_tag=P1s_00002;product=Beta sliding clamp

# Read in the GO file and get all the gene codes
gene_codes = []
with open(go_file_path, 'r') as f:
    f.readline() # skip header
    for line in f:
        line = line.split('\t')
        gene_codes.append(line[1])


# Find all mutations in each gene
gene_length, gene_SNPs = [], []
counter = 0
for code in gene_codes:
    counter += 1
    #print(counter)
    # Specific to PSB data
    with open(gff_file_path, 'r') as f:
        for line in f:
            if line.startswith('psb-scaff03') and code in line:
                info = line.split('\t')[8]
                info = info.split(';')

                start, end = int(line.split('\t')[3]), int(line.split('\t')[4])

                gene_length.append(end-start)

                num_of_mutations = 0
                for i in data_positions:
                    if int(i) >= start and int(i) <= end:
                        num_of_mutations += 1
                    elif int(i) > end:
                        break

                gene_SNPs.append(num_of_mutations)

# convert lists to np arrays
gene_length = np.array(gene_length)
gene_SNPs = np.array(gene_SNPs)

# Find a best fit line
a, b = np.polyfit(gene_length, gene_SNPs, 1)
std_SNPs = np.std(gene_SNPs)

# Plot the length of the gene vs the number of mutations it has
plt.scatter(gene_length, gene_SNPs, alpha=0.5)
#add line of best fit to plot
plt.plot(gene_length, a*gene_length+b, c='r')
plt.plot(gene_length, a*gene_length+b+std_SNPs, c='y')
plt.plot(gene_length, a*gene_length+b-std_SNPs, c='y')
plt.xlabel('Gene Length (nucleotides)')
plt.ylabel('Number of positions with mutations')
plt.ylim(0,100)
plt.xlim(0, 8000)
plt.title('Is there a correlation between number of mutations and the length of the gene?')
plt.show()

