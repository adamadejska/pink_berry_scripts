####################################################################
# Use the background set of GO terms from all genes in the gff file
# to check if experimental set are overrepresented in certain terms. 
#
# How many SNPs are in certain genes of interest? 
# Given a gene code, find how many SNPs and where they are. 
# SRB
####################################################################

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np


gff_file_path = '/home/ada/Desktop/Shraiman_lab/srb_data/otuB_mp_noG.gff'
go_file_path = '/home/ada/Desktop/Shraiman_lab/GO_analysis/SRB_Uniprot_GO_all_genes.tsv'
data = pd.read_csv('/home/ada/Desktop/Shraiman_lab/srb_data/srb_snp_data_2021_chr.csv', index_col=0)

# Index = names of bacterial samples
data_index = data.index.values.tolist()
data_positions = data.columns.tolist()
data_chr = [int(i) for i in data.iloc[0,:].tolist()]


# Example GO file entry
# P03004	P03004	DNAA_ECOLI	Chromosomal replication initiator protein DnaA	dnaA b3702 JW3679	Escherichia coli (strain K12)	DNA replication [GO:0006260]; DNA replication initiation [GO:0006270]; DNA unwinding involved in DNA replication [GO:0006268]; negative regulation of DNA-templated DNA replication initiation [GO:0032297]; positive regulation of DNA-templated DNA replication initiation [GO:0032298]; regulation of DNA replication [GO:0006275]; regulation of DNA-templated DNA replication initiation [GO:0030174]	cytoplasmic side of plasma membrane [GO:0009898]; cytosol [GO:0005829]; DnaA-DiaA complex [GO:1990102]; DnaA-Dps complex [GO:1990084]; DnaA-HU complex [GO:1990103]; DnaA-L2 complex [GO:1990082]; DnaA-oriC complex [GO:1990101]; plasma membrane [GO:0005886]; replication inhibiting complex [GO:1990078]	ATP binding [GO:0005524]; DNA binding [GO:0003677]; DNA replication origin binding [GO:0003688]; identical protein binding [GO:0042802]; sequence-specific DNA binding [GO:0043565]	cytoplasmic side of plasma membrane [GO:0009898]; cytosol [GO:0005829]; DnaA-DiaA complex [GO:1990102]; DnaA-Dps complex [GO:1990084]; DnaA-HU complex [GO:1990103]; DnaA-L2 complex [GO:1990082]; DnaA-oriC complex [GO:1990101]; plasma membrane [GO:0005886]; replication inhibiting complex [GO:1990078]; ATP binding [GO:0005524]; DNA binding [GO:0003677]; DNA replication origin binding [GO:0003688]; identical protein binding [GO:0042802]; sequence-specific DNA binding [GO:0043565]; DNA replication [GO:0006260]; DNA replication initiation [GO:0006270]; DNA unwinding involved in DNA replication [GO:0006268]; negative regulation of DNA-templated DNA replication initiation [GO:0032297]; positive regulation of DNA-templated DNA replication initiation [GO:0032298]; regulation of DNA replication [GO:0006275]; regulation of DNA-templated DNA replication initiation [GO:0030174]

# Example gff entry
# unitig_0_otuB|quiver	Prodigal:2.60	CDS	28	1134	.	-	0	ID=PBSRB_00001;inference=ab initio prediction:Prodigal:2.60;locus_tag=PBSRB_00001;product=hypothetical protein

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
    print(counter)
    # Specific to SRB data
    with open(gff_file_path, 'r') as f:
        for line in f:
            if line.startswith('unitig') and code in line:
                info = line.split('\t')[8]
                info = info.split(';')

                start, end = int(line.split('\t')[3]), int(line.split('\t')[4])
                chromosome = int(line.split('\t')[0].split('_')[1])

    num_of_mutations = 0
    num_of_singletons = 0
    for i in range(0, len(data_positions)):
        pos = data_positions[i]
        whole_pos = int(data_positions[i].split('.')[0])
        chr = int(data_chr[i])

        if whole_pos >= start and whole_pos <= end and chromosome == chr:
            num_of_mutations += 1
            tmp = np.array(data[pos])
            if np.nansum(tmp) == 1:
                num_of_singletons += 1
        elif whole_pos > end:
            break

    if num_of_mutations > 0:
        gene_length.append(end-start)
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
plt.title('Is there a correlation between number of mutations and the length of the gene?')
plt.show()

