####################################################################
# Use the background set of GO terms from all genes in the gff file
# to check if experimental set are overrepresented in certain terms. 
#
# How many SNPs are in certain genes of interest? 
# Given a gene code, find how many SNPs and where they are. 
# PSB
####################################################################

from collections import Counter
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np


# below  are low counts codes for PSB 10 SNPs or more for 1 kb regions
gene_code = ['P0A8A4', 'A0P8X0', 'D0ZRB2', 'Q83ES6', 'P33941', 'P37877', 'Q81ST8', 'Q8ZNY9', 'D4G3R4', 'O34739', 'P46481', 'Q0VZ70', 'P76052', 'P76149', 'Q0P9C7', 'P23837', 'P45745', 'P69348', 'Q9KJ21', 'P9WN57', 'P76562', 'P75990', 'P0AES0', 'P9WMJ1', 'O87941', 'Q9KJ22', 'P0A734', 'P07061', 'P17054', 'P26647']
gene_codes = ['P0AEZ3','P33919','C0SP86','P60240', 'Q7A801', 'A0A0H2VDN9', 'D0ZRB2', 'E0SGL7', 'A1IGV8', 'Q9KIG4', 'G4NYJ6']
gene_names = ['Septum site-determining protein MinD', 'Putative DNA repair helicase RadD','DNA translocase SftA', 'RNA polymerase-associated protein RapA', 'Type-1 restriction enzyme R protein hsdR',  'Secretory immunoglobulin A-binding protein EsiB', 'E3 ubiquitin-protein ligase SlrP', 'Putative deoxyribonuclease RhsC', 'Alpha-agarase agaA33', 'Serine/threonine-protein kinase PK-1 spk1', 'tRNA3(Ser)-specific nuclease WapA']
gene_codes = ['P1s_03747']
gene_names = ['Alpha-agarase']
gff_file_path = '/home/ada/Desktop/Shraiman_lab/data/psb_scaff03_noG.gff'

# Example gff entry
# psb-scaff03	Prodigal:2.6	CDS	1660	2769	.	+	0	ID=P1s_00002;Name=dnaN;dbxref=COG:COG0592;gene=dnaN;inference=ab initio prediction:Prodigal:2.6,similar to AA sequence:UniProtKB:Q9I7C4;locus_tag=P1s_00002;product=Beta sliding clamp

# Specific to PSB data
for n in range(0,len(gene_codes)):
    gene_code = gene_codes[n]
    print(gene_code)
    with open(gff_file_path, 'r') as f:
        for line in f:
            if line.startswith('psb-scaff03'):
                    if gene_code in line:
                        info = line.split('\t')[8]
                        info = info.split(';')

                        start, end = int(line.split('\t')[3]), int(line.split('\t')[4])
                        print(line)

    data = pd.read_csv('/home/ada/Desktop/Shraiman_lab/data/snp_data_2020.csv', index_col=0)

    # Index = names of bacterial samples
    data_index = data.index.values.tolist()
    data_positions = data.columns.tolist()

    # Each mutation is an array where each entry represents a specific berry (index).
    # Make a dictionary for bacteria and count the number of mutations each have per gene.
    # Inititate the dictionary and set every berry to 0. 
    berry_specific_mutations = {}
    for berry in data_index:
         berry_specific_mutations[berry] = 0

    for i in data_positions:
        if int(i) >= start and int(i) <= end:

            tmp = np.array(data[i])
            for j in range(0, len(data_index)):
                if not np.isnan(tmp[j]):
                    berry_specific_mutations[data_index[j]] += tmp[j]
    
    # Take only the berries that have mutations. We don't want to plot a bunch of zeros.
    valid_berries, num_of_mutations = [], []
    for k,v in berry_specific_mutations.items():
        if v > 0:
            valid_berries.append(k)
            num_of_mutations.append(v)

    num_of_mutations, valid_berries = zip(*sorted(zip(num_of_mutations, valid_berries)))

    # color code the berries based on their geographical locations
    location_file = '/home/ada/Desktop/Shraiman_lab/data/sample_sites.txt'
    loc_to_color = {'A': 'red', 'F': 'orange', 'C': 'green', 'B': 'blue', 'E': 'magenta'}
    colors = []
    for i in valid_berries:
        with open(location_file, 'r') as f:
            for line in f:
                strain, location = line.split()

                if i == strain:
                    colors.append(loc_to_color[location.strip()])


    plt.barh(valid_berries, num_of_mutations, color=colors)
    plt.title(gene_names[n])
    plt.xlabel('Number of mutations per berry')
    plt.ylabel('Berry names')
    plt.grid(axis = 'x')
    plt.show()

