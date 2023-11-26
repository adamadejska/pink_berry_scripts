    ####################################################################
# Use the background set of GO terms from all genes in the gff file
# to check if experimental set are overrepresented in certain terms. 
# PSB
####################################################################

from collections import Counter
import re
import matplotlib.pyplot as plt
import pandas as pd
import random

go_file_path = '/home/ada/Desktop/Shraiman_lab/GO_analysis/SRB_Uniprot_GO_all_genes.tsv'

# Open the file with reference GO terms (all genes from gff file)
# Divide the GO terms into 3 lists: biological processes, cellular component, molecular function

bio_processes_all = []
cell_components_all = []
mol_functions_all = []

id_to_GO_codes = {}

with open(go_file_path, 'r') as f:
    f.readline()  # skip header: From	Entry	Entry Name	Protein names	Gene Names	Gene Ontology (biological process)	Gene Ontology (cellular component)	Gene Ontology (molecular function)	Gene Ontology (GO)	Organism

    for line in f:
        line = line.split('\t')
        id = line[0].strip()
        bp = line[5].split(';')
        cc = line[6].split(';')
        mf = line[7].split(';')

        # Extract the GO codes only for each category
        if len(bp[0]) > 0:
            bp = [re.findall('[GO:[0-9]+\]', i) for i in bp]  # use RegEx to find all GO codes
            bp = [i[0].split(':')[-1][:-1] for i in bp]
            bio_processes_all += bp

        if len(cc[0]) > 0:
            cc = [re.findall('[GO:[0-9]+\]', i) for i in cc]  # use RegEx to find all GO codes
            cc = [i[0].split(':')[-1][:-1] for i in cc]
            cell_components_all += cc

        if len(mf[0]) > 0:
            mf = [re.findall('[GO:[0-9]+\]', i) for i in mf]  # use RegEx to find all GO codes
            mf = [i[0].split(':')[-1][:-1] for i in mf]
            mol_functions_all += mf

        # Also save all info in a dictionary
        id_to_GO_codes[id] = [bp, cc, mf]

bio_processes_all = Counter(bio_processes_all)
cell_components_all = Counter(cell_components_all)
mol_functions_all = Counter(mol_functions_all)


def plot_background_GO_codes(process_counter_dict, process_name):

    # Plot background GO codes as bar graphs (horizontal)
    # Take only the most frequent GO codes (ignore the singular ones)
    code, freq = [], []
    for k,v in process_counter_dict.items():
        if v >= 10:
            code.append(k)
            freq.append(v)

    # Sort the values
    code = [x for _,x in sorted(zip(freq,code))]
    freq.sort()

    # Plot
    plt.barh(code, freq)
    plt.title(process_name + ' GO terms frequency in all SRB genes')
    plt.ylabel('GO term code')
    plt.xlabel('number of occurrences (out of 1539)')
    plt.show()

    return()

#plot_background_GO_codes(bio_processes_all, 'Biological processes')
#plot_background_GO_codes(cell_components_all, 'Cellular components')
#plot_background_GO_codes(mol_functions_all, 'Molecular functions')


# Get the list of genes present in the high SNP density regions.
# Compare the GO codes to the background for any that might 
# be overrepresented.
psb_snp_dense_genes_path = '/home/ada/Desktop/Shraiman_lab/GO_analysis/SRB_dense_regions_1kb_30_SNPs_ormore.txt'

dense_ids = []
# Read in the file and extract the Uniprot IDs
with open(psb_snp_dense_genes_path, 'r') as f:
    for line in f:
        if line.split(' ')[2].strip() != 'NA':
            dense_ids.append(line.split(' ')[2].strip())

# Get rid of duplicate gene IDs
dense_ids = list(set(dense_ids))
n = len(dense_ids)
dense_bp, dense_cc, dense_mf = [], [], []

# Extract the GO info for each ID
for id in dense_ids:
    if len(id_to_GO_codes[id][0][0]) > 0:
        dense_bp += id_to_GO_codes[id][0]

    if len(id_to_GO_codes[id][1][0]) > 0:
        dense_cc += id_to_GO_codes[id][1]

    if len(id_to_GO_codes[id][2][0]) > 0:
        dense_mf += id_to_GO_codes[id][2]

dense_bp, dense_cc, dense_mf = Counter(dense_bp), Counter(dense_cc), Counter(dense_mf)

def compare_experimental_to_background(background_counter, random_counter, experimental_counter, title):
    """
    Given a counter with the background frequency of the GO terms, compare it to
    the experimental. Plot the experimental counters against the background to see if they are significant.
    """

    exp_code, exp_freq = [], []
    for k,v in experimental_counter.items():
        if v > 1:   # Get rid of single hits for graph readability
            exp_code.append(k)
            exp_freq.append(v)

    bg_code, bg_freq = [], []
    for code in exp_code:    # Use only the codes that were found in experimental, we don't care about others
        bg_freq.append(background_counter[code])
        bg_code.append(code)

    rnd_code, rnd_freq = [], []
    for code in exp_code:    # Use only the codes that were found in experimental, we don't care about others
        rnd_freq.append(random_counter[code])
        rnd_code.append(code)

    # Sort by background values
    bg_freq, rnd_freq, exp_freq, exp_code = zip(*sorted(zip(bg_freq, rnd_freq, exp_freq, exp_code)))

    df = pd.DataFrame({'experimental': exp_freq,'background': bg_freq, 'random': rnd_freq}, index=exp_code)
    ax = df.plot.barh( width=0.9)
    plt.title(title)
    plt.ylabel('GO codes')
    plt.xlabel('Number of occurrences (log)')
    plt.xscale('log')
    plt.show()

# Simulate what we would expect to see from a purely randomized set 

# Make a list that includes all codes with the number of copies that reflects the frequency of that code in the full set
# Do it a few times and take an average to get a better idea of the probabilities.
m = 100.0
rand_bio_processes, rand_cell_components, rand_mol_functions = [], [], []

for i in range(0, int(m)):
    rand_genes = random.sample(id_to_GO_codes.keys(), n)

    for g in rand_genes:
        rand_bio_processes += id_to_GO_codes[g][0]
        rand_cell_components += id_to_GO_codes[g][1]
        rand_mol_functions += id_to_GO_codes[g][2]

rand_bio_processes = Counter(rand_bio_processes)
rand_cell_components = Counter(rand_cell_components)
rand_mol_functions = Counter(rand_mol_functions)

for k, v in rand_bio_processes.items():
    rand_bio_processes[k] = v/m
    
for k, v in rand_cell_components.items():
    rand_cell_components[k] = v/m

for k, v in rand_mol_functions.items():
    rand_mol_functions[k] = v/m

# Make graphs for PSB genes in SNP dense regions
compare_experimental_to_background(bio_processes_all, rand_bio_processes, dense_bp, 'Biological processes experimental vs background GO code counts\n SRB genes from SNP dense regions (20+ SNPs per 1 kb)')
compare_experimental_to_background(cell_components_all, rand_cell_components, dense_cc, 'Cellular Components experimental vs background GO code counts\n SRB genes from SNP dense regions (20+ SNPs per 1 kb)')
compare_experimental_to_background(mol_functions_all, rand_mol_functions, dense_mf, 'Molecular functions experimental vs background GO code counts\n SRB genes from SNP dense regions (20+ SNPs per 1 kb)')

"""
# Which genes have some of the low n codes? 
low_n_codes = []
exp_code= []
for k,v in dense_mf.items():
    if v > 1:   # Get rid of single hits for graph readability
        exp_code.append(k)

for code in exp_code:    # Use only the codes that were found in experimental, we don't care about others
    if mol_functions_all[code] <= 2:
        low_n_codes.append(code)

#print(low_n_codes)
# Find those codes in the genes from experimental set.
low_n_codes_ids = []
for id in dense_ids:
    for code in low_n_codes:
        if code in id_to_GO_codes[id][2][0]:   # Check only biological processes
            low_n_codes_ids.append(id)

#print(low_n_codes_ids)
"""