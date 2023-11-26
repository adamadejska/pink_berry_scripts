###############################################################
# Given a list of loci separated by a comma, give a list of genes
# (their UniProt ID if applicable) that the loci are in. 
# PSB ANALYSIS
###############################################################


loci_file_path = '/home/ada/Desktop/Shraiman_lab/GO_analysis/PSB_loci_in_blocks.csv'
gff_file_path = '/home/ada/Desktop/Shraiman_lab/data/psb_scaff03_noG.gff'

out_path = '/home/ada/Desktop/Shraiman_lab/GO_analysis/PSB_genes_in_blocks.csv'
out = open(out_path, 'w')

# Extract all loci and save them as a list of ints
loci = []
with open(loci_file_path, 'r') as f:
    for line in f:
        loci.append(int(line.strip()))

# Specific to PSB data
with open(gff_file_path, 'r') as f:
    out.write('gene_name, ID, dbxref, UniProt_ID\n')
    for line in f:
        if line.startswith('psb-scaff') and 'hypothetical protein' not in line:
            
            start, end = int(line.split('\t')[3]), int(line.split('\t')[4])

            not_found = True
            for i in loci:
                if i <= end and i >= start and not_found:
                    
                    not_found = False

                    info = line.split('\t')[8]
                    info = info.split(';')

                    try:
                        id_index = [i for i, elem in enumerate(info) if 'ID=' in elem]
                        id = info[id_index[0]].split('=')[-1]   # GFF specific ID
                    except IndexError:
                        id = 'NA'

                    try:
                        gene_name_index = [i for i, elem in enumerate(info) if 'gene=' in elem]
                        gene_name = info[gene_name_index[0]].split('=')[-1]  # gene name (short)
                    except IndexError:
                        gene_name = 'NA'

                    try:
                        dbxref_index = [i for i, elem in enumerate(info) if 'dbxref=' in elem]
                        dbxref = info[dbxref_index[0]].split('=')[-1]    # database specific OD (from NIH)
                    except IndexError:
                        dbxref = 'NA'

                    try:
                        uniprot_index = [i for i, elem in enumerate(info) if 'UniProtKB:' in elem]
                        uniprot_id = info[uniprot_index[0]].split(':')[-1]  # UniProt ID (probably the most important)
                    except IndexError:
                        uniprot_id = 'NA'

                    out.write(gene_name + ',' + uniprot_id + '\n')