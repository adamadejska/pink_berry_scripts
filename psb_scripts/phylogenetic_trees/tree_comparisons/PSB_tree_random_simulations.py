####################################################################################
# Create a set of newick trees with a random set of nodes based on the original structure
# of the PSB tree.
# This will be used as an inpuit for Visual TreeCmp to calculate distances between
# the reference tree and the random trees and trees without horizontall transferred regions.
# The discance values will be then plotted as a distribution to understand how different
# the experimental tree really is from the original tree.
####################################################################################

import numpy as np
import random

# Original PSB tree
# (outgroup,(AGGCAGAA-AGAGTAGA,(((PB26,PB34),(PB11,(CTCTCTAC-GTAAGGAG,(CAGAGAGG-TATCCTCT,(CTCTCTAC-AAGGAGTA,(CTCTCTAC-AGAGTAGA,(CAGAGAGG-AAGGAGTA,(CTCTCTAC-TATCCTCT,(CAGAGAGG-ACTGCATA,(CAGAGAGG-AGAGTAGA,((CTCTCTAC-TAGATCGC,(CTCTCTAC-ACTGCATA,(GTAGAGGA-CTAAGCCT,(CTCTCTAC-CTCTCTAT,(CAGAGAGG-CTCTCTAT,CTCTCTAC-CTAAGCCT))))),(CAGAGAGG-GTAAGGAG,(CAGAGAGG-CTAAGCCT,CAGAGAGG-TAGATCGC))))))))))))),((PB44,(PB52,(PB59,PB67))),(TAAGGCGA-TATCCTCT,(((PB_5,(PB37,(PB13,PB85))),(PB69,((PB45,PB77),((GGACTCCT-ACTGCATA,GGACTCCT-TAGATCGC),(PB61,PB78))))),((PB28,(AGGCAGAA-ACTGCATA,AGGCAGAA-CTAAGCCT)),(PB64,(PB48,((PB47,(PB16,(PB32,(AAGAGGCA-CTAAGCCT,AAGAGGCA-CTCTCTAT)))),((AAGAGGCA-GTAAGGAG,GCTACGCT-ACTGCATA),(PB31,((AAGAGGCA-AAGGAGTA,(GCTACGCT-GTAAGGAG,(AAGAGGCA-ACTGCATA,AAGAGGCA-AGAGTAGA))),(PB55,(PB24,(PB_8,(PB39,(AAGAGGCA-TAGATCGC,(PB63,(PB80,((AAGAGGCA-TATCCTCT,(PB40,PB87)),((PB73,(TCCTGAGC-CTCTCTAT,TCCTGAGC-TAGATCGC)),(GGACTCCT-CTAAGCCT,(GGACTCCT-AGAGTAGA,GGACTCCT-GTAAGGAG))))))))))))))))))))))));

# Code for the more complex tree simulations
node_set = ['outgroup','PB26','PB34','PB11','CTCTCTAC-GTAAGGAG','CAGAGAGG-TATCCTCT','CAGAGAGG-AAGGAGTA','CAGAGAGG-ACTGCATA','CAGAGAGG-AGAGTAGA','CTCTCTAC-TATCCTCT','CTCTCTAC-TAGATCGC','CTCTCTAC-ACTGCATA','GTAGAGGA-CTAAGCCT','CTCTCTAC-CTCTCTAT','CAGAGAGG-CTCTCTAT','CTCTCTAC-CTAAGCCT','CAGAGAGG-GTAAGGAG','CAGAGAGG-CTAAGCCT','CAGAGAGG-TAGATCGC','CTCTCTAC-AAGGAGTA','CTCTCTAC-AGAGTAGA','AGGCAGAA-AGAGTAGA','PB44','PB52','PB59','PB67','PB73','TCCTGAGC-CTCTCTAT','TCCTGAGC-TAGATCGC','GGACTCCT-CTAAGCCT','GGACTCCT-AGAGTAGA','GGACTCCT-GTAAGGAG','TAAGGCGA-TATCCTCT','PB13','PB85','PB37','PB_5','PB69','PB45','PB77','GGACTCCT-ACTGCATA','GGACTCCT-TAGATCGC','PB61','PB78','PB64','PB48','PB_8','PB39','AAGAGGCA-TAGATCGC','PB80','PB40','PB87','AAGAGGCA-AAGGAGTA','PB24','GCTACGCT-ACTGCATA','AAGAGGCA-GTAAGGAG','PB47','PB16','PB32','AAGAGGCA-CTAAGCCT','AAGAGGCA-CTCTCTAT']
n = 1000
random_trees = []
#print(len(node_set))

#structure = "(xxx,(xxx,(((xxx,xxx),(xxx,(xxx,(xxx,(xxx,(xxx,(xxx,(xxx,(xxx,(xxx,((xxx,(xxx,(xxx,(xxx,(xxx,xxx))))),(xxx,(xxx,xxx))))))))))))),((xxx,(xxx,(xxx,xxx))),(xxx,(((xxx,(xxx,(xxx,xxx))),(xxx,((xxx,xxx),((xxx,xxx),(xxx,xxx))))),((xxx,(xxx,xxx)),(xxx,(xxx,((xxx,(xxx,(xxx,(xxx,xxx)))),((xxx,xxx),(xxx,((xxx,(xxx,(xxx,xxx))),(xxx,(xxx,(xxx,(xxx,(xxx,(xxx,(xxx,((xxx,(xxx,xxx)),((xxx,(xxx,xxx)),(xxx,(xxx,xxx))))))))))))))))))))))));"

structure = "(xxx,(((xxx,xxx),(xxx,(xxx,(xxx,(xxx,((xxx,(xxx,(xxx,((xxx,(xxx,(xxx,(xxx,(xxx,xxx))))),(xxx,(xxx,xxx)))))),(xxx,xxx))))))),((xxx,(xxx,(xxx,(xxx,(xxx,((xxx,(xxx,xxx)),(xxx,(xxx,xxx)))))))),((xxx,(((xxx,xxx),(xxx,xxx)),(xxx,((xxx,xxx),((xxx,xxx),(xxx,xxx)))))),(xxx,(xxx,(((xxx,(xxx,(xxx,(xxx,(xxx,xxx))))),(xxx,xxx)),(xxx,(xxx,(xxx,(xxx,(xxx,(xxx,xxx)))))))))))));"
structure = structure.split('xxx')

for i in range(0, n):
    random.shuffle(node_set)

    newick_tree = []
    for j in range(0, len(node_set)):
        newick_tree.append(structure[j]) 
        newick_tree.append(node_set[j])
    
    
    newick_tree.append(structure[-1])
    newick_tree_str = "".join(newick_tree)
    random_trees.append(newick_tree_str)

outfile = "/home/ada/Desktop/PinkBerry_scripts_paper/psb_scripts/phylogenetic_trees/tree_comparisons/PSB_random_trees.txt"
out = open(outfile, 'w')
for i in random_trees:
    out.write(i + '\n')


# No 8 clade SNPs tree
# (outgroup,(((PB26,PB34),(PB11,(CTCTCTAC-GTAAGGAG,(CAGAGAGG-TATCCTCT,(CAGAGAGG-AAGGAGTA,((CAGAGAGG-ACTGCATA,(CAGAGAGG-AGAGTAGA,(CTCTCTAC-TATCCTCT,((CTCTCTAC-TAGATCGC,(CTCTCTAC-ACTGCATA,(GTAGAGGA-CTAAGCCT,(CTCTCTAC-CTCTCTAT,(CAGAGAGG-CTCTCTAT,CTCTCTAC-CTAAGCCT))))),(CAGAGAGG-GTAAGGAG,(CAGAGAGG-CTAAGCCT,CAGAGAGG-TAGATCGC)))))),(CTCTCTAC-AAGGAGTA,CTCTCTAC-AGAGTAGA))))))),((AGGCAGAA-AGAGTAGA,(PB44,(PB52,(PB59,(PB67,((PB73,(TCCTGAGC-CTCTCTAT,TCCTGAGC-TAGATCGC)),(GGACTCCT-CTAAGCCT,(GGACTCCT-AGAGTAGA,GGACTCCT-GTAAGGAG)))))))),((TAAGGCGA-TATCCTCT,(((PB13,PB85),(PB37,PB_5)),(PB69,((PB45,PB77),((GGACTCCT-ACTGCATA,GGACTCCT-TAGATCGC),(PB61,PB78)))))),(PB64,(PB48,(((PB_8,(PB39,(AAGAGGCA-TAGATCGC,(PB80,(PB40,PB87))))),(AAGAGGCA-AAGGAGTA,PB24)),(GCTACGCT-ACTGCATA,(AAGAGGCA-GTAAGGAG,(PB47,(PB16,(PB32,(AAGAGGCA-CTAAGCCT,AAGAGGCA-CTCTCTAT)))))))))))));
