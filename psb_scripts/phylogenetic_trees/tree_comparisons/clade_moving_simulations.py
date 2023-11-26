
import numpy as np
import random

###
# MOVE ONE CLADE
# More complex tree (reference) (without one clade) (#,(C,(D,E))) 
# ((AD,(AE,(((((O,(P,R)),S),T),(U,W)),(Y,(Z,((AA,AB),AC)))))),(B,(((F,G),(H,I)),(J,(K,((L,M),N))))));

# Code for the more complex tree simulations
node_set = ['AD','AE','O','P','R','S','T','U','W','Y','Z','AA','AB','AC','B','F','G','H','I','J','K','L','M','N']
n = 500
random_trees = []

structure = '((xxx,(xxx,(((((xxx,(xxx,xxx)),xxx),xxx),(xxx,xxx)),(xxx,(xxx,((xxx,xxx),xxx)))))),(xxx,(((xxx,xxx),(xxx,xxx)),(xxx,(xxx,((xxx,xxx),xxx))))));'
structure = structure.split('xxx')

for i in range(0, n):
    placement = random.randint(0,len(node_set)-1)

    while placement == 15:
        placement = random.randint(0,len(node_set)-1)

    newick_tree = []
    for j in range(0, len(node_set)):
        if j != placement:
            newick_tree.append(structure[j]) 
            newick_tree.append(node_set[j])
        else:
            newick_tree.append(structure[j])
            newick_tree.append('(' + node_set[j] + ',(C,(D,E)))')
    
    
    newick_tree.append(structure[-1])
    newick_tree_str = "".join(newick_tree)
    random_trees.append(newick_tree_str)

outfile = "/home/ada/Desktop/PinkBerry_scripts_paper/psb_scripts/phylogenetic_trees/tree_comparisons/simulation_random_trees_one_clade.txt"
out = open(outfile, 'w')
for i in random_trees:
    out.write(i + '\n')


###
# MOVE TWO CLADES
# More complex tree (reference) (without one clade) (#,(C,(D,E))) (#,(AA,(AB,AC)))
# ((AD,(AE,(((((O,(P,R)),S),T),(U,W)),(Y,Z)))),(B,(((F,G),(H,I)),(J,(K,((L,M),N))))));

# Code for the more complex tree simulations
node_set = ['AD','AE','O','P','R','S','T','U','W','Y','Z','B','F','G','H','I','J','K','L','M','N']
n = 500
random_trees = []

structure = '((xxx,(xxx,(((((xxx,(xxx,xxx)),xxx),xxx),(xxx,xxx)),(xxx,xxx)))),(xxx,(((xxx,xxx),(xxx,xxx)),(xxx,(xxx,((xxx,xxx),xxx))))));'
structure = structure.split('xxx')

for i in range(0, n):
    placement1, placement2= random.sample(range(0, len(node_set)-1), 2)

    while placement1 == 12 or placement2 == 10:
        placement1, placement2= random.sample(range(0, len(node_set)-1), 2)

    newick_tree = []
    for j in range(0, len(node_set)):
        if j != placement1 and j != placement2:
            newick_tree.append(structure[j]) 
            newick_tree.append(node_set[j])
        elif j == placement1:
            newick_tree.append(structure[j])
            newick_tree.append('(' + node_set[j] + ',(C,(D,E)))')
        elif j == placement2:
            newick_tree.append(structure[j])
            newick_tree.append('(' + node_set[j] + ',(AA,(AB,AC)))')
    
    
    newick_tree.append(structure[-1])
    newick_tree_str = "".join(newick_tree)
    random_trees.append(newick_tree_str)

outfile = "/home/ada/Desktop/PinkBerry_scripts_paper/psb_scripts/phylogenetic_trees/tree_comparisons/simulation_random_trees_two_clades.txt"
out = open(outfile, 'w')
for i in random_trees:
    out.write(i + '\n')

###
# MOVE THREE CLADES
# More complex tree (reference) (without one clade) (#,(C,(D,E))) (#,(AA,(AB,AC))) (#,(O,(P,R)))
# ((AD,(AE,(((S,T),(U,W)),(Y,Z)))),(B,(((F,G),(H,I)),(J,(K,((L,M),N))))));

# Code for the more complex tree simulations
node_set = ['AD','AE','S','T','U','W','Y','Z','B','F','G','H','I','J','K','L','M','N']
n = 500
random_trees = []

structure = '((xxx,(xxx,(((xxx,xxx),(xxx,xxx)),(xxx,xxx)))),(xxx,(((xxx,xxx),(xxx,xxx)),(xxx,(xxx,((xxx,xxx),xxx))))));'
structure = structure.split('xxx')

#print(structure)
for i in range(0, n):
    placement1, placement2, placement3 = random.sample(range(0, len(node_set)-1), 3)
    
    while placement1 == 9 or placement2 == 7 or placement3 == 2:
        placement1, placement2, placement3 = random.sample(range(0, len(node_set)-1), 3)

    newick_tree = []
    for j in range(0, len(node_set)):
        if j != placement1 and j != placement2 and j != placement3:
            newick_tree.append(structure[j]) 
            newick_tree.append(node_set[j])
        elif j == placement1:
            newick_tree.append(structure[j])
            newick_tree.append('(' + node_set[j] + ',(C,(D,E)))')
        elif j == placement2:
            newick_tree.append(structure[j])
            newick_tree.append('(' + node_set[j] + ',(AA,(AB,AC)))')
        elif j == placement3:
            newick_tree.append(structure[j])
            newick_tree.append('(' + node_set[j] + ',(O,(P,R)))')
    
    
    newick_tree.append(structure[-1])
    newick_tree_str = "".join(newick_tree)
    random_trees.append(newick_tree_str)

outfile = "/home/ada/Desktop/PinkBerry_scripts_paper/psb_scripts/phylogenetic_trees/tree_comparisons/simulation_random_trees_three_clades.txt"
out = open(outfile, 'w')
for i in random_trees:
    out.write(i + '\n')

###
# MOVE FOUR CLADES
# More complex tree (reference) (without one clade) (#,(C,(D,E))) (#,(AA,(AB,AC))) (#,(O,(P,R))) (#, (N,(L,M)))
# ((AD,(AE,(((S,T),(U,W)),(Y,Z)))),(B,(((F,G),(H,I)),(J,K))));

# Code for the more complex tree simulations
node_set = ['AD','AE','S','T','U','W','Y','Z','B','F','G','H','I','J','K']
n = 500
random_trees = []

structure = '((xxx,(xxx,(((xxx,xxx),(xxx,xxx)),(xxx,xxx)))),(xxx,(((xxx,xxx),(xxx,xxx)),(xxx,xxx))));'
structure = structure.split('xxx')

#print(structure)
for i in range(0, n):
    placement1, placement2, placement3, placement4 = random.sample(range(0, len(node_set)-1), 4)
    
    while placement1 == 9 or placement2 == 7 or placement3 == 2 or placement4 == 14:
        placement1, placement2, placement3, placement4 = random.sample(range(0, len(node_set)-1), 4)

    newick_tree = []
    for j in range(0, len(node_set)):
        if j != placement1 and j != placement2 and j != placement3 and j != placement4:
            newick_tree.append(structure[j]) 
            newick_tree.append(node_set[j])
        elif j == placement1:
            newick_tree.append(structure[j])
            newick_tree.append('(' + node_set[j] + ',(C,(D,E)))')
        elif j == placement2:
            newick_tree.append(structure[j])
            newick_tree.append('(' + node_set[j] + ',(AA,(AB,AC)))')
        elif j == placement3:
            newick_tree.append(structure[j])
            newick_tree.append('(' + node_set[j] + ',(O,(P,R)))')
        elif j == placement4:
            newick_tree.append(structure[j])
            newick_tree.append('(' + node_set[j] + ',(N,(L,M)))')
    
    
    newick_tree.append(structure[-1])
    newick_tree_str = "".join(newick_tree)
    random_trees.append(newick_tree_str)

outfile = "/home/ada/Desktop/PinkBerry_scripts_paper/psb_scripts/phylogenetic_trees/tree_comparisons/simulation_random_trees_four_clades.txt"
out = open(outfile, 'w')
for i in random_trees:
    out.write(i + '\n')