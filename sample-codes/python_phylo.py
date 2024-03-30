# Print a phylogenetic tree of complex tandem repeats from the distance matrix of the repeats output by Riki Kawahara's program

# Parse the arguments
# eg. python3 python_phylo.py (input distance matrix file) (input fasta file) (output pdf file)

import sys  # argv
args = sys.argv
mFile = args[1] # distance matrix file
aFile = args[2] # fasta file
oFile = args[3]

dumpFlag = 0
if(len(args)==5):
    dumpFlag = 1
    dFile = args[4]
    dFP = open(dFile, 'w')

# Calculate the lower triangular matrix for Bio.Phylo.TreeConstruction
import numpy as np
mFP = open(mFile, 'r')
aFP = open(aFile, 'r')
data2D = []
name = []
j=0
while 1:
    # Process the distance matrix file mFP
    tmp = mFP.readline().split(' ')
    if not tmp: break   # End of file
    if bool(tmp == ['']): break # End of file
    data1D = []
    if(0 < len(tmp)):
        for i in range(0,min(j+1,len(tmp))): #lower triangular matrix
            data1D.append( float(tmp[i]) )
        data2D.append(data1D)
    
    # Process an annotation in the fasta file aFP
    # Get annotation
    tmp = aFP.readline().split(' ')
    # Put the unique ID before the annotation
    oneName = str(j)+':'+tmp[1]
    name.append( oneName )
    # Ignore string
    aFP.readline();
    j=j+1
mFP.close()
aFP.close()

# Produce phylogenetic trees
from Bio.Phylo import TreeConstruction, draw
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
DM = TreeConstruction._DistanceMatrix(name, data2D)
tree = DistanceTreeConstructor().upgma(DM) # UPGMA
#tree = DistanceTreeConstructor().nj(DM) # Neighbor join

# Remove clade names with 'Inter#'
for clade in tree.find_clades():
    if 'Inner' in clade.name:
        clade.name = ''
if(dumpFlag == 1):
    dFP.write(str(tree))
    dFP.close()

# Print the tree to the output file
import matplotlib.pyplot as plt
import japanize_matplotlib
from Bio import Phylo
import pylab
def get_label(leaf):    # Show the entire label name
    return leaf.name

fig = plt.figure(figsize=(j/2, j/2), dpi=150)
# Adjust the size of the figure according to the number of nodes
#if(j < 30): # horizontal (landscape) view
#    fig = plt.figure(figsize=(j, j/2), dpi=150)
#else: # vertical view
#    fig = plt.figure(figsize=(j/2, j), dpi=150)
    
axes = fig.add_subplot(1, 1, 1)
plt.gca().axis('off') # Do not display the axes

draw(tree, label_func=get_label, axes=axes, do_show=False)  # Avoid showing the tree and print the tree to the output directly

fig.subplots_adjust(right=0.7) # The right edge is located at 0.7 of the entire canvas. See  http://taustation.com/pyplot-subplot-location-adujustment/

fig.savefig(oFile)


