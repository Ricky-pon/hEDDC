# Print a phylogenetic tree of complex tandem repeats from the distance matrix of the repeats output by Riki Kawahara's program

import sys  # argv
args = sys.argv
File1 = args[1] # distance matrix file 1
File2 = args[2] # distance matrix file 2
outFile = args[3] # outFile
FP1 = open(File1, 'r')
FP2 = open(File2, 'r')
d1 = []
d2 = []

j=0
while 1:
    tmp = FP1.readline().split(' ')
    if not tmp: break   # End of file
    if bool(tmp == ['']): break # End of file
    if(0 < len(tmp)):
        for i in range(0,min(j+1,len(tmp))): #lower triangular matrix
            if(i != j): # Ignore the value in the diagonal
                d1.append( float(tmp[i]) )
    j=j+1
FP1.close()

j=0
while 1:
    tmp = FP2.readline().split(' ')
    if not tmp: break   # End of file
    if bool(tmp == ['']): break # End of file
    if(0 < len(tmp)):
        for i in range(0,min(j+1,len(tmp))): #lower triangular matrix
            if(i != j): # Ignore the value in the diagonal
                d2.append( float(tmp[i]) )
    j=j+1
FP2.close()

#print(d1)
#print(d2)

import numpy as np

x = np.array(d1)
y = np.array(d2)
#print(np.corrcoef(x, y)[1,0]) # https://www.higashisalary.com/entry/numpy-corrcoef
corTitle = "C.C. = " + str( format(np.corrcoef(x, y)[1,0], '.5f') )

import matplotlib.pyplot as plt
import matplotlib
matplotlib.style.use("grayscale")
#matplotlib.style.use("ggplot")
plt.title(corTitle,fontsize=20)
plt.tick_params(labelsize=20)

plt.scatter(x, y, s=3)

limit = max([max(x), max(y)])
plt.xlim(0, limit)
plt.ylim(0, limit)

#plt.xlabel(File1,fontsize=5)
#plt.ylabel(File2,fontsize=5)
plt.savefig(outFile)
