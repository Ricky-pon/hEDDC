# Benchmark data

We use a bench mark data set with five complex tandem repeat (TR, for short) patterns: 

- (AAAG)i (AG)j
- (CAG)i (CAA)j  
- (AAAG)i (AG)j (AAAG)k
- (AAAG)i (AG)j, (AAAG)k (AG)l (AAAG)m  
- (AGGGG)i (AAAAGAAAGAGAGGG)j (AGGGG)k

Individual complex TR pattern has two or more units, and each unit in parentheses is associated with a variable (i,j,k,l and m) that shows the number of unit occurrences.

To generate a variety of complex TRs, random values ranging from 0 to 50 are assigned to each variable independently.
We also modified characters in complex TRs by sequencing errors (substitutions, insertions, and deletions) at the rate of 0\%, 1\%, 3\%, 5\%, and 10\%.
For each of the five complex TR patterns, five datasets with different error rates are generated, so a total of 25 datasets are used.

## Directory: input\_fasta\_files

Each files in the directory contains 20 tandem repeats (TRs), because the exact EDDC algorithm below is very slow in calculating distances and can process 20 complex TRs in a reasonable time. 
The file name denotes units used in TRs, the range of unit occurrences in TRs, and the error rate. 

For example, **AAAG_AG_0_50_0.01.fasta** 
implies that TRs consist of units AAAG and AG, the number of unit occurrences ranges from 0 to 50, and the error rate is set to 0.01 (1\%).

## Directory: eddc\_exact

Let EDDC denote the edit distance with duplication and contraction between a pair of complex TRs S and T.
Since each dataset has different complex TRs with different lengths, the distance between two complex TRs largely depends on their lengths.
To facilitate comparison, the EDDC is normalized by dividing it by the square root of the product of the lengths of S and T.
Each file in the directory describes the matrix of the normalized EDDCs between 20 TRs in the input fasta file.

For example, **AAAG_AG_0_50_0.0_distance.txt** describes the 20x20 matrix of the normalized EDDCs between 20 TRs in **AAAG_AG_0_50_0.01.fasta**.

## Directory: eddc\_heuristics

Using the EDDC heuristic algorithm, we estimate the EDDC between a pair of complex TRs. 
Each file in the directory shows the matrix of the estimated EDDCs between 20 TRs in the input fasta file. 
The estimated values are also normalized by the square root of the product of the lengths of S and T.

## Directory output\_pdf\_files

### Dot plots and Pearson correlation coefficient

For each of the 25 datasets, the Pearson correlation coefficient of the normalized distances output by the exact and heuristic EDDC algorithms is used to measure how the two algorithms are consistent.

For example, for the data in **AAAG_AG_0_50_0.01.fasta**, the dot plot in **AAAG_AG_0_50_0.0_cor.pdf**, in which "cor" is appended to the file name, is generated. 

In the dot plot, a dot shows a pair of the normalized distances calculated by the heuristic EDDC algorithm in the x-axis and the exact EDDC algorithm in the y-axis. The Pearson correlation coefficient is shown at the top of each plot.

### Phylogenetic trees of complex TRs

According to the distance matrix of normalized (estimated, respectively) EDDC output by the exact (heuristic) algorithm for each dataset, a phylogenetic tree of complex TRs is generated using unweighted pair group method with arithmetic mean (UPGMA).

For example, the phylogenetic tree shown in **AAAG_AG_0_50_0.0_phylo_exact.pdf** is calcuated from the distance matrix in **AAAG_AG_0_50_0.0_distance.txt** in directory **eddc\_exact**.

Similarly, the tree in **AAAG_AG_0_50_0.0_phylo_heu.pdf** is from **AAAG_AG_0_50_0.0_distance.txt** in directory **eddc\_heuristics**.


## Scripts for generating dot plots and phylogenetic (UPGMA) trees

We generated the dot plots and phylogenetic trees using bash script main.sh that calls python scripts cor.py and phylo.py.
