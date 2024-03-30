# Benchmark data

We use bench mark datasets with five complex tandem repeat (TR, for short) patterns: 

- (AAAG)i (AG)j
- (ACC)i (GTT)j  
- (AAAG)i (AG)j (AAAG)k
- (AAAG)i (AG)j, (AAAG)k (AG)l (AAAG)m  
- (AGGGG)i (AAAAGAAAGAGAGGG)j (AGGGG)k

Individual complex TR pattern has two or more units, and each unit in parentheses is associated with a variable (i,j,k,l and m) that shows the number of unit occurrences.

To generate a variety of complex TRs, random values ranging from 10 to 50 are assigned to each variable independently.
We also changed the characters in the complex TR by incurring sequence errors (substitutions, insertions, and deletions) at rates of 0%, 1%, 3%, 5%, and 10%.
For each of the five complex TR patterns, five datasets with different error rates are generated, so a total of 25 datasets are used.

Directory `sample_data` has these 25 datasets.
Each file in the directory contains 20 tandem repeats (TRs), because the program named `EDDC_exact` is very slow in calculating distances and can process 20 complex TRs in a reasonable time. 

The name of each file denotes units used in TRs, the range of unit occurrences in TRs, and the error rate. 
For example,`AAAG_AG_10_50_0.01.fasta` implies that TRs consist of units AAAG and AG, the number of unit occurrences ranges from 10 to 50, and the error rate is set to 0.01 (1\%).

# Sample codes to run `hEDDC`  

Recall that we run `hEDDC` by calling

`./eddc-heuristics/target/release/hEDDC --reads <READS> --units <UNITS>`

Then, `hEDDC` estimates the edit distance with duplication and contraction (EDDC, for short) between all pairs of complex TRs, say S and T, in the input file `<READS>`.
Since each dataset has different complex TRs with different lengths, the distance between two complex TRs largely depends on their lengths.
To facilitate comparison, `hEDDC` normalizes the EDDC by dividing it by the square root of the product of the lengths of S and T.

`sample_code_hEDDC.sh` is a sample bash script that feeds 25 datasets in the directory `sample_data` and outputs the matrix of the estimated EDDCs between 20 TRs in each of 25 datasets into the directory `sample_data_distance_hEDDC`. 
For example, 

`sample_data_distance_hEDDC/AAAG_AG_10_50_0.01_distance.txt` 

describes the 20x20 matrix of the normalized EDDCs between 20 TRs in `sample_data/AAAG_AG_10_50_0.01.fasta`.

Similarly, for computing exact EDDC values, `sample_code_EDDC_exact.sh` is a sample bash script that feeds 25 datasets in the directory `sample_data` and outputs the matrix of the *exact* EDDC values between 20 TRs in each of 25 datasets into the directory `sample_data_distance_EDDC_exact`. 


# Visualization of dot plots and phylogenetic trees

### Dot plots of exact and estimated distances and Pearson correlation coefficient

For each of the 25 datasets, script `python_cor.py` calculates the Pearson correlation coefficient of the normalized distances in the directories `sample_data_distance_hEDDC` and `sample_data_distance_EDDC_exact` that are respectively output by `hEDDC` and `EDDC_exact` in order to measure how the two algorithms are consistent.

`phylo_cor.sh` calls `python_cor.py` and outputs the PDF of the dot plot with Pearson correlation coefficient for each of the 25 datasets into the directory `phylo_cor_PDF`.

For example, for the data in `AAAG_AG_10_50_0.01.fasta`, the dot plot in `AAAG_AG_10_50_0.01_correlation.pdf` is generated. 

In the dot plot, a dot shows a pair of the normalized distances calculated by `hEDDC` in the x-axis and `EDDC_exact` in the y-axis. The Pearson correlation coefficient is shown at the top of each plot.

### Phylogenetic trees of complex TRs

From the distance matrix of normalized estimated EDDC (exact EDDC, respectively) output by `hEDDC` (`EDDC_exact`) for each dataset, script `python_phylo.py` generates a phylogenetic tree of complex TRs using unweighted pair group method with arithmetic mean (UPGMA).

`phylo_cor.sh` also calls `python_phylo.py` and outputs the PDF of the phylogenetic trees from the distance matrix output by `hEDDC` and `EDDC_exact` for each of the 25 datasets into the directory `phylo_cor_PDF`.

For example, the phylogenetic tree shown in `AAAG_AG_10_50_0.01_phylo_heuristics.pdf` is calcuated from the distance matrix in `sample_data_distance_hEDDC/AAAG_AG_10_50_0.01_distance.txt`.

Similarly, the tree in `AAAG_AG_10_50_0.01_phylo_exact.pdf` is from `sample_data_distance_EDDC_exact/AAG_AG_10_50_0.01_distance.txt`.


### Scripts for generating dot plots and phylogenetic (UPGMA) trees

We generated the dot plots and phylogenetic trees using bash script `phylo_cor.sh` that calls python scripts `python_cor.py` and `python_phylo.py`.
