# hEDDC

hEDDC implements an efficient heuristic algorithm to estimate the edit distance with duplication and contraction of units (EDDC, for short) between complex tandem repeats with different units (motifs).
Examples of complex tandem repeats are 

`(AAAG)i (AG)j`, `(ACC)i (GTT)j`, and `(AAAG)i (AG)j (AAAG)k`, 

where each unit in parentheses is associated with a variable (i,j, and k) that shows the number of unit occurrences.
Given a set of complex tandem repeats in a fasta format, hEDDC outputs a (normalized) distance matrix between all pairs of complex tandem repeats in the input fasta file.
The output distance matrix can be used to generate a phylogenetic tree of the input complex tandem repeats.

## Requirements

Rust(>1.59)

## Install

After installing rust language, go to the directory `eddc-heuristics` and type

`cargo run --release`

Then, you will have an executable binary:

`./eddc-heuristics/target/release/hEDDC` 

hEDDC is a heuristic algorithm for estimating the accurate EDDC, but its estimation is highly correlated with the accurate EDDC, which is shown by substantial benchmark data. 

If one wants to calcualte the exact EDDC value, go to the directory `eddc-exact` and type

`cargo run --release`

Then, you can obtain:

`./eddc-exact/target/release/EDDC_exact` ,

which outputs the exact EDDC value, though `EDDC_exact` is orders of magnitude slower than `hEDDC`.

## How to use

To call `./eddc-heuristics/target/release/hEDDC`, you need to specify two fasta file, reads (strings with complex tandem repeats) and frequent units that occur in complex tandem repeats. 
They should be both (possibly multi-line) fasta files.
An example fasta file of reads starts with:

    > (AAAG)15(AG)15
    AAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAG

and an example fasta file of units is:

    > frequent unit.
    AAAG
    > frequent unit.
    AG

Given two fasta files, call

`./eddc-heuristics/target/release/hEDDC --reads <READS> --units <UNITS>`

The normalized distance matrix bewteen complex tandem repeats in `<READS>` will be output to the stdout. 
Similarly, to call `./eddc-exact/target/release/EDDC_exact`, you need to specify arguments of the same type.
 
This binary is fully parallelized, and the number of the threads would be automatically estimated.

## Sample codes to run `hEDDC` using benchmark data

See the directory `sample-data`.

## Visualization of phylogenetic trees of complex tandem repeats

See the directory `sample-data`.

## Authors

- Riki Kawahara <4956309606@edu.k.u-tokyo.ac.jp> implemented codes in  `eddc-heuristics` and `eddc-exact`.
- Bansho Masutani <ban-m@g.ecc.u-tokyo.ac.jp> codes in `string-decomposer`.
- Shinichi Morishita <moris@edu.k.u-tokyo.ac.jp> data and codes in `sample-codes`.
