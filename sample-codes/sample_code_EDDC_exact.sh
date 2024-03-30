#!/bin/bash

run_heuristics=../eddc-exact/target/release/EDDC_exact
dataDir=sample_data
distanceDir=sample_data_distance_EDDC_exact

MaxReads=20  # Number of reads
echo "Number of reads is"  $MaxReads

for PAT in {"AAAG_AG","ACC_GTT","AAAG_AG_AAAG","AAAG_AG_AGGG_AG_AAAG","AGGGG_AAAAGAAAGAGAGGG_AGGGG"};do

    for ErrRate in {"0.0","0.01","0.03","0.05","0.1"}; do
        ID=$PAT"_10_50_"$ErrRate
        echo $ID

        readFile=./$dataDir/$ID".fasta"
        unitFile=./$dataDir/$ID"_units.fasta"
        distanceFile=./$distanceDir/$ID"_distance.txt"

        time -p $run_heuristics --reads $readFile --units $unitFile > $distanceFile

    done
done
