#!/bin/bash

run_eddc_exact=../eddc-exact/target/release/align_tandem_repeats_exact
dataDir=sample_data
distanceDir=sample_data_distance_eddc_exact

MaxReads=20  # Number of reads
echo "Number of reads is"  $MaxReads

for PAT in {"AAAG_AG","CAG_CAA","AAAG_AG_AAAG","AAAG_AG_AGGG_AG_AAAG","AGGGG_AAAAGAAAGAGAGGG_AGGGG"};do

    for ErrRate in {"0.0","0.01","0.03","0.05","0.1"}; do
        ID=$PAT"_0_50_"$ErrRate
        echo $ID

        readFile=./$dataDir/$ID".fasta"
        unitFile=./$dataDir/$ID"_units.fasta"
        distanceFile=./$distanceDir/$ID"_distance.txt"

        time -p $run_eddc_exact --reads $readFile --units $unitFile > $distanceFile

    done
done
