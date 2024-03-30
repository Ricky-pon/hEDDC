#!/bin/bash

dirFasta=sample_data
dirDistance_hEDDC=sample_data_distance_hEDDC
dirDistance_EDDC_exact=sample_data_distance_EDDC_exact
dir_pdf_cor_phylo=phylo_cor_PDF

for PAT in {"AAAG_AG","ACC_GTT","AAAG_AG_AAAG","AAAG_AG_AGGG_AG_AAAG","AGGGG_AAAAGAAAGAGAGGG_AGGGG"};do

    for ErrRate in {"0.0","0.01","0.03","0.05","0.1"}; do
        #echo $PAT $ErrRate
        ID=$PAT"_10_50_"$ErrRate
        echo $ID
    
        heuDistance=$dirDistance_hEDDC/${ID}"_distance.txt"
        exactDistance=$dirDistance_EDDC_exact/${ID}"_distance.txt"
        pdfCor=$dir_pdf_cor_phylo/${ID}"_correlation.pdf"
        python3 python_cor.py $heuDistance $exactDistance $pdfCor
        
        fastaFile=$dirFasta/${ID}".fasta"
        
        pdfPhyloHeu=$dir_pdf_cor_phylo/${ID}"_phylo_heuristics.pdf"
        python3 python_phylo.py $heuDistance $fastaFile $pdfPhyloHeu
        
        pdfPhyloExact=$dir_pdf_cor_phylo/${ID}"_phylo_exact.pdf"
        python3 python_phylo.py $exactDistance $fastaFile $pdfPhyloExact
        
    done
done

exit 0
