#!/bin/bash

dir_sample_data=sample_data
dir_hEDDC=sample_data_distance_hEDDC
dir_eddc_exact=sample_data_distance_eddc_exact
dir_pdf_cor_phylo=pdf_cor_phylo

for PAT in {"AAAG_AG","AAAG_AG_AAAG","AAAG_AG_AGGG_AG_AAAG","AGGGG_AAAAGAAAGAGAGGG_AGGGG","CAG_CAA"};do

    #for ErrRate in {"0.0"}; do
    for ErrRate in {"0.0","0.01","0.03","0.05","0.1"}; do
        #echo $PAT $ErrRate
        ID=$PAT"_0_50_"$ErrRate
        echo $ID
    
        heuDistance=$dir_hEDDC/${ID}"_distance.txt"
        exactDistance=$dir_eddc_exact/${ID}"_distance.txt"
        pdfCor=$dir_pdf_cor_phylo/${ID}"_correlation.pdf"
        python3 cor.py $heuDistance $exactDistance $pdfCor
        
        fastaFile=$dir_sample_data/${ID}".fasta"
        
        pdfPhyloHeu=$dir_pdf_cor_phylo/${ID}"_phylo_heuristics.pdf"
        python3 phylo.py $heuDistance $fastaFile $pdfPhyloHeu
        
        pdfPhyloExact=$dir_pdf_cor_phylo/${ID}"_phylo_exact.pdf"
        python3 phylo.py $exactDistance $fastaFile $pdfPhyloExact
        
    done
done
exit 0
