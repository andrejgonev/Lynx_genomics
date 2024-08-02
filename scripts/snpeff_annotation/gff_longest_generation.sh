#!/bin/bash

while IFS= read -r line; do
    # Check if the third column matches the word "gene"
    if [[ $(echo "$line" | awk '{print $3}') == "gene" ]]; then
        echo "$line" >> snpeff_annotation/LYLY1_2A.FA_longestisoform.gff3
        GENE=$(echo "$line" | cut -f1 -d';' | cut -f9 | sed 's/ID=//')
        TRANSCRIPT=$(grep "${GENE}" snpeff_annotation/LYLY1_2A.transcripts_longest.fa | sed 's/>//')
        grep -w "${TRANSCRIPT}" LYLY1_2A.FA.gff3 >> snpeff_annotation/LYLY1_2A.FA_longestisoform.gff3
    else
        echo "Skypping this line, not a gene"
    fi
done < LYLY1_2A.FA.gff3
