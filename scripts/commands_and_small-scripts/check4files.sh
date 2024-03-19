for row in $(cat Barcodes/barcodes_filtered | tr ' ' ':'); do
    
    fastq_dir=$(echo $row | cut -d':' -f2)
    r1_fastq=$(echo $row | cut -d':' -f3)
    r2_fastq=$(echo $row | cut -d':' -f4)

ls $fastq_dir/$r1_fastq
ls $fastq_dir/$r2_fastq

done