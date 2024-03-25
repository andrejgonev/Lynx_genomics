for row in $(cat data/barcodes/barcodes_filtered | tr ' ' ':'); do
    
    fastq_dir=$(echo $row | cut -d':' -f2)
    r1_fastq=$(echo $row | cut -d':' -f3)
    r2_fastq=$(echo $row | cut -d':' -f4)
    if [[ ${r1_fastq} == *.fastq.gz ]]; then
    r1_fastp=${r1_fastq/.fastq.gz/.fastp.fastq.gz}
    r2_fastp=${r2_fastq/.fastq.gz/.fastp.fastq.gz}
    elif [[ ${r1_fastq} == *.fq.gz ]]; then
    r1_fastp=${r1_fastq/.fq.gz/.fastp.fastq.gz}
    r2_fastp=${r2_fastq/.fq.gz/.fastp.fastq.gz}
    else
    echo "Error: fastq files should have extension .fastq.gz or .fq.gz"
    exit 1
    fi
    rm ${fastq_dir}/fastp/${r1_fastp} ${fastq_dir}/fastp/${r2_fastp} ${fastq_dir}/fastp/${r1_fastq/.fastq.gz/_fastp.html} \
       ${fastq_dir}/fastp/${r1_fastq/.fastq.gz/_fastp.json} ${fastq_dir}/fastp/${r1_fastq/.fastq.gz/_unpaired.fastq.gz} \
       ${fastq_dir}/fastp/${r2_fastq/.fastq.gz/_unpaired.fastq.gz} ${fastq_dir}/fastp/${r1_fastq/.fastq.gz/_failed.fastq.gz}
done