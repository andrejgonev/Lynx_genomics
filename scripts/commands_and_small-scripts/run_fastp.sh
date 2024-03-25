### I put this on the terminal to call Enrico's script that runs fastp. I use each row in my filtered barcodes file as a separate command, so each read is a separate slurm job

for row in $(cat barcodes/barcodes_filtered | tr ' ' ':'); do
    
    fastq_dir=$(echo $row | cut -d':' -f2)
    r1_fastq=$(echo $row | cut -d':' -f3)
    r2_fastq=$(echo $row | cut -d':' -f4)
    fastqid=$(echo $r1_fastq | sed 's/_1.fastq.gz$//')
    echo "sbatch scripts/sbatch_fastp_fqdir_r1fq_r2fq.sh ${fastq_dir} ${r1_fastq} ${r2_fastq}"
    sbatch -o logs/fastp/${fastqid}.out -e logs/fastp/${fastqid}.err scripts/sbatch_fastp_fqdir_r1fq_r2fq.sh ${fastq_dir} ${r1_fastq} ${r2_fastq}

done