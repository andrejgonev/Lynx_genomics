#!/bin/bash

# list of bams to process
bams=$(ls /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/mLynLyn1.2_ref_bams/*er.bam | grep -vE "ca_0249|ca_0253")

# reference genome directory
ref_genome=/mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/reference_genomes/lynx_lynx_mLynLyn1.2/mLynLyn1.2.revcomp.scaffolds.fa

# sample list directory
sample_list=/mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/agonev/barcodes/samples_sex.txt

# output directory
output_dir=/mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/mLynLyn1.2_ref_gvcfs/

# loop through bams and submit jobs
for bam in $bams; do
    sample_name=$(basename ${bam} _mLynLyn_ref.sorted.rg.merged_sorted.rmdup.indelrealigner.bam)
    echo "sample_name: ${sample_name}"
    sbatch -o /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/agonev/logs/variant_calling/deepvariant/${sample_name}.out -e /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/agonev/logs/variant_calling/deepvariant/${sample_name}.err \
    -c 32 --mem=64GB -t 06:00:00 --gres=gpu:a100:1 /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/agonev/scripts/variant_calling_deepvariant.sh ${bam} ${ref_genome} ${sample_list} ${output_dir} ${sample_name}
done