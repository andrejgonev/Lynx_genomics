#!/bin/bash

# bam folder
bam_folder=/mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/mLynLyn1.2_ref_bams

# sbatch QualiMap of all samples
bams=$(ls ${bam_folder}/*er.bam)

# sbatch QualiMap of specific samples I want to check out
#bams=$(ls ${bam_folder}/*er.bam | grep -E "c_ll_ya_0140|c_ll_ya_0141|c_ll_ya_0142|c_ll_ya_0143|c_ll_ya_0145|c_ll_ya_0146|c_ll_ya_0147")

#!/bin/bash

for bam in $bams ; do

    sampleid=$(basename "$bam" | sed 's/_mLynLyn_ref.sorted.rg.merged_sorted.rmdup.indelrealigner.bam$//')
    echo "sbatch qualimap of ${sampleid}"
    sbatch -o logs/alignment/qualimap/qualimap-${sampleid}.out -e logs/alignment/qualimap/qualimap-${sampleid}.err \
    scripts/sbatch_qualimap_bam_out.sh ${bam} data/qualimap

done