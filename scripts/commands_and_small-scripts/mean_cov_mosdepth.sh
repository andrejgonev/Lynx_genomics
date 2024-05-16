#!/bin/bash
#SBATCH -e /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/agonev/logs/mosdepth/slurm-%j.err
#SBATCH -o /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/agonev/logs/mosdepth/slurm-%j.out

# This script generates a bed file that contains the global mean depth and per chromosome from a set of BAM files. 

# Example of the mosdepth command:
#   mosdepth -t 4 -n --by 10000 <output_prefix> <input_bam>

# Usage of the script simultaneously for several bam files:
#    for input_bam in $(ls /bams/directory/*.bam); do 
#        job_id=(sbatch -c 16 -mem=10GB -t 00:15:00 mean_cov_mosdepth.sh <${input_bam}> <output_directory> | awk '{print $4}')
#        echo "${job_id} ${input_bam}" >> /logs/directory/job_ids_depth_filtering_bed_generation.txt
#    done

# In my case, I first created a bed file with all the chromosomes of the lynx genome, excluding the Y chromosome and all the small scaffolds.
# I did this by running the following command:
#grep "Chr" /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/reference_genomes/mLynLyn1.2.revcomp.scaffolds.fa.fai | grep -v "ChrY" | awk '{print $1, $2}'> /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/reference_genomes/lynx_lynx_mLynLyn1.2/mLynLyn1.2.revcomp.scaffolds.big_chromosomes_noY.temp
#awk '{print $1, 0, $2}' /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/reference_genomes/lynx_lynx_mLynLyn1.2/mLynLyn1.2.revcomp.scaffolds.big_chromosomes_noY.temp > /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/reference_genomes/lynx_lynx_mLynLyn1.2/mLynLyn1.2.revcomp.scaffolds.big_chromosomes_noY.bed
#rm /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/reference_genomes/lynx_lynx_mLynLyn1.2/mLynLyn1.2.revcomp.scaffolds.big_chromosomes_noY.temp
# change the spaces with tabs in /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/reference_genomes/lynx_lynx_mLynLyn1.2/mLynLyn1.2.revcomp.scaffolds.big_chromosomes_noY.bed using sed
# sed -i 's/ /\t/g' /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/reference_genomes/lynx_lynx_mLynLyn1.2/mLynLyn1.2.revcomp.scaffolds.big_chromosomes_noY.bed

# To Run:
# for input_bam in $(ls /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/mLynLyn1.2_ref_bams/*ner.bam); do 
#   job_id=$(sbatch -c 16 --mem=10GB -t 00:15:00 /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/agonev/scripts/mean_cov_mosdepth.sh ${input_bam} /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/mLynLyn1.2_ref_bams/mosdepth | awk '{print $4}')
#     echo "${job_id} ${input_bam}" >> /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/agonev/logs/mosdepth/job_ids_mean_cov_mosdepth.txt
# done




# Load the mosdepth module
module load mosdepth

# Define the input bam file basename:
basename_bam=$(basename ${1} .bam)

# Define the output directory:
output_dir=${2}

# Run mosdepth:
mosdepth -t 4 -n --by /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/reference_genomes/lynx_lynx_mLynLyn1.2/mLynLyn1.2.revcomp.scaffolds.big_chromosomes_noY.bed ${output_dir}/${basename_bam}_mosdepth ${1}

