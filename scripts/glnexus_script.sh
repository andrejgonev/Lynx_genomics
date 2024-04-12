#!/bin/bash
#SBATCH --output=/mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/agonev/logs/glnexus/slurm-%j.out
#SBATCH --error=/mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/agonev/logs/glnexus/slurm-%j.err

# In this script, the gVCFs obtained with DeepVariant are merged into a single VCF file using glnexus.

# Usage of the script: sbatch -c 30 --mem=100G -t 03:00:00 glnexus_script.sh <gvcfs_list> <output_file>

# Load the required modules:
module load cesga/2020 glnexus/1.4.1
module load cesga/2020 bzip2/1.0.8
module load cesga/2020 samtools/1.14

# Define the list of gVCFS to be merged:
gvcfs_list=${1}

# Output file (and directory):
output=${2}

# Remove the following directory if it exists (it is created when running GLNexus, and it has to be removed before running it again):
rm -rf /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/agonev/scripts/GLnexus.DB

# Run glnexus:
glnexus_cli --config DeepVariantWGS --list ${gvcfs_list} --mem-gbytes 95 -t 30 | bcftools view - | bgzip -@ 30 -c > ${output}