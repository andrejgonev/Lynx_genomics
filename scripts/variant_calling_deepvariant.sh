#!/bin/bash

# This is a slightly modified version of Lucia's script for variant calling using DeepVariant.

# This script calls variants from a set of bam files to obtain gVCFs with DeepVariant, using the model WGS.

# Example of the deepvariant command:
#   deepvariant_gpu --model_type=WGS --ref=/path/to/reference/genome \
#   --reads=/path/to/bam --output_gvcf=/path/to/output/output.g.vcf.gz --num_shards=32 --output_vcf=/path/to/output/output.vcf.gz

# Usage of the script simultaneously for several bam files:
#    for input_bam in $(ls /bams/directory/*.bam); do 
#        job_id=(sbatch -c 32 --mem=50GB -t 06:00:00 --cpus-per-task=32 --gres=gpu:a100:1 variant_calling_deepvariant.sh \
#               <${input_bam}> <ref_genome> <files_list> <output_directory> | awk '{print $4}')
#        echo "${job_id} ${input_bam}" >> /logs/directory/job_ids_deepvariant.txt
#    done

# Load the deepvariant module:
module load cesga/2020 deepvariant/1.6.0

# Define the input bam file:
bam=${1}

# Define reference genome:
ref_genome=${2}

# Define list of samples (FASTQ files):
list_file=${3}

# Define the output directory:
output_dir=${4}

# Define sample name of bam:
sample_name=${5}
echo "sample_name: ${sample_name}"

# Define sample sex:
sex=$(grep ${sample_name} ${list_file} | cut -f2)
echo "sample_sex: ${sex}"

# Run deepvariant (changing the command according to the individuals sex):
if [[ ${sex} == "male" ]]; then
    deepvariant_gpu --model_type=WGS \
                    --ref=${ref_genome} \
                    --reads=${bam} \
                    --output_gvcf=${output_dir}/${sample_name}mLynLyn1.2_ref.g.vcf.gz \
                    --output_vcf=${output_dir}/${sample_name}_mLynLyn1.2_ref.vcf.gz \
                    --par_regions_bed="/mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/reference_genomes/lynx_pardinus_mLynLyn1.2/mLynLyn1.2.PAR1_sexChr.bed" \
                    --haploid_contigs="mLynLyn1.2_ChrX,mLynLyn1.2_ChrY,mLynLyn1.2_ChrY_unloc_*" \
                    --num_shards=32
else
    deepvariant_gpu --model_type=WGS \
                    --ref=${ref_genome} \
                    --reads=${bam} \
                    --output_gvcf=${output_dir}/${sample_name}_mLynLyn1.2_ref.g.vcf.gz \
                    --output_vcf=${output_dir}/${sample_name}_mLynLyn1.2_ref.vcf.gz \
                    --num_shards=32
fi