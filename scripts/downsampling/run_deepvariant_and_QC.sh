# list of bams to process
bams=$(ls /mnt/netapp2/Store_csebdjgl/lynx_genome/lynx_data/mLynLyn1.2_ref_bams/downsampling/*.bam)

# reference genome directory
ref_genome=/mnt/netapp2/Store_csebdjgl/reference_genomes/lynx_lynx_mLynLyn1.2/mLynLyn1.2.revcomp.scaffolds.fa

# sample list directory
sample_list=/mnt/netapp2/Store_csebdjgl/agonev/data/barcodes/samples_sex.txt

# output directory
output_dir=/mnt/netapp2/Store_csebdjgl/lynx_genome/lynx_data/mLynLyn1.2_ref_gvcfs/downsampling

# loop through bams and submit jobs
for bam in $bams; do
    sample_name=$(basename ${bam} | cut -d'_' -f1,2,3,4)
    echo "sample_name: ${sample_name}"
    sbatch -o logs/variant_calling/deepvariant/downsampling/${sample_name}.out -e logs/variant_calling/deepvariant/downsampling/${sample_name}.err \
    -c 32 --mem=200GB -t 12:00:00 --gres=gpu:a100:1 scripts/variant_calling_deepvariant.sh ${bam} ${ref_genome} ${sample_list} ${output_dir} ${sample_name}
done

# QC

module load samtools

for i in /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/mLynLyn1.2_ref_gvcfs/*_mLynLyn1.2_ref.vcf.gz; do
  bcftools stats ${i} > /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/mLynLyn1.2_ref_gvcfs/bcftools_stats/$(basename "${i}" .vcf.gz)_stats.txt
done

module load multiqc 
multiqc *_stats.txt

# grep all the .out files that contain the word "error"
grep -l "error" *.out > error_files.txt


bams=$(ls /mnt/netapp2/Store_csebdjgl/lynx_genome/lynx_data/mLynLyn1.2_ref_bams/downsampling/*.bam | grep -E 'c_ll_ca_0245|c_ll_ca_0248|c_ll_ca_0254')