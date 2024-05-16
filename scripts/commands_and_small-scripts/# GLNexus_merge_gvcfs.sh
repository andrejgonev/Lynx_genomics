# Check quality of gvcfs

module load samtools

for i in /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/mLynLyn1.2_ref_gvcfs/*_mLynLyn1.2_ref.vcf.gz; do
  bcftools stats ${i} > /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/mLynLyn1.2_ref_gvcfs/bcftools_stats/$(basename "${i}" .vcf.gz)_stats.txt
done

module load multiqc 
multiqc *_stats.txt

## Merge gvcfs with GLNexus

#mkdir -p /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/mLynLyn1.2_ref_vcfs

ls /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/mLynLyn1.2_ref_gvcfs/*g.vcf.gz > /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/mLynLyn1.2_ref_gvcfs/gvcfs_list.txt

sbatch -c 32 --mem=100G -t 03:00:00 /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/agonev/scripts/glnexus_script.sh /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/mLynLyn1.2_ref_gvcfs/gvcfs_list.txt /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/mLynLyn1.2_ref_vcfs/c_ll_105_mLynLyn1.2_ref.vcf.gz 
