vcf_dir=/mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/mLynLyn1.2_ref_vcfs
invcf=${vcf_dir}/c_ll_105_mLynLyn1.2_ref.vcf.gz
ref_dir=/mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/reference_genomes/lynx_lynx_mLynLyn1.2
ref=${ref_dir}/mLynLyn1.2.revcomp.scaffolds.fa
mask=${ref_dir}/repeats_lowcomplexity_regions.bed

sbatch --mem=12GB -t 03:00:00 scripts/variant_filter_1to4.sh ${ref} ${invcf} ${mask}



#zcat c_ll_105_mLynLyn1.2_ref.vcf.gz | grep -v "#" | wc -l
#cat c_ll_105_mLynLyn1.2_ref.filter4.vcf | grep -v "#" | wc -l