# INBREEDING

# First we convert the VCF file to GEN file format
module load bcftools
bcftools convert -t ^mLynLyn1.2_ChrX -g c_ll_mLynLyn1.2_ref_downsampled --3N6 --tag PL c_ll_mLynLyn1.2_ref_downsampled.filter4.rd.miss30_polarized.vcf.gz
# 2310011 records written, 0 skipped: 0/0/0/0 no-ALT/non-biallelic/filtered/duplicated
# I have less SNPs here than in the original VCF file, because I removed the X chromosome


# remove the third column from the GEN file, keep all the rest
gunzip -c c_ll_mLynLyn1.2_ref_downsampled.gen.gz
cut -d' ' --complement -f3 c_ll_mLynLyn1.2_ref_downsampled.gen > c_ll_mLynLyn1.2_ref_downsampled_polarized.gen
rm c_ll_mLynLyn1.2_ref_downsampled.gen
gzip c_ll_mLynLyn1.2_ref_downsampled_polarized.gen

# RZooRoH computes allele frequencies from the gen file, and if we don't separate the gen file by population, it will compute the inbreeding coefficient for the whole dataset
# Split vcf file by population
# I used the split_vcf.sh script to split the polarized vcf file by population (check downsampling/filtering/VCF_filtering_pipeline.sh)

# Now I convert the VCF files to GEN files
module load bcftools
for pop in bal cau crp kir nor pol; do
    bcftools convert -t ^mLynLyn1.2_ChrX -g temp."$pop" --3N6 --tag PL c_ll_mLynLyn1.2_ref_downsampled.filter4.rd.miss30_polarized."$pop"_pop.vcf
done

# Remove the third column from the GEN files. Keep the first two lines and keep only the first column from the sample files
for pop in bal cau crp kir nor pol; do
    gunzip temp."$pop".gen.gz
    cut -d' ' --complement -f3 temp."$pop".gen > c_ll_mLynLyn1.2_ref_downsampled_polarized."$pop".gen
    tail -n +3 temp."$pop".sample | cut -d' ' -f1 > c_ll_mLynLyn1.2_ref_downsampled_polarized."$pop".sample
done

rm temp.*