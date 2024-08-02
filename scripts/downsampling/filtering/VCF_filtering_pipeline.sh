ls /mnt/netapp2/Store_csebdjgl/lynx_genome/lynx_data/mLynLyn1.2_ref_gvcfs/downsampling/*g.vcf.gz > /mnt/netapp2/Store_csebdjgl/lynx_genome/lynx_data/mLynLyn1.2_ref_gvcfs/downsampling/gvcfs_list.txt

sbatch -c 32 --mem=100G -t 03:00:00 scripts/glnexus_script.sh /mnt/netapp2/Store_csebdjgl/lynx_genome/lynx_data/mLynLyn1.2_ref_gvcfs/downsampling/gvcfs_list.txt /mnt/netapp2/Store_csebdjgl/lynx_genome/lynx_data/mLynLyn1.2_ref_vcfs/downsampling/c_ll_mLynLyn1.2_ref_downsampled.vcf.gz 


vcf_dir=/mnt/netapp2/Store_csebdjgl/lynx_genome/lynx_data/mLynLyn1.2_ref_vcfs/downsampling
invcf=${vcf_dir}/c_ll_mLynLyn1.2_ref_downsampled.vcf
ref_dir=/mnt/netapp2/Store_csebdjgl/reference_genomes/lynx_lynx_mLynLyn1.2
ref=${ref_dir}/mLynLyn1.2.revcomp.scaffolds.fa
mask=${ref_dir}/repeats_lowcomplexity_regions.bed

sbatch --mem=12GB -t 03:00:00 scripts/variant_filter_1to4.sh ${ref} ${invcf} ${mask}



###### split vcf

#!/bin/bash
#SBATCH --time=1:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=16

ref_dir=/mnt/netapp2/Store_csebdjgl/reference_genomes/lynx_lynx_mLynLyn1.2
ref=${ref_dir}/mLynLyn1.2.revcomp.scaffolds.fa
vcf_dir=/mnt/netapp2/Store_csebdjgl/lynx_genome/lynx_data/mLynLyn1.2_ref_vcfs/downsampling/
invcf=${vcf_dir}/c_ll_mLynLyn1.2_ref_downsampled.filter4.vcf


## Create an index file for the input vcf if it doesn't exist
 gatk IndexFeatureFile \
     -I ${invcf}

populations=("bal" "cau" "crp" "kir" "nor" "pol")

for pop in "${populations[@]}"; do
    if [ ${pop} == "bal" ]; then
        samples=($(grep "ll_ba" data/downsampling/sample.list))
    elif [ ${pop} == "cau" ]; then
        samples=($(grep "ll_ca" data/downsampling/sample.list | grep -vE "ca_0245|ca_0248"))
    elif [ ${pop} == "crp" ]; then
        samples=($(grep "ll_cr" data/downsampling/sample.list))
    elif [ ${pop} == "kir" ]; then
        samples=($(grep -E "ll_ki" data/downsampling/sample.list))
    elif [ ${pop} == "nor" ]; then
        samples=($(grep "ll_no" data/downsampling/sample.list))
    elif [ ${pop} == "pol" ]; then
        samples=($(grep "ll_po" data/downsampling/sample.list))        
    fi
    echo "-- creating vcf of ${pop} --"
    gatk SelectVariants \
        -R ${ref} \
        -V ${invcf} \
        $(for sample in ${samples[@]}; do echo "-sn ${sample}";done) \
        -O ${invcf/.vcf/.${pop}_pop.vcf}
done


## Depth filters

ls /mnt/netapp2/Store_csebdjgl/lynx_genome/lynx_data/mLynLyn1.2_ref_bams/downsampling/*.bam | grep -vE "ca_0249|ca_0253" > path/to/folder/agonev/data/c_ll_105.bamlist 

inbams=($(cat data/c_ll_105.bamlist))
for bam in ${inbams[*]}; do
    echo "calculating depth of $bam"
    sbatch scripts/sbatch_mosdepth_10k_bam_outdir.sh \
        ${bam} data/variant_filtering/depth
done




# Missingness
sbatch -o logs/variant_filtering/missing/missingness_beds.out -e logs/variant_filtering/missing/missingness_beds.err scripts/missingness_bed_table_downsampled.sh

#not done yet!



## Apply filters

vcf_dir=/mnt/netapp2/Store_csebdjgl/lynx_genome/lynx_data/mLynLyn1.2_ref_vcfs/downsampling/
invcf=${vcf_dir}/c_ll_mLynLyn1.2_ref_downsampled.filter4.vcf

    # apply the filter based on read depth
    echo "filtering ${invcf} : removing 10k windows with an excess of read depth"
    bedtools subtract -header \
        -a ${invcf} \
        -b data/downsampling/variant_filtering/depth/bal.rd_filter.bed | 
    bedtools subtract -header \
        -a stdin \
        -b data/downsampling/variant_filtering/depth/cau.rd_filter.bed | 
    bedtools subtract -header \
        -a stdin \
        -b data/downsampling/variant_filtering/depth/crp.rd_filter.bed |
    bedtools subtract -header \
        -a stdin \
        -b data/downsampling/variant_filtering/depth/kir.rd_filter.bed |
    bedtools subtract -header \
        -a stdin \
        -b data/downsampling/variant_filtering/depth/nor.rd_filter.bed |
    bedtools subtract -header \
        -a stdin \
        -b data/downsampling/variant_filtering/depth/pol.rd_filter.bed \
    > ${vcf_dir}/c_ll_mLynLyn1.2_ref_downsampled.filter4.rd_fil.vcf


# after applying the filter on the filter4.vcf I went from 2,661,404 to 2,643,151 (lost 18,253). I checked it with Enrico and it seems ok (we checked the printed summary of the number of lost windows). 
#* I lost these snps because I removed the 10k windows that had a higher depth than 1.5 times the mode of the summed depth distribution per population. This is the amount of snps that were present in the windows that were flagged for removal.

# i applied a general bcftools filter -i 'F_MISSING < 0.25 to the rd filtered vcf and lost 6,737 snps (from 2,643,151 to 2,636,414). I deleted this one, and will apply the filter for each population separately.
# I ran the missingness_beds.sh script and got the missingness bed files for each population. I will now apply the filter for each population separately. I will remove 30% missing data for each population.

vcf_dir=/mnt/netapp2/Store_csebdjgl/lynx_genome/lynx_data/mLynLyn1.2_ref_vcfs/downsampling/
invcf=${vcf_dir}/c_ll_mLynLyn1.2_ref_downsampled.filter4.rd_fil.vcf

    # apply the filter based on missingness
    echo "filtering ${invcf} : removing 30% missing data for each population"
    bedtools subtract -header \
        -a ${invcf} \
        -b data/downsampling/variant_filtering/missing/bal.f_miss.0_3.bed | 
    bedtools subtract -header \
        -a stdin \
        -b data/downsampling/variant_filtering/missing/cau.f_miss.0_3.bed | 
    bedtools subtract -header \
        -a stdin \
        -b data/downsampling/variant_filtering/missing/crp.f_miss.0_3.bed |
    bedtools subtract -header \
        -a stdin \
        -b data/downsampling/variant_filtering/missing/kir.f_miss.0_3.bed |
    bedtools subtract -header \
        -a stdin \
        -b data/downsampling/variant_filtering/missing/nor.f_miss.0_3.bed |
    bedtools subtract -header \
        -a stdin \
        -b data/downsampling/variant_filtering/missing/pol.f_miss.0_3.bed \
    > ${vcf_dir}/c_ll_mLynLyn1.2_ref_downsampled.filter4.rd_fil.miss30.vcf


# if I use this method for 25%, I go from 2,643,151 to 2,365,073 snps. I would lose 278,078 snps.
# I decited to apply a 30% missingness filter. After applying, I got to 2,492,033 snps. I lost 151,118 snps.


# finally i removed sites with AF=0.00

bcftools view -e "INFO/AF=0.00" ${vcf_dir}/c_ll_mLynLyn1.2_ref_downsampled.filter4.rd_fil.miss30.vcf > ${vcf_dir}/c_ll_mLynLyn1.2_ref_downsampled.filter4.rd_fil.miss30.variant.vcf

# Final number of SNPs: 2492029. I lost 4 SNPs with this. 
#To make a note for the future, it would be better to eliminate these sites right after I remove the half calls and recalculate the INFO fields.