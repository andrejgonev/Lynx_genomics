vcf_dir=/mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/mLynLyn1.2_ref_vcfs

populations=("bal" "cau" "crp" "lva" "nor" "pol" "tva" "mng" "wru" "yak" "pyk")

for pop in "${populations[@]}"; do
    # apply the filter based on missingness
    echo "filtering ${pop} vcf: removing SNPs with an excess of missing data"
    bedtools subtract -header \
        -a ${vcf_dir}/c_ll_105_mLynLyn1.2_ref.filter4.${pop}_pop.vcf \
        -b data/variant_filtering/missingness/${pop}.miss_filter.bed \
    > ${vcf_dir}/c_ll_105_mLynLyn1.2_ref.filter4.${pop}_pop.miss20.vcf

    # apply the filter based on read depth
    echo "filtering ${pop} vcf: removing 10k windows with an excess of read depth"
    bedtools subtract -header \
        -a ${vcf_dir}/c_ll_105_mLynLyn1.2_ref.filter4.${pop}_pop.miss20.vcf \
        -b data/variant_filtering/depth/${pop}.rd_filter.bed \
    > ${vcf_dir}/c_ll_105_mLynLyn1.2_ref.filter4.${pop}_pop.miss20.rd_fil.vcf
done