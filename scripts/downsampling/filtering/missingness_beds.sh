#!/bin/bash
#SBATCH --time=6:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=24

module load bcftools

vcf_dir=/mnt/netapp2/Store_csebdjgl/lynx_genome/lynx_data/mLynLyn1.2_ref_vcfs/downsampling/
data_dir=/mnt/netapp2/Store_csebdjgl/agonev/data/downsampling

# total number of SNPs is the same in each population
total_snps=$(grep -v "#" ${vcf_dir}/c_ll_mLynLyn1.2_ref_downsampled.filter4.vcf | wc -l)

# Populations
populations=("bal" "cau" "crp" "kir" "nor" "pol")

# create a table for plotting
echo "population total_snps f_miss filtered_snps" > ${data_dir}/variant_filtering/missing/miss_table.txt

for pop in "${populations[@]}"; do
    vcf=${vcf_dir}/c_ll_mLynLyn1.2_ref_downsampled.filter4.${pop}_pop.vcf
    # work with different proportions of missing data
    for f_miss in 0.1 0.15 0.2 0.25 0.3; do
        echo "extracting SNPs in ${pop} with F_MISSING > ${f_miss}"

        ## to extract the bed files:
        # change f_miss string for better file name
        f_miss_x=$(echo $f_miss | tr '.' '_')
        # extract bed of SNPs with more missing data than proportion defined by f_miss
        bcftools filter -i "F_MISSING > ${f_miss}" ${vcf} | grep -v "#" | awk '{print $1, $2-1, $2}' | tr ' ' '\t' \
            > ${data_dir}/variant_filtering/missing/${pop}.f_miss.${f_miss_x}.bed
        
        ## to fill the table for plotting:
        # count the number of SNPs that have more missing data than proportion
        filtered_snps=$(cat ${data_dir}/variant_filtering/missing/${pop}.f_miss.${f_miss_x}.bed | wc -l)
        # add a row to table for plotting
        echo "${pop} ${total_snps} ${f_miss} ${filtered_snps}" >> ${data_dir}/variant_filtering/missing/miss_table.txt
        
        ## to print a summary:
        if [ ${pop} == "bal" ]; then
            nsamples=$(grep "ll_ba" ${data_dir}/sample.list | wc -l)
        elif [ ${pop} == "cau" ]; then
            nsamples=$(grep -E "ll_ca" ${data_dir}/sample.list | grep -vE "ca_0245|ca_0248" | wc -l)
        elif [ ${pop} == "crp" ]; then
            nsamples=$(grep -E "ll_cr" ${data_dir}/sample.list | wc -l)
        elif [ ${pop} == "kir" ]; then
            nsamples=$(grep -E "ll_ki" ${data_dir}/sample.list | wc -l)
        elif [ ${pop} == "nor" ]; then
            nsamples=$(grep -E "ll_no" ${data_dir}/sample.list | wc -l)
        elif [ ${pop} == "pol" ]; then
            nsamples=$(grep -E "ll_po" ${data_dir}/sample.list | wc -l)

        fi

        n_miss=$(echo "${nsamples} * ${f_miss}" | bc -l)
        perc_left=$(echo "${filtered_snps} / ${total_snps} * 100" | bc -l | cut -c1-5)
        echo "${pop}: ${filtered_snps} SNPs have a proportion of missing genotypes greater than ${f_miss} (>${n_miss} samples) = ${perc_left}% of total SNPs filtered"
    done
done