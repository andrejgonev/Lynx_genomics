#!/bin/bash
#SBATCH --time=1:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=16

ref_dir=//mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/reference_genomes/lynx_lynx_mLynLyn1.2
ref=${ref_dir}/mLynLyn1.2.revcomp.scaffolds.fa
vcf_dir=/mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/mLynLyn1.2_ref_vcfs
invcf=${vcf_dir}/c_ll_105_mLynLyn1.2_ref.filter4.vcf


## Create an index file for the input vcf if it doesn't exist
# gatk IndexFeatureFile \
#     -I ${invcf}

populations=("bal" "cau" "crp" "lva" "nor" "pol" "tva" "mng" "wru" "yak" "pyk")

for pop in "${populations[@]}"; do
    if [ ${pop} == "bal" ]; then
        samples=($(grep "ll_ba" data/sample.list))
    elif [ ${pop} == "cau" ]; then
        samples=($(grep "ll_ca" data/sample.list))
    elif [ ${pop} == "crp" ]; then
        samples=($(grep "ll_cr" data/sample.list))
    elif [ ${pop} == "lva" ]; then
        samples=($(grep "ll_la" data/sample.list))
    elif [ ${pop} == "nor" ]; then
        samples=($(grep "ll_no" data/sample.list))
    elif [ ${pop} == "pol" ]; then
        samples=($(grep "ll_po" data/sample.list))
    elif [ ${pop} == "tva" ]; then
        samples=($(grep "ll_tu" data/sample.list))
    elif [ ${pop} == "mng" ]; then
        samples=($(grep -E "ll_ka|ll_og|ll_to" data/sample.list))
    elif [ ${pop} == "wru" ]; then
        samples=($(grep -E "ll_ki|ll_ur" data/sample.list))
    elif [ ${pop} == "yak" ]; then
        samples=($(grep "ll_ya" data/sample.list))
    elif [ ${pop} == "pyk" ]; then
        samples=($(grep "ll_vl" data/sample.list))
    fi
    echo "-- creating vcf of ${pop} --"
    gatk SelectVariants \
        -R ${ref} \
        -V ${invcf} \
        $(for sample in ${samples[@]}; do echo "-sn ${sample}";done) \
        -O ${invcf/.vcf/.${pop}_pop.vcf}
done