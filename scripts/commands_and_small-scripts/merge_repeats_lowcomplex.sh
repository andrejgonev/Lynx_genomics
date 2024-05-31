### Option 1

ref_dir=/mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/reference_genomes/lynx_lynx_mLynLyn1.2

cat <(grep -v "#" ${ref_dir}/Repeats.4jb.gff3 | awk -F'\t' '{OFS="\t"; print $1, $4-1, $5}') \
    <(grep -v "#" ${ref_dir}/repetitive_regions/low_complex/mLynLyn1.2.revcomp.scaffolds.fa.out.gff | awk -F'\t' '{OFS="\t"; print $1, $4-1, $5}') |
    sort -k1,1 -k2,2n -k3,3n |
    bedtools merge -i - \
    > ${ref_dir}/repeats_lowcomplexity_regions.bed
    

# To calculate the length of these regions I run: 
awk '{sum+=$3-$2} END {print sum}' ${ref_dir}/repeats_lowcomplexity_regions.bed # 1047269226 

# To get the total length of the genome: 
awk '{sum+=$2} END {print sum}' ${ref_dir}/mLynLyn1.2.revcomp.scaffolds.fa.fai # 2432111198

# To get the percentage of the genome that is repetitive:
echo "scale=2; 1047269226 / 2432111198 * 100" | bc  # 43.06%


### Option 2 (Lucia)

# # generating the merged GFF
# zcat Repeats.4jb.gff3.gz > repeats_temp.gff
# cat repeats_temp.gff repetitive_regions/low_complex/mLynPar1.2.scaffolds.revcomp.scaffolds.fa.out.gff | sort -k1,1 -k4,4n > repeats_lowcomplexity_mLynPar1.2.scaffolds.revcomp.gff3
# rm repeats_temp.gff

# # generating the merged BED:
# awk -F'\t' '{OFS="\t"; print $1, $4-1, $5}' repeats_lowcomplexity_mLynPar1.2.scaffolds.revcomp.gff3 > repeats_lowcomplexity_mLynPar1.2.scaffolds.revcomp.bed

# # merging the repetitive from both files in the BED:
# bedtools merge -i repeats_lowcomplexity_mLynPar1.2.scaffolds.revcomp.bed > repeats_lowcomplexity_mLynPar1.2.scaffolds.revcomp_merged.bed

# # Calculate the total length of intersected regions:
# awk '{sum+=$3-$2} END {print sum}' repeats_lowcomplexity_mLynPar1.2.scaffolds.revcomp_merged.bed 