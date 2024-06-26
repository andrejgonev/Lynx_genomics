## Synteny

## Fasta and BED generation

# Define general paths
REF=/GRUPOS/grupolince/reference_genomes/lynx_lynx_mLynLyn1.2/mLynLyn1.2.revcomp.scaffolds.fa
CHR=($(cat /GRUPOS/grupolince/reference_genomes/lynx_lynx_mLynLyn1.2/mLynLyn1.2.revcomp.scaffolds.fa.fai | cut -f 1 | uniq))
#VCF=/GRUPOS/grupolince/mLynPar1.2_ref_vcfs/c_lp_all_novogene_sept23_mLynPar1.2_ref.filter5_QUAL20_rd.miss.vcf.gz

############################ RUFUS ############################
# Define paths
INBAM=/GRUPOS/grupolince/mLynLyn1.2_ref_bams/polarization/lr1_mLynLyn1.2_ref_sorted_rg_merged_sorted_rmdup_indelrealigner.bam
OUTFASTA=/GRUPOS/grupolince/agonev/polarization/fasta/lr1_mLynLyn1.2_ref_sorted_rg_merged_sorted_rmdup_indelrealigner_pu2fa.fa
INTERSECTFASTA=/GRUPOS/grupolince/agonev/polarization/lr1_mLynLyn1.2_ref_intersect.fa
INTERSECTBED=/GRUPOS/grupolince/mLynLyn1.2_ref_bams/polarization/lr1_mLynLyn1.2_ref_intersect.bed

# get fasta from bam per chr
rm $OUTFASTA

for i in ${CHR[@]:0:20}; do
  echo "getting FASTA for chromosome $i"
  samtools mpileup -s -q30 -f $REF $INBAM -r $i | /GRUPOS/grupolince/Chrom-Compare/pu2fa -c $i -C 100 >> $OUTFASTA
done

echo "getting intersection FASTAs"

bedtools getfasta -fi $OUTFASTA -bed $VCF -fo $INTERSECTFASTA # we only get 16257 N (Lorena had 222076)
bedtools getfasta -fi $OUTFASTA -bed $VCF -fo $INTERSECTBED -tab


############################ MARBLED CAT ############################
# Define paths
INBAM=/GRUPOS/grupolince/mLynLyn1.2_ref_bams/polarization/pm_mLynLyn1.2_ref_sorted_rg_merged_sorted_rmdup_indelrealigner.bam
OUTFASTA=/GRUPOS/grupolince/agonev/polarization/fasta/pm_mLynLyn1.2_ref_sorted_rg_merged_sorted_rmdup_indelrealigner_pu2fa.fa
INTERSECTFASTA=/GRUPOS/grupolince/agonev/polarization/pm_mLynLyn1.2_ref_intersect.fa
INTERSECTBED=/GRUPOS/grupolince/agonev/polarization/pm_mLynLyn1.2_ref_intersect.bed

# get fasta from bam per chr
rm $OUTFASTA

for i in ${CHR[@]:0:20}; do
  echo "$i"
  samtools mpileup -s -q30 -f $REF $INBAM -r $i | /GRUPOS/grupolince/Chrom-Compare/pu2fa -c $i -C 100 >> $OUTFASTA
done

bedtools getfasta -fi $OUTFASTA -bed $VCF -fo $INTERSECTFASTA    # we only get 30736 N (Lorena had 222076)
bedtools getfasta -fi $OUTFASTA -bed $VCF -fo $INTERSECTBED -tab


############################ CAT ############################
# Define paths
INBAM=/GRUPOS/grupolince/mLynLyn1.2_ref_bams/polarization/fc_mLynLyn1.2_ref_sorted_rg_merged_sorted_rmdup_indelrealigner.bam
OUTFASTA=/GRUPOS/grupolince/agonev/polarization/fasta/fc_mLynLyn1.2_ref_sorted_rg_merged_sorted_rmdup_indelrealigner_pu2fa.fa
INTERSECTFASTA=/GRUPOS/grupolince/agonev/polarization/fc_mLynLyn1.2_ref_intersect.fa
INTERSECTBED=/GRUPOS/grupolince/agonev/polarization/fc_mLynLyn1.2_ref_intersect.bed

# get fasta from bam per chr
rm $OUTFASTA

for i in ${CHR[@]:0:20}; do
  echo "$i"
  samtools mpileup -s -q30 -f $REF $INBAM -r $i | /GRUPOS/grupolince/Chrom-Compare/pu2fa -c $i -C 100 >> $OUTFASTA
done

bedtools getfasta -fi $OUTFASTA -bed $VCF -fo $INTERSECTFASTA # we only get 26685 N
bedtools getfasta -fi $OUTFASTA -bed $VCF -fo $INTERSECTBED -tab