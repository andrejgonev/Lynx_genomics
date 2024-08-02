## Synteny

## Fasta and BED generation

# Define general paths
REF=/GRUPOS/grupolince/reference_genomes/lynx_lynx_mLynLyn1.2/mLynLyn1.2.revcomp.scaffolds.fa
CHR=($(cat /GRUPOS/grupolince/reference_genomes/lynx_lynx_mLynLyn1.2/mLynLyn1.2.revcomp.scaffolds.fa.fai | cut -f 1 | uniq))
VCF=/GRUPOS/grupolince/mLynLyn1.2_ref_vcfs/c_ll_mLynLyn1.2_ref_downsampled.filter4.rd.miss30.filtered.vcf

############################ RUFUS ############################
# Define paths
INBAM=/GRUPOS/grupolince/mLynLyn1.2_ref_bams/polarization/lr1_mLynLyn1.2_ref_sorted_rg_merged_sorted_rmdup_indelrealigner.bam
OUTFASTA=/GRUPOS/grupolince/agonev/polarization/fasta/lr1_mLynLyn1.2_ref_sorted_rg_merged_sorted_rmdup_indelrealigner_pu2fa.fa
INTERSECTFASTA=/GRUPOS/grupolince/agonev/polarization/lr1_mLynLyn1.2_ref_intersect.fa
INTERSECTBED=/GRUPOS/grupolince/agonev/polarization/lr1_mLynLyn1.2_ref_intersect.bed

# get fasta from bam per chr
rm $OUTFASTA

for i in ${CHR[@]:0:20}; do
  echo "getting FASTA for chromosome $i"
  samtools mpileup -s -q30 -f $REF $INBAM -r $i | /GRUPOS/grupolince/Chrom-Compare/pu2fa -c $i -C 100 >> $OUTFASTA
done

echo "getting intersection FASTAs"

bedtools getfasta -fi $OUTFASTA -bed $VCF -fo $INTERSECTFASTA # we only get 23130 N
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

bedtools getfasta -fi $OUTFASTA -bed $VCF -fo $INTERSECTFASTA    # we only get 45060 N
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

bedtools getfasta -fi $OUTFASTA -bed $VCF -fo $INTERSECTFASTA # we only get 41900 N
bedtools getfasta -fi $OUTFASTA -bed $VCF -fo $INTERSECTBED -tab


# i removed all the scaffolds from the bedfiles with grep -v "scaffold"
# remember that i have the ChrY in the vcf

cut -f2 pm_mLynLyn1.2_ref_intersect.bed | paste fc_mLynLyn1.2_ref_intersect.bed - > fc_pm_sinteny_intersect.bed
cut -f2 lr1_mLynLyn1.2_ref_intersect.bed | paste fc_pm_sinteny_intersect.bed - > temp_fc_pm_lr.bed

zgrep -v "#" /GRUPOS/grupolince/mLynLyn1.2_ref_vcfs/c_ll_mLynLyn1.2_ref_downsampled.filter4.rd.miss30.filtered.vcf | grep -v "scaffold\|ChrY" | cut -f4 | paste -d'\0' temp_fc_pm_lr.bed - | awk -F"\t|:|-" '{printf ("%s\t%s\t%s\t%s=%s%s%s\n", $1,$2,$3,"fc_pm_lr_ll",$4,$5,$6)}' > fc_pm_lr_ll_sinteny_intersect.bed
rm temp_fc_pm_lr.bed

### parsimony 


awk '{                                       
split($0,a,":");
split(a[1],b,"=");
split(b[2],c,"");
if (c[1]==c[2] && c[1]==c[3]) printf ("%s\t%s\t%s\t%s\t%s\t%s\n", $1,$2,$3,$4,c[3],c[4]); # all equal
else if (c[1]==c[2] && c[1]!=c[3]) printf ("%s\t%s\t%s\t%s\t%s\t%s\n", $1,$2,$3,$4,"N",c[4]); #fc pm vs lr N
else if (c[1]==c[3] && c[1]!=c[2]) printf ("%s\t%s\t%s\t%s\t%s\t%s\n", $1,$2,$3,$4,c[3],c[4]); #fc lr* vs pm
else if (c[2]==c[3] && c[1]!=c[3]) printf ("%s\t%s\t%s\t%s\t%s\t%s\n", $1,$2,$3,$4,c[3],c[4]); #pm lr * vs fc

else if ((c[1]=="N") && c[2]==c[3]) printf ("%s\t%s\t%s\t%s\t%s\t%s\n", $1,$2,$3,$4,c[3],c[4]); # fc lr * all equal
else if ((c[1]=="N") && c[2]!=c[3]) printf ("%s\t%s\t%s\t%s\t%s\t%s\n", $1,$2,$3,$4,"N",c[4]); # lr vs fc N

else printf ("%s\t%s\t%s\t%s\t%s\t%s\n", $1,$2,$3,$4,"N",c[4]); #others null
}' fc_pm_lr_ll_sinteny_intersect.bed > outgroup_parsimony_ancestral_state_fc_pm_lr_ll_variants.bed


# get the unpolarizable 
awk '$5=="N" {print $0}' outgroup_parsimony_ancestral_state_fc_pm_lr_ll_variants.bed >  unpolarizable_fc_pm_lr_ll_variants.bed

# get the polarizable
awk '$5!="N" {print $0}' outgroup_parsimony_ancestral_state_fc_pm_lr_ll_variants.bed >  polarizable_fc_pm_lr_ll_variants.bed


# Generate file with inconsistent sites (wrongly polarised - ancestral state needs to be changed - or unpolarizable) --> for the ancestral reference fasta generation
awk '$5!=$6 {print $0}' outgroup_parsimony_ancestral_state_fc_pm_lr_ll_variants.bed > inconsistent_ancestral_state_fc_pm_lr_ll_variants.bed         
    


# Create the tsv file to annotate the VCF and index it
awk 'BEGIN {OFS="\t"} {print $1, $3, $5}' polarizable_fc_pm_lr_ll_variants.bed | bgzip -c > aa_annotation_polarizable_snps_fc_pm_lr_ll.tsv.gz
tabix -s1 -b2 -e2 aa_annotation_polarizable_snps_fc_pm_lr_ll.tsv.gz            # -b2 and -e2 to indicate that it is 1-based


## POLARIZING THE VCF

cd /GRUPOS/grupolince/mLynLyn1.2_ref_vcfs

# I generate a VCF file excluding the SNPs in scaffolds and ChrY_unloc (I should've done this from the beginning)
zgrep -v "_scaffold_\|ChrY" c_ll_mLynLyn1.2_ref_downsampled.filter4.rd.miss30.filtered.vcf > c_ll_mLynLyn1.2_ref_downsampled.filter4.rd.miss30_bigChr.vcf
bgzip c_ll_mLynLyn1.2_ref_downsampled.filter4.rd.miss30_bigChr.vcf

# First, we are going to keep only the sites in my VCF that are polarizable.
bedtools intersect -a c_ll_mLynLyn1.2_ref_downsampled.filter4.rd.miss30_bigChr.vcf.gz -b /GRUPOS/grupolince/agonev/polarization/polarizable_fc_pm_lr_ll_variants.bed -header > c_ll_mLynLyn1.2_ref_downsampled.filter4.rd.miss30_bigChr_polarizable.vcf
bgzip c_ll_mLynLyn1.2_ref_downsampled.filter4.rd.miss30_bigChr_polarizable.vcf

# Then we add the ancestral allele column to the INFO field
/opt/bcftools-1.6/bcftools annotate -a /GRUPOS/grupolince/agonev/polarization/aa_annotation_polarizable_snps_fc_pm_lr_ll.tsv.gz -c CHROM,POS,INFO/AA c_ll_mLynLyn1.2_ref_downsampled.filter4.rd.miss30_bigChr_polarizable.vcf.gz -h <(echo '##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral Allele">') -Oz -o c_ll_mLynLyn1.2_ref_downsampled.filter4.rd.miss30_bigChr_polarizable_aafilled.vcf.gz


## To polarize the SNPs in the VCF (alleles will be switched whenever the ancestral allele matches the alternative one, and genotypes will be properly recoded as well, based on the new INFO/AA column). I take the code from Dani's improved_polarization.Rmd, line 4372

# The following code was originally provided by Pierre Lindenbaum and modified by José Luis Castro.
java -jar /opt/jvarkit/dist/vcffilterjdk.jar -e 'if(variant.getNAlleles()!=2 || !variant.hasAttribute("AA")) return true; 
final String aa = variant.getAttributeAsString("AA",""); 
if(!variant.getAlleles().get(1).getDisplayString().equalsIgnoreCase(aa)) return true; 
VariantContextBuilder vb=new VariantContextBuilder(variant); 

Allele oldalt = variant.getAlleles().get(1);
Allele oldref = variant.getAlleles().get(0); 
Allele ref= Allele.create(oldalt.getDisplayString(),true); 
Allele alt= Allele.create(oldref.getDisplayString(),false);

vb.alleles(Arrays.asList(ref,alt)); 

List genotypes= new ArrayList<>(); 
for(Genotype g: variant.getGenotypes()) 
  { 
  if(!g.isCalled()) 
  { genotypes.add(g); continue;} 
  GenotypeBuilder gb = new GenotypeBuilder(g); 
  List alleles = new ArrayList<>(); 
  for(Allele a:g.getAlleles()) { 
    if(a.equals(oldalt)) { a=ref;} 
    else if(a.equals(oldref)) { a=alt;} 
    alleles.add(a); 
    } 
  if(g.hasPL()) { 
    int pl[] = g.getPL(); 
    int pl2[] = new int[pl.length]; 
    for(int i=0;i< pl.length;i++) pl2[i]=pl[(pl.length-1)-i]; 
    gb.PL(pl2); 
    } 
  if(g.hasAD()) 
    { int ad[] = g.getAD(); 
    int ad2[] = new int[ad.length]; 
    for(int i=0;i< ad.length;i++) ad2[i]=ad[(ad.length-1)-i];
    gb.AD(ad2); 
  } 
  genotypes.add(gb.alleles(alleles).make()); 
  }

vb.attribute("AF",1.0d - Double.parseDouble(variant.getAttributeAsString("AF",""))); vb.attribute("AC",variant.getGenotypes().stream().flatMap(G->G.getAlleles().stream()).filter(A->A.equals(oldref)).count()); 
vb.genotypes(genotypes); 
return vb.make();' -o c_ll_mLynLyn1.2_ref_downsampled.filter4.rd.miss30_bigChr_polarizable_aafilled_polarized.vcf.gz c_ll_mLynLyn1.2_ref_downsampled.filter4.rd.miss30_bigChr_polarizable_aafilled.vcf.gz

#incompatble 

cut -f2 lr1_mLynLyn1.2_ref_intersect.bed | paste fc_pm_sinteny_intersect.bed - > temp_fc_pm_lr.bed
zgrep -v "#" /GRUPOS/grupolince/mLynLyn1.2_ref_vcfs/c_ll_mLynLyn1.2_ref_downsampled.filter4.rd.miss30.filtered.vcf | grep -v "scaffold\|ChrY" | cut -f4,5 | sed 's/\t//' | paste -d'\0' temp_fc_pm_lr.bed - | awk -F"\t|:|-" '{printf ("%s\t%s\t%s\t%s=%s%s%s\n", $1,$2,$3,"fc_pm_lr_llRefAlt",$4,$5,$6)}' > temp_incompatibilities.bed

awk '{                                       
split($0,a,":");
split(a[1],b,"=");
split(b[2],c,"");
if (c[1]==c[2] && c[1]==c[3]) printf ("%s\t%s\t%s\t%s %s\n", $1,$2,$3,$4,c[3]); # all equal
else if (c[1]==c[2] && c[1]!=c[3]) printf ("%s\t%s\t%s\t%s %s\n", $1,$2,$3,$4,"N"); #fc pm vs lr N
else if (c[1]==c[3] && c[1]!=c[2]) printf ("%s\t%s\t%s\t%s %s\n", $1,$2,$3,$4,c[3]); #fc lr* vs pm
else if (c[2]==c[3] && c[1]!=c[3]) printf ("%s\t%s\t%s\t%s %s\n", $1,$2,$3,$4,c[3]); #pm lr * vs fc

else if ((c[1]=="N") && c[2]==c[3]) printf ("%s\t%s\t%s\t%s %s\n", $1,$2,$3,$4,c[3]); # fc lr * all equal
else if ((c[1]=="N") && c[2]!=c[3]) printf ("%s\t%s\t%s\t%s %s\n", $1,$2,$3,$4,"N"); # lr vs fc N

else printf ("%s\t%s\t%s\t%s %s\n", $1,$2,$3,$4,"N"); #others null
}' temp_incompatibilities.bed | awk -F ' ' '$5!="N" {print $0}' > temp_incompatibilities_polarizable.bed



# Compare both of the ll alleles with the other species alleles to see if there's any incompatibility
sed 's/ //g' temp_incompatibilities_polarizable.bed |  
awk '{                                       
split($0,a,":");
split(a[1],b,"=");
split(b[2],c,"");
if (c[4]==c[6]) printf ("%s\t%s\t%s\t%s\t%s\n", $1,$2,$3,$4,"ok_ref"); # ref allele matches any of the other species alleles
else if (c[5]==c[6]) printf ("%s\t%s\t%s\t%s\t%s\n", $1,$2,$3,$4,"ok_alt"); # alt allele matches any of the other species alleles
else printf ("%s\t%s\t%s\t%s\t%s\t%s\n", $1,$2,$3,$4,"incompatible","N"); # none of the alleles matches any of the other species alleles
}' | grep "incompatible" > allele_incompatibities_polarization.bed 

# 5653 incompatible sites

rm temp*.bed 

# To change the ancestral state for an N in that file:

awk 'NR==FNR{a[$1,$2,$3]=$0;next}{if(($1,$2,$3) in a){print $1"\t"$2"\t"$3"\t"$4"\t""N""\t"$6}else{print $0}}' allele_incompatibities_polarization.bed inconsistent_ancestral_state_fc_pm_lr_ll_variants.bed > ancestral_state_changes_fc_pm_lr_ll_variants.bed

#I also need to remove those variants from the polarized VCF

bedtools intersect -a /GRUPOS/grupolince/mLynLyn1.2_ref_vcfs/c_ll_mLynLyn1.2_ref_downsampled.filter4.rd.miss30_bigChr_polarizable_aafilled_polarized.vcf.gz -b allele_incompatibities_polarization.bed -v -header > /GRUPOS/grupolince/mLynLyn1.2_ref_vcfs/c_ll_mLynLyn1.2_ref_downsampled.filter4.rd.miss30_polarized.vcf
bgzip /GRUPOS/grupolince/mLynLyn1.2_ref_vcfs/c_ll_mLynLyn1.2_ref_downsampled.filter4.rd.miss30_polarized.vcf

# I have 2,358,779 SNPs in the final polarized VCF


##Ancestral reference fasta generation
# GATK FastaAlternateReferenceMaker takes a VCF and the reference genome as inputs and replaces the reference variant for the alternative one for the positions in the VCF.

# For that, I first need to generate a "pseudovcf" (it will only have useful information about the alleles and the SNP position) file from the BED file containing the positions that need to be changed (the alternative allele column containing the ancestral state).

# Generate the VCF
awk 'BEGIN {printf("##fileformat=VCFv4.2\n##INFO=<ID=END,Number=1,Type=Integer>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");} {printf("%s\t%s\t.\t%s\t%s\t.\t.\tEND=%d\n",$1,int($2)+1,$6,$5,$3);}' ancestral_state_changes_fc_pm_lr_ll_variants.bed > ancestral_state_changes_fc_pm_lr_ll_variants.vcf

# Generate the ancestral fasta
java -jar /opt/GATK/GenomeAnalysisTK.jar -T FastaAlternateReferenceMaker \
-R /GRUPOS/grupolince/reference_genomes/lynx_lynx_mLynLyn1.2/mLynLyn1.2.revcomp.scaffolds.fa \
-V ancestral_state_changes_fc_pm_lr_ll_variants.vcf \
-o ancestral_mLynLyn1.2.revcomp.scaffolds.fa