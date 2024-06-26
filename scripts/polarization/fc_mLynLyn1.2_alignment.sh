#!/bin/bash
# Define variable paths
FASTP_DIR=/GRUPOS/grupolince/rawdata/felis_catus_SRA/fastp
BAM_DIR=/GRUPOS/grupolince/mLynLyn1.2_ref_bams/polarization
REF=/GRUPOS/grupolince/reference_genomes/lynx_lynx_mLynLyn1.2/mLynLyn1.2.revcomp.scaffolds.fa

# align reads from files 'SRR5055389_1.fastq.gz' and 'SRR5055389_2.fastq.gz':
bwa mem $REF $FASTP_DIR/SRR5055389_1.fastp.fastq.gz $FASTP_DIR/SRR5055389_2.fastp.fastq.gz -t 20 | samtools view -hbS -@ 20 - -o $BAM_DIR/fc_mLynLyn1.2_ref.bam
echo "Reads aligned to the reference genome"

# sort the resulting bam file:
samtools sort -@ 20 $BAM_DIR/fc_mLynLyn1.2_ref.bam -o $BAM_DIR/fc_mLynLyn1.2_ref_sorted.bam
echo "Bam file sorted"

# add read groups to the sorted bam file:
java -jar /opt/picard-tools/picard.jar AddOrReplaceReadGroups I=$BAM_DIR/fc_mLynLyn1.2_ref_sorted.bam O=$BAM_DIR/fc_mLynLyn1.2_ref_sorted_rg.bam RGID=cinammon RGLB=SRR5055389fc_lib RGPL=Illumina RGPU=cinammon RGSM=SRR5055389fc VALIDATION_STRINGENCY=SILENT
echo "Read groups added"

# There's only one fastqid!
# renaming its only sorted_rg BAM to merged_sorted:
mv $BAM_DIR/fc_mLynLyn1.2_ref_sorted_rg.bam $BAM_DIR/fc_mLynLyn1.2_ref_sorted_rg_merged_sorted.bam

# marking duplicates of merged bam:
java -jar /opt/picard-tools/picard.jar MarkDuplicates METRICS_FILE=$BAM_DIR/fc_rmdup.txt I=$BAM_DIR/fc_mLynLyn1.2_ref_sorted_rg_merged_sorted.bam O=$BAM_DIR/fc_mLynLyn1.2_ref_sorted_rg_merged_sorted_rmdup.bam MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=800
echo "Duplicates marked"

# index the duplicate marked bam for gatk:
samtools index $BAM_DIR/fc_mLynLyn1.2_ref_sorted_rg_merged_sorted_rmdup.bam
echo "Bam file indexed"

# identify indel realignment targets:
java -jar /opt/GATK/GenomeAnalysisTK.jar -T RealignerTargetCreator -nt 20 -R $REF -I $BAM_DIR/fc_mLynLyn1.2_ref_sorted_rg_merged_sorted_rmdup.bam -o $BAM_DIR/fc_realignertargetcreator.intervals
echo "Indel realignment targets identified"

# realign indels:
java -jar /opt/GATK/GenomeAnalysisTK.jar -T IndelRealigner -R $REF -targetIntervals $BAM_DIR/fc_realignertargetcreator.intervals -I $BAM_DIR/fc_mLynLyn1.2_ref_sorted_rg_merged_sorted_rmdup.bam -o $BAM_DIR/fc_mLynLyn1.2_ref_sorted_rg_merged_sorted_rmdup_indelrealigner.bam
echo "Indels realigned"

# index the indel realigned bam for future analyses:
samtools index $BAM_DIR/fc_mLynLyn1.2_ref_sorted_rg_merged_sorted_rmdup_indelrealigner.bam
echo "Final bam file indexed"