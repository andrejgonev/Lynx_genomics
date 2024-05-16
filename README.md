# Balkan lynx genomic study

This is 

## Preprocessing

To ensure higher quality of the sequences, all FASTQ genome files used in the study need to be preprocessed first. This encompasses quality control, adapter sequence detection and trimming, poly-G tail trimming, base correction in overlapped regions etc. For this, I used [fastp](https://github.com/OpenGene/fastp) ([Chen et al. 2018](https://doi.org/10.1093/bioinformatics/bty560)). Before running the software, I prepared a barcodes file [barcodes_filtered](data/barcodes/barcodes_filtered) which contains the sample codes of all samples that I will use, the directory paths to each of their reads, as well as the file names of the R1 and R2 (paired-end sequencing). Then I ran a simple script, adapted from Enrico, which is run_fastp:

```
for row in $(cat data/barcodes/barcodes_filtered | tr ' ' ':'); do
    
    fastq_dir=$(echo $row | cut -d':' -f2)
    r1_fastq=$(echo $row | cut -d':' -f3)
    r2_fastq=$(echo $row | cut -d':' -f4)
    fastqid=$(echo $r1_fastq | sed 's/_1.fastq.gz$//')
    echo "sbatch scripts/sbatch_fastp_fqdir_r1fq_r2fq.sh ${fastq_dir} ${r1_fastq} ${r2_fastq}"
    sbatch -o logs/fastp/${fastqid}.out -e logs/fastp/${fastqid}.err scripts/sbatch_fastp_fqdir_r1fq_r2fq.sh ${fastq_dir} ${r1_fastq} ${r2_fastq}

done
```

For each line of the barcodes file, this small script will submit a job to the cluster using the Enrico's [sbatch_fastp_fqdir_r1fq_r2fq.sh](scripts/sbatch_fastp_fqdir_r1fq_r2fq.sh) script. Fastp was ran with the following flags:

```
--dont_overwrite (protect the existing files not to be overwritten by fastp)
--trim_poly_g (detect the polyG in read tails and trim them)
--length_required 30 (reads shorter than length_required will be discarded)
--correction (enable base correction in overlapped regions)
--detect_adapter_for_pe (enable adapter sequence auto-detection)
--thread 6 (worker thread number)
```


Some genomes were preprocessed before by Enrico and were not subjected to this step.

Once this is done, we can move on to the alignment and bam generation.



## Alignments

The samples will be aligned to the new high-quality Eurasian lynx reference genome (link). 

To run the alignments and sorting, each individual needs to have its configuration file in yaml. This file contains information such as: sample name, all associated fastq IDs (here I use the preprocessed fastq files with fastp), with the corresponding R1s and R2s, the path to them, the reference genome and its path, the output directory for the bams, and the software used for each of the steps. (explain them)

Some yaml configuration files were already done by Enrico for another [study](https://github.com/Enricobazzi/Lynxtrogression_v2/tree/main/config/alignment), so I used those, but changed all the paths and also the reference genomes, as we will now use the new Eurasian lynx reference genome. To generate the remaining yaml files, I prepared a table with the necessary information ([input_data.txt](config/test/input_data.txt)) and used a bash [script](config/test/group_individual_yamls.sh) to extract information from it and output it in a separate yaml file for each sample. 

Before running the alignment, I created a sequence dictionary of the reference genome, using [samtools](https://github.com/samtools/samtools).

```
bwa index mLynLyn1.2.revcomp.scaffolds.fa  # This was done before!
samtools faidx mLynLyn1.2.revcomp.scaffolds.fa
samtools dict mLynLyn1.2.revcomp.scaffolds.fa -o mLynLyn1.2.revcomp.scaffolds.dict
```

To run the script [sbatch_alignment_of_sample_from_yaml](scripts/sbatch_alignment_of_sample_from_yaml.sh), I did:

```
for yaml in $(ls config/alignment/); do

    sampleid=$(echo $yaml | sed 's/.mLynLyn.alignment.yaml$//')
    echo $sampleid
    sbatch -o logs/alignment/${sampleid}.out -e logs/alignment/${sampleid}.err scripts/sbatch_alignment_of_sample_from_yaml.sh config/alignment/$yaml

done
```

*** modify the next paragraph a bit, this is copied from enrico's git

The script makes use of [bwa v](), [samtools v](), [picard v]() and [gatk v]() and does the following: 
 - Align each of the sample's R1-R2 fastq pairs to the reference genome using BWA-MEM and Samtools view
 - samtools sort to sort the reads in the bam file
 - picardtools AddOrReplaceReadGroups to add read groups to the reads in the bam file
 - samtools merge to merge the bam files from the multiple R1-R2 pairs of the sample if necessary
 - picardtools MarkDuplicates to mark duplicate reads in the bam
 - realign indels with GATK which comprises:
 - gatk RealignerTargetCreator to identify realign targets
 - gatk IndelRealigner to realign the indels
 - samtools index to index the indel realigned bam for downstream analyses

 ## Alignment Quality Control

 Finally, before moving to the next steps, we need to assess the quality of the alignments. For this purpose, I used [QualiMap v]() (ref) with the [sbatch_qualimap_bam_out](scripts/sbatch_qualimap_bam_out.sh). The script was called with [run_qualimap](scripts/commands_and_small-scripts/run_qualimap.sh):

 ```
 # bam folder
bam_folder=/mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/mLynLyn1.2_ref_bams

# sbatch QualiMap of all samples
bams=$(ls ${bam_folder}/*er.bam)

for bam in $bams ; do

    sampleid=$(basename "$bam" | sed 's/_mLynLyn_ref.sorted.rg.merged_sorted.rmdup.indelrealigner.bam$//')
    echo "sbatch qualimap of ${sampleid}"
    sbatch -o logs/alignment/qualimap/qualimap-${sampleid}.out -e logs/alignment/qualimap/qualimap-${sampleid}.err \
    scripts/sbatch_qualimap_bam_out.sh ${bam} data/qualimap

done
```

To have a more global overview of all the QC results, we can summarise all the reports into one by using [multiqc]() (ref)

```
module load multiqc 

cd data/qualimap
multiqc .
```
A useful tutorial for navigating through multiqc reports can be found [here](https://youtu.be/qPbIlO_KWN0) 

This will help you identify any poor quality or weird samples, that should be discarded in the downstream analyses. I discarded two samples (c_ll_ca_0247 and c_ll_ca).

# Variant calling and filtering

After our alignments are done and we are satisfied with their quality, we can proceed to call the variants and build the VCFs.

## Variant calling

[table]
Number of SNPs after GLNexus merging: 10977787 (this applies a soft Qual >10 filter)
After first filter (low complexity and repeats) : 5472753
After second filter (non-biallelic sites and INDELs): 4238598
After third filter (invariant sites): 4237961
After fourth filter (QUAL>=30): 3757500




## Depth filtering 

I created a bam list like this: 

```
ls path/to/bams/mLynLyn1.2_ref_bams/*er.bam | grep -vE "ca_0249|ca_0253" > path/to/folder/agonev/data/c_ll_105.bamlist 
```

script [sbatch_mosdepth_10k_bam_outdir](scripts/sbatch_mosdepth_10k_bam_outdir.sh). I ran it like : 

```
inbams=($(cat data/c_ll_105.bamlist))

for bam in ${inbams[*]}; do
    echo "calculating depth of $bam"
    sbatch scripts/sbatch_mosdepth_10k_bam_outdir.sh \
        ${bam} data/variant_filtering/depth
done
```


[1] "bal fail: 663"
[1] "cau fail: 671"
[1] "crp fail: 592"
[1] "mng fail: 560"
[1] "wru fail: 561"
[1] "lva fail: 653"
[1] "nor fail: 635"
[1] "pol fail: 619"
[1] "tva fail: 579"
[1] "pyk fail: 605"
[1] "yak fail: 626"

[1] "all fail: 442"
[1] "bal fail and others pass: 0"


