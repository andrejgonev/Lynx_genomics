# Balkan lynx genomic study

This is 

# 1. Preprocessing, QC and Alignment

## 1.1. Preprocessing

To ensure higher quality of the sequences, all FASTQ genome files used in the study need to be preprocessed first. This encompasses quality control, adapter sequence detection and trimming, poly-G tail trimming, base correction in overlapped regions etc. For this, I used [fastp](https://github.com/OpenGene/fastp) ([Chen et al. 2018](https://doi.org/10.1093/bioinformatics/bty560)). Before running the software, I prepared a barcodes file [barcodes_filtered](data/barcodes/barcodes_filtered) which contains the sample codes of all samples that I will use, the directory paths to each of their reads, as well as the file names of the R1 and R2 (paired-end sequencing). Then I ran a simple script (run_fastp), adapted from Enrico:

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



## 1.2. Alignments

The samples will be aligned to the new high-quality Eurasian lynx reference genome (link). 

To run the alignments and sorting, each individual needs to have its configuration file in yaml. This file contains information such as: sample name, all associated fastq IDs (here I use the preprocessed fastq files with fastp), with the corresponding R1s and R2s, the path to them, the reference genome and its path, the output directory for the bams, and the software used for each of the steps. (explain them)

Some yaml configuration files were already done by Enrico for another [study](https://github.com/Enricobazzi/Lynxtrogression_v2/tree/main/config/alignment), so I used those, but changed all the paths and also the reference genomes, as we will now use the new Eurasian lynx reference genome. To generate the remaining yaml files, I prepared a table with the necessary information ([input_data.txt](config/test/input_data.txt)) and used a bash [script](config/test/group_individual_yamls.sh) to extract information from it and output it in a separate yaml file for each sample. 

Before running the alignment, I indexed and created a sequence dictionary of the reference genome, using [bwa](https://github.com/lh3/bwa) and [samtools](https://github.com/samtools/samtools).

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
 - Align each of the sample's R1-R2 fastq pairs to the reference genome using *BWA-MEM* and *Samtools view*
 - *Samtools sort* to sort the reads in the bam file
 - *picardtools AddOrReplaceReadGroups* to add read groups to the reads in the bam file
 - *Samtools merge* to merge the bam files from the multiple R1-R2 pairs of the sample if necessary
 - *picardtools MarkDuplicates* to mark duplicate reads in the bam
 - realign indels with GATK which comprises:
     - *gatk RealignerTargetCreator* to identify realign targets
    - *gatk IndelRealigner* to realign the indels
 - *Samtools index* to index the indel realigned bam for downstream analyses

 ### Alignment Quality Control

 Finally, before moving to the next steps, we need to assess the quality of the alignments. For this purpose, I used [QualiMap v2.2.1](http://qualimap.conesalab.org/) ([Okonechnikov et al. 2016](https://doi.org/10.1093/bioinformatics/btv566)) with the [sbatch_qualimap_bam_out](scripts/sbatch_qualimap_bam_out.sh). The script was called with [run_qualimap](scripts/commands_and_small-scripts/run_qualimap.sh):

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

To have a more global overview of all the QC results, we can summarise all the reports into one by using [multiqc](https://multiqc.info/) ([Ewels et al. (2016)](https://academic.oup.com/bioinformatics/article/32/19/3047/2196507)).

```
module load multiqc 

cd data/qualimap
multiqc .
```

A useful tutorial for navigating through multiqc reports can be found [here](https://youtu.be/qPbIlO_KWN0) 

This will help you identify any poor quality or weird samples, that should be discarded in the downstream analyses. I discarded two samples (c_ll_ca_0247 and c_ll_ca_0253).

# 2. Variant calling and filtering

After our alignments are done and we are satisfied with their quality, we can proceed to call the variants and build the VCFs.

## 2.1. Variant calling

For the variant calling, I used Google's [DeepVariant v1.6.0](https://github.com/google/deepvariant) deep learning caller with their WGS model.

...

For the pseudoautosomal regions, I followed [Lucia](https://github.com/luciamayorf/Variant_calling_and_filtering)'s method and established a standard PAR1 region spanning 7 Mb at the beginning and end of the X chromosome. She chose this length based on literature on both the Iberian and Eurasian lynx: [6 Mb](https://doi.org/10.1186/s13059-016-1090-1)-[6.5Mb](https://doi.org/10.1186%2Fs12864-017-3946-5) and [10 Mb](https://doi.org/10.1111/eva.13302).

To do this, I created a bed file ("mLynLyn1.2.PAR1_sexChr.bed"), where I put the length of the PARs (7 Mb) from both ends of the chromosome (0-7000000; 117087176-124087175), using the indexed fasta file to get the chromosome size. These regions will be treated as autosomal in all my male samples when doing the calling. 

I performed the calling with the script [variant_calling_deepvariant.sh](scripts/variant_calling_deepvariant.sh), which i ran for my samples like this: 

```
# list of bams to process
bams=$(ls /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/mLynLyn1.2_ref_bams/*er.bam | grep -vE "ca_0249|ca_0253")

# reference genome directory
ref_genome=/mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/reference_genomes/lynx_lynx_mLynLyn1.2/mLynLyn1.2.revcomp.scaffolds.fa

# sample list directory
sample_list=/mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/agonev/barcodes/samples_sex.txt

# output directory
output_dir=/mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/mLynLyn1.2_ref_gvcfs/

# loop through bams and submit jobs
for bam in $bams; do
    sample_name=$(basename ${bam} _mLynLyn_ref.sorted.rg.merged_sorted.rmdup.indelrealigner.bam)
    echo "sample_name: ${sample_name}"
    sbatch -o /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/agonev/logs/variant_calling/deepvariant/${sample_name}.out -e /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/agonev/logs/variant_calling/deepvariant/${sample_name}.err \
    -c 32 --mem=64GB -t 06:00:00 --gres=gpu:a100:1 /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/agonev/scripts/variant_calling_deepvariant.sh ${bam} ${ref_genome} ${sample_list} ${output_dir} ${sample_name}
done
```

#### Variant calling QC

As per Lucia's method (paranoia), I also did a QC of the VC to check if everything went well. I did the same as her:

```
module load samtools

for i in /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/mLynLyn1.2_ref_gvcfs/*_mLynLyn1.2_ref.vcf.gz; do
  bcftools stats ${i} > /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/mLynLyn1.2_ref_gvcfs/bcftools_stats/$(basename "${i}" .vcf.gz)_stats.txt
done

module load multiqc 
multiqc *_stats.txt
```


## 2.2. gVCF merging

The merging of the genome VCFs (gVCF) into a single combined VCF file was done with [GLnexus v1.4.1](https://github.com/dnanexus-rnd/GLnexus). This first creates a BCF file (a binary version), which then needs to be converted to a VCF. I ran GLnexus with the "DeepVariantWGS" configuration, which already applies some soft quality filters (AQ > 10). 

To do so, first I created a list of all the gVCFs that will be merged:

```
ls /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/mLynLyn1.2_ref_gvcfs/*g.vcf.gz > /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/mLynLyn1.2_ref_gvcfs/gvcfs_list.txt
```

Then, I ran the script [glnexus_script.sh](scripts/glnexus_script.sh): 

```
sbatch -c 32 --mem=100G -t 03:00:00 /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/agonev/scripts/glnexus_script.sh /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/mLynLyn1.2_ref_gvcfs/gvcfs_list.txt /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/mLynLyn1.2_ref_vcfs/c_ll_105_mLynLyn1.2_ref.vcf.gz 
``` 

#### VCF QC

I didn't do this, put it only for future reference as [Lucia did it](https://github.com/luciamayorf/Variant_calling_and_filtering/blob/main/README.md#vcf-qc) with her merged vcf.


```
module load samtools
module load multiqc

bcftools index -t filename.vcf.gz
bcftools stats filename.vcf.gz > filename_ref_stats.txt

multiqc filename_ref_stats.txt
```



### Removing half-calls

*This is something that I discovered later on in the process, but it is good to address it now, to avoid any unnecessary headaches in the downstream analyses.* 

GLnexus merging can sometimes result in "half-called" loci. These are all the genotypes that appear as ./0 and ./1 (refer to this [blog post](https://github.com/dnanexus-rnd/GLnexus/wiki/Reading-GLnexus-pVCFs) for a more detailed explanation). Perhaps there is an option in GLnexus, but I got rid of these half-calls like this: 

```

bcftools +setGT input.vcf -- -t ./x -n . > output_temp.vcf     # Recode all half-calls into missing genotypes ./.
bcftools +fill-tags output_temp.vcf -Ov -o output.vcf -- -t AC,AF,AN    # Recalculate the AC, AF and AN fields of the VCF INFO
rm output_temp.vcf

```


## 2.3. Variant filtering

### 2.3.1. Identify repetitive and low complexity regions

The [CNAG](https://www.cnag.eu/) team had previously identified all repetitive regions as part of the reference genome assembly, using [RepeatModeler](https://www.repeatmasker.org/RepeatModeler/). These regions are found in the reference genome directory, in the file **Repeats.4jb.gff3**. 

As it was done for the Iberian lynx (and for Enrico's [Lynxtrogression](https://github.com/Enricobazzi/Lynxtrogression_v2/blob/main/variant_filtering.md#find-repetitive-and-low-complexity-regions) project), I masked the genome both with the intersperse repeats in the Repeats.4jb.gff3 and with the low complexity regions, which I calculated using [RepeatMasker](https://www.repeatmasker.org/RepeatMasker/). To do this, I sbatched the script [repeatmasker_lowcomplex.sh](scripts/repeatmasker_lowcomplex.sh):

```
sbatch /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/agonev/scripts/repeatmasker_lowcomplex.sh
```

To combine the repeats with the low complexity regions, I did: 

```
ref_dir=/mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/reference_genomes/lynx_lynx_mLynLyn1.2

cat <(grep -v "#" ${ref_dir}/Repeats.4jb.gff3 | awk -F'\t' '{OFS="\t"; print $1, $4-1, $5}') \
    <(grep -v "#" ${ref_dir}/repetitive_regions/low_complex/mLynLyn1.2.revcomp.scaffolds.fa.out.gff | awk -F'\t' '{OFS="\t"; print $1, $4-1, $5}') |
    sort -k1,1 -k2,2n -k3,3n |
    bedtools merge -i - \
    > ${ref_dir}/repeats_lowcomplexity_regions.bed

```

The resulting bed file **repeats_lowcomplexity_regions.bed** will be used for filtering.

Finally, to see what portion of the genome I have masked, I checked the 

```
# To get the total length of the genome: 
awk '{sum+=$2} END {print sum}' ${ref_dir}/mLynLyn1.2.revcomp.scaffolds.fa.fai # 2432111198

# To calculate the length of these regions I run: 
awk '{sum+=$3-$2} END {print sum}' ${ref_dir}/repeats_lowcomplexity_regions.bed # 1047269226 

# To get the percentage of the genome that is masked:
echo "scale=2; 1047269226 / 2432111198 * 100" | bc  # 43.06%

```

With this step, I masked 1,047,269,226 bp, which is 43% of the total genome (2,432,111,198).

### 2.3.2. First round of filtering

The first round of filtering was done in a 4-step procedure: 
- Filter 1: removing variants in repetitive/low complexity regions.
- Filter 2: removing non-biallic sites and indels.
- Filter 3: removing invariant sites (substitutions from the reference genome, AF=1).
- Filter 4: removing variants with a low quality score (QUAL >= 30).

To apply these filters, I used the script [variant_filter_1to4](scripts/variant_filter_1to4.sh), which makes use of [bedtools](https://bedtools.readthedocs.io/en/latest/), [gatk](https://gatk.broadinstitute.org/hc/en-us) and [bcftools](https://samtools.github.io/bcftools/bcftools.html). 

To run this script, I did: 

```
vcf_dir=/mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/mLynLyn1.2_ref_vcfs
invcf=${vcf_dir}/c_ll_105_mLynLyn1.2_ref.vcf.gz
ref_dir=/mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/reference_genomes/lynx_lynx_mLynLyn1.2
ref=${ref_dir}/mLynLyn1.2.revcomp.scaffolds.fa
mask=${ref_dir}/repeats_lowcomplexity_regions.bed

sbatch --mem=12GB -t 03:00:00 scripts/variant_filter_1to4.sh ${ref} ${invcf} ${mask}
```

Summary of the first round of filtering:


| Filter                               | N of variants  | Filtered variants |
|:------------------------------------:|:--------------:|:-----------------:|
| 0. GLNexus merging (Qual >10)        |     10,977,787 |   0               |
| 1. Low complexity and repeats        |     5,472,753  |   5,505,034       |
| 2. Non-biallelic sites and INDELs    |     4,238,598  |   1,234,155       |
| 3. Invariant sites                   |     4,237,961  |   637             |
| 4. Quality filter (QUAL>=30)         |     3,757,500  |   480,461         |



### 2.3.3. Second round of filtering

The next two filters are applied on population, so to do that, I divided the VCF (filter4) into eleven population VCFs: 

- bal (Balkan Lynx)
- cau (Caucasian Lynx)
- crp (Carpathian Lynx)
- lva (Lynx in Latvia)
- nor (Lynx in Norway)
- pol (Lynx in NE Poland)
- tva (Lynx in Tuva region, Russia)
- mng (Lynx in Mongolia)
- wru (Lynx in Western Russia: Kirov and Ural regions)
- yak (Lynx in Yakutia, Russia)
- pyk (Lynx in Primorsky Krai, Russia)

To do this split, I used the [split_vcfs](scripts/commands_and_small-scripts/split_vcfs.sh) script.

### Depth filtering 

To avoid including possible paralogs in the analysis, one can eliminate genomic regions where an excess of sequencing reads align to the reference genome. Same as [Enrico](https://github.com/Enricobazzi/Lynxtrogression_v2/blob/main/variant_filtering.md#calculate-read-depth-filters-in-10k-bp-window), I decided to use 10kbp windows, where the mean read depth was calculated for each sample's BAM file using [mosdepth v0.3.2](https://github.com/brentp/mosdepth). 

I first created a list of bams like this: 
```
ls path/to/bams/mLynLyn1.2_ref_bams/*er.bam | grep -vE "ca_0249|ca_0253" > path/to/folder/agonev/data/c_ll_105.bamlist 
```

Then, I ran the script [sbatch_mosdepth_10k_bam_outdir](scripts/sbatch_mosdepth_10k_bam_outdir.sh): 

```
inbams=($(cat data/c_ll_105.bamlist))

for bam in ${inbams[*]}; do
    echo "calculating depth of $bam"
    sbatch scripts/sbatch_mosdepth_10k_bam_outdir.sh \
        ${bam} data/variant_filtering/depth
done
```

This then generates a bed file which I can use with the script [make_rdfilter_beds](data/make_rdfilter_beds.R) R script (a modified version of Enrico's [python](https://github.com/Enricobazzi/Lynxtrogression_v2/blob/main/src/variant_filtering/make_rdfilter_beds.py) script), to get the windows for populations whose sum of read depth values exceeds 1.5 times the mode of the values. These windows are then marked as "failed" for that population, and printed in a [table](). The script also outputs one bed file per population which contain the coordinates for the windows that do not pass the filter, draws plots and prints a summary:

```
bal fail: 663
cau fail: 671
crp fail: 592
mng fail: 560
wru fail: 561
lva fail: 653
nor fail: 635
pol fail: 619
tva fail: 579
pyk fail: 605
yak fail: 626

all fail: 442
bal fail and others pass: 0
```

<table>
    <tr>
        <td><img src="data/variant_filtering/depth/bal_depth_distribution.png" alt="Balkan" style="width: 90%;" /></td>
        <td><img src="data/variant_filtering/depth/cau_depth_distribution.png" alt="Caucasian" style="width: 90%;" /></td>
        <td><img src="data/variant_filtering/depth/crp_depth_distribution.png" alt="Carpathian" style="width: 90%;" /></td>
        <td><img src="data/variant_filtering/depth/mng_depth_distribution.png" alt="Mongolia" style="width: 90%;" /></td>
        <td><img src="data/variant_filtering/depth/wru_depth_distribution.png" alt="Western Russia" style="width: 90%;" /></td>
        <td><img src="data/variant_filtering/depth/lva_depth_distribution.png" alt="Latvia" style="width: 90%;" /></td>
    </tr>
    <tr>
        <td><img src="data/variant_filtering/depth/nor_depth_distribution.png" alt="Norway" style="width: 90%;" /></td>
        <td><img src="data/variant_filtering/depth/pol_depth_distribution.png" alt="Poland" style="width: 90%;" /></td> 
        <td><img src="data/variant_filtering/depth/tva_depth_distribution.png" alt="Tuva" style="width: 90%;" /></td>
        <td><img src="data/variant_filtering/depth/pyk_depth_distribution.png" alt="Primorski Krai" style="width: 90%;" /></td>
        <td><img src="data/variant_filtering/depth/yak_depth_distribution.png" alt="Yakutia" style="width: 90%;" /></td>       
    </tr>
</table>



### Missingness filtering

After many different methods tested, I decided to approach the missingness by filtering out the excessively missing variants in each population. This filter, like the read depth filtering, will be calculated on each population separately (this time on their VCFs) by using [bcftools filter](https://samtools.github.io/bcftools/howtos/filtering.html).

For this, I sbatch the script [missingness_bed_table](scripts/missingness_bed_table.sh), which generates a separate BED file for each of the defined threshold values for maximum missingness allowed (from 10 to 30 percent in 5% increments). This BED file contains the SNPs that would be excluded.

Besides, it also generates a [table](data/variant_filtering/missingness/miss_table.txt) which is useful for plotting and visualising the missingness. I plotted using the [miss_plot](data/variant_filtering/missingness/plot_miss.R) R script.

![Missingness_plot](data/variant_filtering/missingness/missing_plot.png). 

Finally, after visually inspecting the missingness, I'd like to apply the filter that would retain 95% of the SNPs in the original VCF (25% if I use all populations, tbd).

*This step is not necessary, but it's useful if one decides to change the missingness threshold at any point and explore more options (also visualize the data). If you know the missingness % you want to filter out, then it is possible to just use the `bcftools filter` command and pipe it in a new vcf.* 