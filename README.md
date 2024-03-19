# Balkan lynx genomic study

This is 

## Preprocessing

All FASTQ genome files used in the study need to be preprocessed first. This encompasses quality control, adapter sequence detection and trimming, poly-G tail trimming, base correction in overlapped regions etc. For this, I used fastp (Chen et al. 2018). Before running the software, I prepared a barcodes file (barcodes_filtered) which contains the sample codes of all samples that I will use, the directory paths to each of their reads, as well as the file names of the r1 and r2 (paired-end sequencing). Then I ran a simple script, adapted from Enrico, which is torun_fastp. 
For each line of the barcodes file, this small script will submit a job to the cluster using the Enrico's sbatch_fastp_fqdir_r1fq_r2fq.sh script. The settings for the fastp can be seen from the script. 

Some genomes were preprocessed before by Enrico.

Once this is done, we can move on to the alignment and bam generation.



## Alignments

The samples will be aligned to the new Eurasian lynx reference genome (link). 

To run the alignments and sorting, each individual needs to have its configuration file in yaml. This file contains information such as: sample name, all associated fastq ids (here I use the preprocessed fastq files with fastp), with the corresponding r1s and r2s, the path to them, the reference genome and its path, the output directory for the bams, and the softwared used for each of the steps. (explain them)

Some yaml configuration files were already done by Enrico for another study, so I used those, but changed all the paths and also the reference genomes, as we will now use the new Eurasian lynx reference genome. To generate the remaining yaml files, I prepared a table with the necessary information (input_data.txt) and used a bash script (link it) to extract information from it and output it in a separate yaml file for each sample. 

