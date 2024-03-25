#!/bin/bash
#SBATCH --time=3:00:00
#SBATCH --mem=12G
#SBATCH --cpus-per-task=6

# Remember to prepare the reference genome in advance! Like this:
# ${bwa} index ${reference_genome}
# ${samtools} faidx ${reference_genome}
# ${samtools} dict ${reference_genome} -o ${reference_genome/.fa/.dict}

# prepare the environment in CESGA ft3
module load cesga/2020 gcccore/system 
module load bwa/0.7.17
module load samtools/1.9
module load picard/2.25.5
module load gatk/3.7-0-gcfedb67

bwa index /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/reference_genomes/lynx_lynx_mLynLyn1.2/mLynLyn1.2.revcomp.scaffolds.fa
samtools faidx /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/reference_genomes/lynx_lynx_mLynLyn1.2/mLynLyn1.2.revcomp.scaffolds.fa
samtools dict /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/reference_genomes/lynx_lynx_mLynLyn1.2/mLynLyn1.2.revcomp.scaffolds.fa -o /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/reference_genomes/lynx_lynx_mLynLyn1.2/mLynLyn1.2.revcomp.scaffolds.dict

