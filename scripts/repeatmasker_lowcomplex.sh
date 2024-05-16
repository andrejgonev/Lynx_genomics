#!/bin/bash

#SBATCH --job-name=repeatmasker_ll
#SBATCH --output=/mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/agonev/logs/repetitive_regions/repeatmasker_lowcomplex.out
#SBATCH --error=/mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/agonev/logs/repetitive_regions/repeatmasker_lowcomplex.err
#SBATCH --time=1-00:00:00
#SBATCH --mem=40G
#SBATCH --cpus-per-task=36

module load repeatmasker

RepeatMasker -pa 9 -gff -xsmall -noint -dir /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/reference_genomes/lynx_lynx_mLynLyn1.2/repetitive_regions/low_complex/ /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/reference_genomes/lynx_lynx_mLynLyn1.2/mLynLyn1.2.revcomp.scaffolds.fa