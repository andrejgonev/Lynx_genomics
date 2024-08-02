# Final pipeline for both species

# I installed the phyluce software in a conda environment (phyluce-1.7.3), as instructed in the manual.

compute -c 12 --mem 24

conda activate phyluce-1.7.3

cd /mnt/netapp2/Store_csebdjgl/agonev/data/

mkdir -p UCE/LYLY1_2A
cd UCE/
mkdir LYPA1_2A

# use fato2bit to convert the fasta genomes to 2bit format
# I downloaded the fato2bit tools in my softwares folder, following the instructions on https://hgdownload.soe.ucsc.edu/downloads.html#source_downloads

export PATH="/mnt/netapp2/Store_csebdjgl/agonev/software/twoBit:$PATH"

# Eurasian lynx
faToTwoBit /mnt/netapp2/Store_csebdjgl/agonev/software/snpEff/data/LYLY1_2A_anc/sequences.fa LYLY1_2A/LYLY1_2A.2bit
twoBitInfo LYLY1_2A/LYLY1_2A.2bit sizes.tab

#Iberian lynx
faToTwoBit /mnt/netapp2/Store_csebdjgl/lucia/polarization/ancestral_mLynPar1.2.scaffolds.revcomp.scaffolds.fa LYPA1_2A/LYPA1_2A.2bit
twoBitInfo LYPA1_2A/LYPA1_2A.2bit sizes.tab

# get the Tetrapods-UCE-5Kv1 probes from https://www.ultraconserved.org/
# I downloaded it locally, unziped it and moved it to the UCE folder with scp

# Align the probes to the genomes

phyluce_probe_run_multiple_lastzs_sqlite \
    --db lynx.sqlite \
    --output lynx-genome-lastz \
    --scaffoldlist LYLY1_2A LYPA1_2A \
    --genome-base-path ./ \
    --probefile Tetrapods-UCE-5Kv1.fasta \
    --cores 12


# create a genomes.conf file with the paths to the 2bit files. It should look like this:
[scaffolds]
LYLY1_2A:/mnt/netapp2/Store_csebdjgl/agonev/data/UCE/LYLY1_2A/LYLY1_2A.2bit
LYPA1_2A:/mnt/netapp2/Store_csebdjgl/agonev/data/UCE/LYPA1_2A/LYPA1_2A.2bit

# extract the sequences from the genomes

phyluce_probe_slice_sequence_from_genomes \
--lastz lynx-genome-lastz \
--conf genomes.conf \
--flank 500 \
--name-pattern "Tetrapods-UCE-5Kv1.fasta_v_{}.lastz.clean" \
--output lynx-genome-fasta

# Check if the sequence is correctly extracted from the genome fasta file
# fasta_file="your_fasta_file.fasta"
# start=179097896
# end=179098016

# awk -v start="$start" -v end="$end" 'BEGIN{sequence=""} /^>/{if(NR>1)exit} !/^>/{sequence=sequence$0} END{print substr(sequence, start, end-start+1)}' $fasta_file

# I checked the reference genome and the sequence is correctly extracted.


# Extract SNPs in the UCEs
vcf_file="/mnt/netapp2/Store_csebdjgl/lynx_genome/lynx_data/mLynLyn1.2_ref_vcfs/downsampling/polarization/c_ll_mLynLyn1.2_ref_downsampled.filter4.rd.miss30_polarized_annotated.vcf.gz"

# Create a bed file with the coordinates of the uces
# grep ">" lyly-genome-fasta/lyly1_2a.fasta | cut -d'|' -f2,5 | cut -d':' -f2,3 | sed 's/-/\t/g;s/|match:/\t/' > lyly-genome-fasta/lyly1_2a_uces.bed
# bedtools intersect -a $vcf_file -b lyly-genome-fasta/lyly1_2a_uces.bed > temp.vcf

# Try with bcftools instead
grep ">" lyly-genome-fasta/lyly1_2a.fasta | cut -d'|' -f2,5 | cut -d':' -f2,3 | sed 's/|match//' > lyly-genome-fasta/lyly1_2a_uces.txt
bcftools view -R lyly-genome-fasta/lyly1_2a_uces.txt $vcf_file > lyly_snps_uce.vcf # 191 snps
# this works so we keep this and delete the bed tools one (also the bed file). Because bed files are 0 based, we needed to convert it.


# # check if any of these 191 snps fall outside the gene boundaries (to get the snps in the UCNEs)
# gff_file="/mnt/netapp2/Store_csebdjgl/reference_genomes/lynx_lynx_mLynLyn1.2/LYLY1_2A.FA.gff3"
# awk '$3 == "gene"' $gff_file | awk -F"\t" '{print $1"\t"$4-1"\t"$5}' > lyly-genome-fasta/lyly1_2a_genes.bed

# bedtools intersect -a lyly_snps_uce.vcf -b lyly-genome-fasta/lyly1_2a_genes.bed -v > noncoding_snps.vcf #101 snps, which is 0.004% of the total snps (Dani had 0.02%)

# # because using bedtools intersect doesn't give the vcf header, we'll get it from the previous vcf and merge it with the noncoding snps vcf
# grep "^#" lyly_snps_uce.vcf > header.vcf
# cat header.vcf noncoding_snps.vcf > lyly_ucne_only.vcf
# rm header.vcf noncoding_snps.vcf


# # Get the counts of snps in the UCNEs
# ncols=$(grep -m1 '^#CHROM' lyly_ucne_only.vcf | tr '\t' '\n' | wc -l)

# for i in $(seq 10 $ncols); do
#     individual=$(grep -m1 '^#CHROM' lyly_ucne_only.vcf | cut -f"$i")
#     line="$individual"
#     echo "Processing individual $individual"
#     count=$(cut -f"$i" lyly_ucne_only.vcf | cut -d':' -f1 | tr '/' '+' | bc | awk '{ sum += $1 } END { print sum }')
#     line="$line\t$count"
#     echo -e "$line" >> counts_UCNE.tsv
# done


# Because UCNEs are often in the introns, we should intersect with the exons only to get the counts of snps in the UCNEs
cd counts

gff_file="/mnt/netapp2/Store_csebdjgl/reference_genomes/lynx_lynx_mLynLyn1.2/LYLY1_2A.FA.gff3"
awk '$3 == "exon"' $gff_file | awk -F"\t" '{print $1"\t"$4-1"\t"$5}' > lyly1_2a_exons.bed

bedtools intersect -a lyly_snps_uce.vcf -b lyly1_2a_exons.bed -v > noncoding_snps.vcf #176

# because using bedtools intersect doesn't give the vcf header, we'll get it from the previous vcf and merge it with the noncoding snps vcf
grep "^#" lyly_snps_uce.vcf > header.vcf
cat header.vcf noncoding_snps.vcf > lyly_ucne_no_exons.vcf
rm header.vcf noncoding_snps.vcf

# Get the counts of snps in the UCNEs
ncols=$(grep -m1 '^#CHROM' lyly_ucne_no_exons.vcf | tr '\t' '\n' | wc -l)

for i in $(seq 10 $ncols); do
    individual=$(grep -m1 '^#CHROM' lyly_ucne_no_exons.vcf | cut -f"$i")
    line="$individual"
    echo "Processing individual $individual"
    count=$(cut -f"$i" lyly_ucne_no_exons.vcf | cut -d':' -f1 | tr '/' '+' | bc | awk '{ sum += $1 } END { print sum }')
    line="$line\t$count"
    echo -e "$line" >> counts_UCNE_no_exons.tsv
done


### DO THE SAME FOR UCNEs from UCNEBase


cd /mnt/netapp2/Store_csebdjgl/agonev/data/UCE
mkdir -p UCNEbase

cd UCNEbase/

# get the human ucne fasta from ucnebase and treat it as our probe
wget https://epd.expasy.org/ucnebase/data/download/fasta/hg19_UCNEs.fasta.gz

# Align the probes to the genomes

phyluce_probe_run_multiple_lastzs_sqlite \
    --db lynx.sqlite \
    --output lynx-genome-lastz \
    --scaffoldlist LYLY1_2A LYPA1_2A \
    --genome-base-path ../ \
    --probefile hg19_UCNEs_probes.fasta \
    --cores 12

# create a genomes.conf file with the paths to the 2bit files. It should look like this:
[scaffolds]
LYLY1_2A:/mnt/netapp2/Store_csebdjgl/agonev/data/UCE/LYLY1_2A/LYLY1_2A.2bit
LYPA1_2A:/mnt/netapp2/Store_csebdjgl/agonev/data/UCE/LYPA1_2A/LYPA1_2A.2bit

# extract the sequences from the genomes

phyluce_probe_slice_sequence_from_genomes \
--lastz lynx-genome-lastz \
--conf genomes.conf \
--flank 500 \
--name-pattern "hg19_UCNEs_probes.fasta_{}.lastz.clean" \
--output lynx-genome-fasta

# This part did not work. The fasta file is not designed to be used as a probe. Some things are missing in the header. 

# Count the snps in the UCNEs from UCNEbase
vcf_file="/mnt/netapp2/Store_csebdjgl/lynx_genome/lynx_data/mLynLyn1.2_ref_vcfs/downsampling/polarization/c_ll_mLynLyn1.2_ref_downsampled.filter4.rd.miss30_polarized_annotated.vcf.gz"
mkdir counts
cd counts

awk -F"\t" '{print $2"\t"$4-1"\t"$5}' ../lynx-genome-lastz/hg19_UCNEs_probes.fasta_LYLY1_2A.lastz.clean > ucne_coords.bed

bedtools intersect -a $vcf_file -b ucne_coords.bed > snp_ucnebase_noheader.vcf

zgrep "^#" $vcf_file > header.vcf
cat header.vcf snp_ucnebase_noheader.vcf > snp_ucnebase.vcf
rm header.vcf snp_ucnebase_noheader.vcf

# Get the counts of snps in the UCNEs
ncols=$(grep -m1 '^#CHROM' snp_ucnebase.vcf | tr '\t' '\n' | wc -l)

for i in $(seq 10 $ncols); do
    individual=$(grep -m1 '^#CHROM' snp_ucnebase.vcf | cut -f"$i")
    line="$individual"
    echo "Processing individual $individual"
    count=$(cut -f"$i" snp_ucnebase.vcf | cut -d':' -f1 | tr '/' '+' | bc | awk '{ sum += $1 } END { print sum }')
    line="$line\t$count"
    echo -e "$line" >> counts_UCNEbase.tsv
done

### Finally I will do a check to see if the elements in the UCNEbase overlap with the UCNEs from using the Tetrapods probes
cd /mnt/netapp2/Store_csebdjgl/agonev/data/UCE/

# create a csv file with the coordinates of the UCNEs from UCNEbase
cut -f4,5 UCNEbase/lynx-genome-lastz/hg19_UCNEs_probes.fasta_LYLY1_2A.lastz.clean > file1.csv
# create a csv file with the coordinates of the UCNEs from the Tetrapods probes
cut -f2,3 lyly-genome-lastz/Tetrapods-UCE-5Kv1.fasta_LYLY1_2A.lastz.clean > file2.csv

# I will check the overlap the check_overlap_UCNEs.sh script

# Number of overlapping segments: 2020
# Number of exclusive segments in file1: 2332
# Number of exclusive segments in file2: 1445