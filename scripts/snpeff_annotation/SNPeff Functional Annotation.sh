## SNPeff Functional Annotation

module load snpeff/5.0
cp /opt/cesga/2020/software/Core/snpeff/5.0/snpEff.config /mnt/netapp2/Store_csebdjgl/agonev/software/snpEff


# add to config file
# Lynx_lynx, version mLynLyn1.2
LYLY1_2A.genome : Eurasian lynx

# Longest isoform selection
# To build the database, SnpEff needs all the files to be in the same folder (named after the reference genome code in the configuration file). It needs the following files: the reference genome FASTA, the GFF file, and a coding and protein sequences fasta).

# To get one effect per variant, we need to select one transcript and one protein per gene. We will base out decision on the longest peptide that it's transcribed (and its correspondent transcript).

# We have the FASTA with the protein sequence of the longest peptide available (provided by CNAG with the new reference genomes). We will start with this file to build the longest transcript file, after which we will generate the new gff3 that contains only the longest isoform.

cd /mnt/netapp2/Store_csebdjgl/reference_genomes/lynx_lynx_mLynLyn1.2
# For the transcripts file
for i in $(zgrep ">"  LYLY1_2A.longestpeptide.fa.gz | sed '/^>/ s/\(.*\)P/\1T/;s/>//'); do
    awk -v id="${i}" -v RS='>' -v ORS='' '$1 == id {print ">"$0}' LYLY1_2A.transcripts.fa
done > snpeff_annotation/LYLY1_2A.transcripts_longest.fa

# To check that they match:
diff <(zgrep ">" LYLY1_2A.longestpeptide.fa.gz | sed '/^>/ s/\(.*\)P/\1T/' | sed 's/>//' | sort) <(grep ">" snpeff_annotation/LYLY1_2A.transcripts_longest.fa | sed 's/>//' | sort)



## To generate the new GFF3 file, I use the custom script gff_longest_generation: (but be in the ref genome directory)
sbatch -t 04:00:00 --mem=500MB /mnt/netapp2/Store_csebdjgl/agonev/scripts/gff_longest_generation.sh 
# JOB ID: 7938316


# Now that we have all the files ready, we move them to the same folder and rename them so that SnpEff can build the database

# create a directory inside the software's dependencies whose name matches the code
mkdir /mnt/netapp2/Store_csebdjgl/agonev/software/snpEff/data/LYLY1_2A

# copy the reference genome FASTA to the new folder and rename it following the manuals instructions
cp /mnt/netapp2/Store_csebdjgl/reference_genomes/lynx_lynx_mLynLyn1.2/mLynLyn1.2.revcomp.scaffolds.fa /mnt/netapp2/Store_csebdjgl/agonev/software/snpEff/data/LYLY1_2A/
mv mLynLyn1.2.revcomp.scaffolds.fa sequences.fa

# copy and rename the gff3 file
cp /mnt/netapp2/Store_csebdjgl/reference_genomes/lynx_lynx_mLynLyn1.2/snpeff_annotation/LYLY1_2A.FA_longestisoform.gff3 /mnt/netapp2/Store_csebdjgl/agonev/software/snpEff/data/LYLY1_2A/
mv /mnt/netapp2/Store_csebdjgl/agonev/software/snpEff/data/LYLY1_2A/LYLY1_2A.FA_longestisoform.gff3 /mnt/netapp2/Store_csebdjgl/agonev/software/snpEff/data/LYLY1_2A/genes.gff

# copy and rename the transcripts file 
cp /mnt/netapp2/Store_csebdjgl/reference_genomes/lynx_lynx_mLynLyn1.2/snpeff_annotation/LYPA1_2A.transcripts_longest_CESGA.fa /mnt/netapp2/Store_csebdjgl/agonev/software/snpEff/data/LYLY1_2A/
mv /mnt/netapp2/Store_csebdjgl/reference_genomes/lynx_lynx_mLynLyn1.2/snpeff_annotation/LYPA1_2A.transcripts_longest_CESGA.fa /mnt/netapp2/Store_csebdjgl/agonev/software/snpEff/data/LYLY1_2A/cds.fa

# move and rename the protein file (the code for the protein needs to be the same as the code for the transcript, so we need to change the last "P" for a "T").
zcat /mnt/netapp2/Store_csebdjgl/reference_genomes/lynx_lynx_mLynLyn1.2/LYLY1_2A.longestpeptide.fa.gz | sed '/^>/ s/\(.*\)P/\1T/' > /mnt/netapp2/Store_csebdjgl/agonev/software/snpEff/data/LYLY1_2A/protein.fa



# build the database using the following command in a compute node
compute -c 6 --mem 24
module load snpeff/5.0

java -jar $EBROOTSNPEFF/snpEff.jar build -d -gff3 -v LYLY1_2A -c /mnt/netapp2/Store_csebdjgl/agonev/software/snpEff/snpEff.config -dataDir /mnt/netapp2/Store_csebdjgl/agonev/software/snpEff/data


## Running SNPeff

SNPeff_dir="/mnt/netapp2/Store_csebdjgl/agonev/software/snpEff"
vcf_dir="/mnt/netapp2/Store_csebdjgl/lynx_genome/lynx_data/mLynLyn1.2_ref_vcfs/downsampling"

java -Xmx16g -jar $EBROOTSNPEFF/snpEff.jar LYLY1_2A \
 -v -c $SNPeff_dir/snpEff.config \
 -s "$SNPeff_dir/data/LYLY1_2A/snpeff_c_ll_mLynLyn1.2_ref_downsampled.filter4.rd.miss30.html" \
 -csvStats "$SNPeff_dir/data/LYLY1_2A/snpeff_c_ll_mLynLyn1.2_ref_downsampled.filter4.rd.miss30.csv" \
 $vcf_dir/c_ll_mLynLyn1.2_ref_downsampled.filter4.rd.miss30.filtered.vcf > $vcf_dir/functional_annotation/snpeff/c_ll_mLynLyn1.2_ref_downsampled.filter4.rd.miss30_annotated.vcf.gz

 #####

 # convert 
 # copy the ancestral genome to the CESGA server
scp ancestral_mLynLyn1.2.revcomp.scaffolds.fa cseyeeba@ft3.cesga.es:/mnt/netapp2/Store_csebdjgl/agonev/data/polarization

# convert the ancestral genome to the same format as the reference genome
module load samtools/1.10
samtools faidx ancestral_mLynLyn1.2.revcomp.scaffolds.fa