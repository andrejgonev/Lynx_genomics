 ############################################################################################################

 # copy the ancestral genome to the CESGA server
scp ancestral_mLynLyn1.2.revcomp.scaffolds.fa cseyeeba@ft3.cesga.es:/mnt/netapp2/Store_csebdjgl/agonev/data/polarization

# copy the polarized vcf to the CESGA server
scp c_ll_mLynLyn1.2_ref_downsampled.filter4.rd.miss30_polarized.vcf.gz cseyeeba@ft3.cesga.es:/mnt/netapp2/Store_csebdjgl/lynx_genome/lynx_data/mLynLyn1.2_ref_vcfs/downsampling/polarization
gunzip c_ll_mLynLyn1.2_ref_downsampled.filter4.rd.miss30_polarized.vcf.gz

### 

# First we install bioperl 
cpan Module::Bio # Maybe this is not necessary
cpan App::cpanminus # Install the app to install Bio::Perl
cpanm Bio::Perl --force # Install Bio::Perl, force because it didn't work peacefully

# Source the .bashrc to load the new paths and varaibles
source ~/.bashrc

# We create the script necessary to run the CNAG script, SeqOp.pm
mkdir -p /mnt/netapp2/Store_csebdjgl/agonev/software/perlmods
cd /mnt/netapp2/Store_csebdjgl/agonev/software/perlmods
nano SeqOp.pm # paste the script from CNAG and change the use lib to the correct path containing bioperl: "/home/csic/eye/eba/perl5/lib/perl5/Bio";

# Then in the same directory, we created the script to extract the info from GFF and FASTA
nano CDS2seq.v2.pl
# we added the directory to the SeqOP.pm script in the use lib line: ""/mnt/netapp2/Store_csebdjgl/agonev/software/perlmods";"

#Finally we ran the script as follows:
mkdir Test_CNAG
perl /mnt/netapp2/Store_csebdjgl/agonev/scripts/CDS2seq.v2.pl -i ../LYLY1_2A.FA.gff3 -f ../mLynLyn1.2.revcomp.scaffolds.fa

#Check the difference between the CDS files

diff <(zcat LYLY1_2A.cds.fa.gz) Test_CNAG/LYLY1_2A.FA.cds.fa # no differences!!!!!

# Check the difference between the protein files
diff <(zcat LYLY1_2A.pep.fa.gz) Test_CNAG/LYLY1_2A.FA.pep.fa # The problem with one more amino acid persists. They have one more amino acid than us.

# Get the longest peptide from the newly created pep.fa file

# put the transcripts next to the protein sequence names, tab separated
awk '/^>/ {if (seq) print seq; printf "%s ", $0; seq=""; next} {seq=seq""$0} END {if (seq) print seq}' LYLY1_2A.FA.pep.fa > LYLY1_2A_tab.pep.fa

# To get the longest peptide for each gene, we need to get the longest peptide for each gene
for i in $(sed 's/\(>LYLY1_2A[0-9]*\)P[0-9]*/\1/' LYLY1_2A_tab.pep.fa | cut -f1 -d' ' | sort -u) ; do 
    grep "${i}" LYLY1_2A_tab.pep.fa | sed 's/partial_cds //' | awk '(length($2) > max_length) { max_length = length($2); max_line = $0 } END { print max_line }'  >> LYLY1_2A_tab.longest.pep.fa
done

# Check the difference between the longest peptide files
diff <(zcat LYLY1_2A.longestpeptide.fa.gz | awk '/^>/ {if (seq) print seq; printf "%s ", $0; seq=""; next} {seq=seq""$0} END {if (seq) print seq}' | sort -k1) <(sed 's/partial_cds //' Test_CNAG/LYLY1_2A_tab.longest.pep.fa | sort -k1) | less -S
# Only 75 differences, because we selected the longest peptide for each gene, but we have one amino acid less then the original file. Sometimes this is enough to cause a diferent peptide to be longer.

# Finally, reverse it back to a fasta format to translate
awk '{print $1"\n"$2}'  Test_CNAG/LYLY1_2A_tab.longest.pep.fa > Test_CNAG/LYLY1_2A.longest.pep.fa


# Now do the same with the ancestral reference genome: /mnt/netapp2/Store_csebdjgl/agonev/data/polarization/ancestral_mLynLyn1.2.revcomp.scaffolds.fa

# change the chromosome names to only match the original reference genome
sed -E '/^>/ s/^>([0-9]+ )?(mLynLyn[^:]*)(_rc)?(:1)?$/>\2/' ancestral_mLynLyn1.2.revcomp.scaffolds.fa > ancestral_mLynLyn1.2.revcomp.scaffolds_chr.fa

# Since it worked, I'll overwrite the original file
mv ancestral_mLynLyn1.2.revcomp.scaffolds_chr.fa ancestral_mLynLyn1.2.revcomp.scaffolds.fa

# Now we run the script to get the CDS and protein sequence fastas
perl /mnt/netapp2/Store_csebdjgl/agonev/scripts/CDS2seq.v2.pl -i /mnt/netapp2/Store_csebdjgl/reference_genomes/lynx_lynx_mLynLyn1.2/LYLY1_2A.FA.gff3 -f ancestral_mLynLyn1.2.revcomp.scaffolds.fa

# checking the length of the proteins and cds

awk '/^>/ {if (seq) print seq; printf "%s ", $0; seq=""; next} {seq=seq""$0} END {if (seq) print seq}' LYLY1_2A.FA.pep.fa | sed 's/partial_cds //' | awk '{print $1, length($2)}' > pep_length_LYLY.txt
zcat /mnt/netapp2/Store_csebdjgl/reference_genomes/lynx_lynx_mLynLyn1.2/LYLY1_2A.pep.fa.gz | awk '/^>/ {if (seq) print seq; printf "%s ", $0; seq=""; next} {seq=seq""$0} END {if (seq) print seq}' | sed 's/partial_cds //'  | awk '{print $1, length($2)}' > pep_length_LYLY_ref.txt


awk '/^>/ {if (seq) print seq; printf "%s ", $0; seq=""; next} {seq=seq""$0} END {if (seq) print seq}' LYLY1_2A.FA.cds.fa | sed 's/partial_cds //' | awk '{print $1, length($2)}' > cds_length_LYLY.txt
zcat /mnt/netapp2/Store_csebdjgl/reference_genomes/lynx_lynx_mLynLyn1.2/LYLY1_2A.cds.fa.gz | awk '/^>/ {if (seq) print seq; printf "%s ", $0; seq=""; next} {seq=seq""$0} END {if (seq) print seq}' | sed 's/partial_cds //'  | awk '{print $1, length($2)}' > cds_length_LYLY_ref.txt

# when we are finished with the inspections, we can either remove or move the files to a checks folder
mkdir checks
mv pep_length_LYLY.txt pep_length_LYLY_ref.txt cds_length_LYLY.txt cds_length_LYLY_ref.txt checks/

# check if there are any stop codons in the middle of the translated protein
grep -E '[^ ]\*[^ ]' LYLY1_2A.FA.pep.fa | less -S

# count the number of such occurences 
grep -o -E '[^ ]\*[^ ]' LYLY1_2A.FA.pep.fa | wc -l # 28

# Generate the longest peptide file

awk '/^>/ {if (seq) print seq; printf "%s ", $0; seq=""; next} {seq=seq""$0} END {if (seq) print seq}' LYLY1_2A.FA.pep.fa > LYLY1_2A_tab.pep.fa

for i in $(sed 's/\(>LYLY1_2A[0-9]*\)P[0-9]*/\1/' LYLY1_2A_tab.pep.fa | cut -f1 -d' ' | sort -u) ; do 
    grep "${i}" LYLY1_2A_tab.pep.fa | sed 's/partial_cds //' | awk '(length($2) > max_length) { max_length = length($2); max_line = $0 } END { print max_line }'  >> LYLY1_2A_tab.longest.pep.fa
done

awk '{print $1"\n"$2}' LYLY1_2A_tab.longest.pep.fa > LYLY1_2A.longest.pep.fa

# check the number of proteins that don't have a stop codon at the end
awk '/^>/ {if (seq) print seq; printf "%s ", $0; seq=""; next} {seq=seq""$0} END {if (seq) print seq}' LYLY1_2A.FA.pep.fa | grep -vE '[^ ]\*$' | wc -l # 749
awk '/^>/ {if (seq) print seq; printf "%s ", $0; seq=""; next} {seq=seq""$0} END {if (seq) print seq}' LYLY1_2A.longest.pep.fa | grep -vE '[^ ]\*$' | wc -l # 701

grep -o -E '[^ ]\*[^ ]' LYLY1_2A.longest.pep.fa


# Get the longest transcript

# For the transcripts file
for i in $(grep ">"  LYLY1_2A.longest.pep.fa | sed 's/>//'); do
    awk -v id="${i}" -v RS='>' -v ORS='' '$1 == id {print ">"$0}' LYLY1_2A.FA.cds.fa
done > LYLY1_2A.longest.cds.fa

# The longest cds file has the product name, I need to change it to the actual transcript. This will be done once I have the gff file with the longest isoform only.

# Edit the SnpEff config file to include the new ancestral genome
nano /mnt/netapp2/Store_csebdjgl/agonev/software/snpEff/snpEff.config
# Add the following line:
# Lynx_lynx, version mLynLyn1.2 ancestral
LYLY1_2A_anc.genome : Eurasian lynx ancestral

# Make a directory for the ancestral genome
mkdir -p /mnt/netapp2/Store_csebdjgl/agonev/software/snpEff/data/LYLY1_2A_anc

# ## To generate the new GFF3 file, I use the custom script gff_longest_generation: (but be in the ref genome directory)
# sbatch -t 04:00:00 -c 1 --mem=1G /mnt/netapp2/Store_csebdjgl/agonev/scripts/gff_longest_generation.sh 
# # We ran it interactively in a compute node. 
# Doesn't really work well!

# Because there are cases with multiple transcripts that translate to the same longest product, we ran another version of the script which only keeps the first one. 
sbatch -t 04:00:00 -c 1 --mem=1G /mnt/netapp2/Store_csebdjgl/agonev/scripts/gff_longest_single_transcript.sh
# JOB ID: 8030850
 
# copy and rename the gff3 file
cp /mnt/netapp2/Store_csebdjgl/agonev/data/polarization/LYLY1_2A.longest.isoform_single_transcript.gff3 /mnt/netapp2/Store_csebdjgl/agonev/software/snpEff/data/LYLY1_2A_anc/genes.gff

####### Checking the gff3 files #######

# to check that the new gff3 files have the same number of genes as the cds file
awk '$3 == "transcript"' LYLY1_2A.longest.isoform.gff3 | wc -l # 24996 - should have more than products because some products are coded by multiple transcripts
awk '$3 == "gene"' LYLY1_2A.longest.isoform.gff3 | wc -l # 22525 - has more genes than expected
awk '$3 == "gene"' LYLY1_2A.longest.isoform.gff3 | cut -f9 | sort | uniq | wc -l # 22429 - has the expected number of genes.
# also check the unique products in the same file, and then compare with the single transcript file


# check all the transcripts that don't have a gene in the line before them
grep -w -B1 "transcript" LYLY1_2A.longest.isoform.gff3 | grep -w -A1 "CDS" | less -S
cut -f3 LYLY1_2A.longest.isoform.gff3 | grep -B1 "transcript" | grep "CDS" | wc -l # 2472 (24996 - 2472 = 22524, which ~corresponds to the # genes observed)

# get the genes that are dupicate
awk '$3 == "gene"' LYLY1_2A.longest.isoform.gff3 | sort | uniq -cd | less -S # 86 duplicates

# Filter the transcript product tsv file to only contain the longest products found in LYLY1_2A.longest.pep.fa . The products are in the second column of the tsv file
awk -F"\t" 'NR==FNR{a[$1];next} $1 in a' <(grep ">" LYLY1_2A.longest.pep.fa | sed 's/>//') /mnt/netapp2/Store_csebdjgl/agonev/data/polarization/LYLY1_2A.FA_longestisoform.tsv > /mnt/netapp2/Store_csebdjgl/agonev/data/polarization/LYLY1_2A.FA_longestisoform_filtered.tsv

# We proceed with the single transcript only!

#######################

# Copy the ancestral genome to the snpEff data folder
cp /mnt/netapp2/Store_csebdjgl/agonev/data/polarization/ancestral_mLynLyn1.2.revcomp.scaffolds.fa /mnt/netapp2/Store_csebdjgl/agonev/software/snpEff/data/LYLY1_2A_anc/sequences.fa

# Copy the longest transcript file to the snpEff data folder
cp /mnt/netapp2/Store_csebdjgl/agonev/data/polarization/LYLY1_2A.longest.cds.fa /mnt/netapp2/Store_csebdjgl/agonev/software/snpEff/data/LYLY1_2A_anc/cds.fa

# Copy the longest peptide file to the snpEff data folder
cp /mnt/netapp2/Store_csebdjgl/agonev/data/polarization/LYLY1_2A.longest.pep.fa /mnt/netapp2/Store_csebdjgl/agonev/software/snpEff/data/LYLY1_2A_anc/protein.fa

# Before building the database, I need to change the cds file to have the transcript name instead of the product name

cd /mnt/netapp2/Store_csebdjgl/agonev/software/snpEff/data/LYLY1_2A_anc 


# Step 1: Create a mapping file from IDs in cds.fa (the product names) to transcripts in genes.gff
awk '/^>/ {print substr($1,2)}' cds.fa | grep -wFf - genes.gff | cut -f9 | cut -f4,1 -d';' | sed 's/ID=//; s/product=//' | awk -F '[;]' '{print $2"\t"$1}' > id_to_transcript.tsv

# Step 2: Use awk to replace IDs in cds.fa with transcripts, based on the mapping file
awk 'BEGIN{FS=OFS="\t"} 
    NR==FNR {map[$1]=$2; next} 
    /^>/ {
        split($0, parts, " ", seps); 
        sub(/^>/, "", parts[1]); 
        if(parts[1] in map) 
            print ">" map[parts[1]] seps[1] parts[2]; 
        else 
            print ">" parts[1] seps[1] parts[2]; 
    } 
    !/^>/ {print}' id_to_transcript.tsv cds.fa > cds_replaced.fa

# After checking that everything went well and all the partial cds are retained, we can overwrite the original cds.fa file (because we have it in the polarization folder: LYLY1_2A.longest.cds.fa)
mv cds_replaced.fa cds.fa


# we have to do the same for the protein.fa file, change all the product names to transcript names

awk 'BEGIN{FS=OFS="\t"} 
    NR==FNR {map[$1]=$2; next} 
    /^>/ {
        split($0, parts, " ", seps); 
        sub(/^>/, "", parts[1]); 
        if(parts[1] in map) 
            print ">" map[parts[1]] seps[1] parts[2]; 
        else 
            print ">" parts[1] seps[1] parts[2]; 
    } 
    !/^>/ {print}' id_to_transcript.tsv protein.fa > protein_replaced.fa

# we can overwrite the original protein.fa file (because we have it in the polarization folder)
mv protein_replaced.fa protein.fa

# build the database using the following command in a compute node
compute -c 6 --mem 24

java -jar $EBROOTSNPEFF/snpEff.jar build -d -gff3 -v LYLY1_2A_anc -c /mnt/netapp2/Store_csebdjgl/agonev/software/snpEff/snpEff.config -dataDir /mnt/netapp2/Store_csebdjgl/agonev/software/snpEff/data | tee -a snpeff_build_database.out

#check the database
java -Xmx16g -jar $EBROOTSNPEFF/snpEff.jar dump LYLY1_2A_anc. -c /mnt/netapp2/Store_csebdjgl/agonev/software/snpEff/snpEff.config | less -S 

## Running SNPeff

SNPeff_dir="/mnt/netapp2/Store_csebdjgl/agonev/software/snpEff"
vcf_dir="/mnt/netapp2/Store_csebdjgl/lynx_genome/lynx_data/mLynLyn1.2_ref_vcfs/downsampling/polarization"

java -Xmx16g -jar $EBROOTSNPEFF/snpEff.jar LYLY1_2A_anc \
 -v -c $SNPeff_dir/snpEff.config \
 -s "$SNPeff_dir/data/LYLY1_2A_anc/snpeff_c_ll_mLynLyn1.2_ref_downsampled.filter4.rd.miss30_polarized.html" \
 -csvStats "$SNPeff_dir/data/LYLY1_2A_anc/snpeff_c_ll_mLynLyn1.2_ref_downsampled.filter4.rd.miss30_polarized.csv" \
 $vcf_dir/c_ll_mLynLyn1.2_ref_downsampled.filter4.rd.miss30_polarized.vcf > $vcf_dir/c_ll_mLynLyn1.2_ref_downsampled.filter4.rd.miss30_polarized_annotated.vcf


# Build the amino acid change file
cd /mnt/netapp2/Store_csebdjgl/lynx_genome/lynx_data/mLynLyn1.2_ref_vcfs/downsampling/polarization

grep -v "#" c_ll_mLynLyn1.2_ref_downsampled.filter4.rd.miss30_polarized_annotated.vcf | grep "missense_variant" | cut -f7,11 -d'|' | sed 's/|/\t/;s/p.//' > transcript_mutation.tsv

# Check which are the 34 proteins with errors. (and if they overlap with the 21 proteins that have a stop codon in the middle)