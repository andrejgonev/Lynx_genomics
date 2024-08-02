#!/bin/bash

vcfFile="/mnt/netapp2/Store_csebdjgl/lynx_genome/lynx_data/mLynLyn1.2_ref_vcfs/downsampling/polarization/c_ll_mLynLyn1.2_ref_downsampled.filter4.rd.miss30_polarized_annotated.vcf"
output_file="/mnt/netapp2/Store_csebdjgl/lynx_genome/lynx_data/mLynLyn1.2_ref_vcfs/downsampling/polarization/effect_counts.tsv"

# Extract unique effects
unique_effects=$(grep -v "^#" "$vcfFile" | cut -f8 | cut -d'|' -f2 | sort | uniq)

# Extract individual names (columns 10 onwards) and create the header
header="Individual"
for effect in $unique_effects; do
    header="$header\t$effect"
done
echo -e "$header" > $output_file

# Extract individual names
# individuals=$(grep '^#CHROM' "$vcfFile" | cut -f10- | tr '\t' '\n')

ncols=$(grep -m1 '^#CHROM' "$vcfFile" | tr '\t' '\n' | wc -l)

for i in $(seq 10 $ncols); do
    individual=$(grep -m1 '^#CHROM' "$vcfFile" | cut -f"$i")
    line="$individual"
    echo "Processing individual $individual"

    for effect in $unique_effects; do
        count=$(grep "|$effect|" "$vcfFile" | cut -f"$i" | cut -d':' -f1 | tr '/' '+' | bc | awk '{ sum += $1 } END { print sum }')
        echo $effect $count
        line="$line\t$count"
    done
    echo -e "$line" >> $output_file
done