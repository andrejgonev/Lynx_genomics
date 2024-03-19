#!/bin/bash

# Input file
input_file="input_data.txt"

# Function to print sample information to a separate file
print_sample_info() {
    local sample_name="$1"
    local fastq_ids="$2"
    local fastq_r1s="$3"
    local fastq_r2s="$4"
    local in_fastq_folder="$5"
    local output_file="$sample_name.mLynLyn.alignment.yaml"

    echo "sample_name: $sample_name" > "$output_file"
    echo "fastq_id: '$fastq_ids'" >> "$output_file"
    echo "fastq_r1: '$fastq_r1s'" >> "$output_file"
    echo "fastq_r2: '$fastq_r2s'" >> "$output_file"
    echo "in_fastq_folder: $in_fastq_folder" >> "$output_file"
}

# Initialize variables
prev_sample_name=""
prev_in_fastq_folder=""

# Read input file line by line
while read -r sample_name fastq_id1 fastq_r1_1 fastq_r2_1 in_fastq_folder; do
    # Check if sample_name changed
    if [[ "$sample_name" != "$prev_sample_name" ]]; then
        # Output data for the previous sample, if it's not the first iteration
        if [ -n "$prev_sample_name" ]; then
            print_sample_info "$prev_sample_name" "$fastq_ids" "$fastq_r1s" "$fastq_r2s" "$prev_in_fastq_folder"
        fi
        # Initialize variables for the new sample
        fastq_ids="$fastq_id1"
        fastq_r1s="$fastq_r1_1"
        fastq_r2s="$fastq_r2_1"
        prev_sample_name="$sample_name"
        prev_in_fastq_folder="$in_fastq_folder"
    else
        # Append data to the variables
        fastq_ids+=" $fastq_id1"
        fastq_r1s+=" $fastq_r1_1"
        fastq_r2s+=" $fastq_r2_1"
    fi
done < "$input_file"

# Output data for the last sample
if [ -n "$prev_sample_name" ]; then
    print_sample_info "$prev_sample_name" "$fastq_ids" "$fastq_r1s" "$fastq_r2s" "$prev_in_fastq_folder"
fi

echo "Grouped data has been written to individual .mLynLyn.alignment.yaml files"
