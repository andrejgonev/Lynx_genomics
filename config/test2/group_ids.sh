#!/bin/bash

# Enable debugging mode
set -x

# Input file
input_file="sample_data.txt"

# Output file
output_file="grouped_data.txt"

# Function to print sample information
print_sample_info() {
    local sample_name="$1"
    local fastq_ids="$2"
    local fastq_r1s="$3"
    local fastq_r2s="$4"
    local in_fastq_folder="$5"

    echo "sample_name: $sample_name" >> "$output_file"
    echo "fastq_id: '$fastq_ids'" >> "$output_file"
    echo "fastq_r1: '$fastq_r1s'" >> "$output_file"
    echo "fastq_r2: '$fastq_r2s'" >> "$output_file"
    echo "in_fastq_folder: $in_fastq_folder" >> "$output_file"
    echo >> "$output_file"
}

# Check if input file exists
if [ ! -f "$input_file" ]; then
    echo "Error: Input file '$input_file' not found."
    exit 1
fi

# Read input file line by line
while read -r sample_name fastq_id fastq_r1 fastq_r2 in_fastq_folder; do
    # Append data to arrays corresponding to the sample_name
    sample_info["${sample_name}_fastq_ids"]+=" $fastq_id"
    sample_info["${sample_name}_fastq_r1s"]+=" $fastq_r1"
    sample_info["${sample_name}_fastq_r2s"]+=" $fastq_r2"
    sample_info["${sample_name}_in_fastq_folder"]=$in_fastq_folder
done < "$input_file"

# Check if output file can be written
if [ -e "$output_file" ]; then
    echo "Error: Output file '$output_file' already exists. Please remove it or choose a different name."
    exit 1
fi

# Write grouped data to the output file
for sample_name_key in "${!sample_info[@]}"; do
    if [[ $sample_name_key == *_fastq_ids ]]; then
        sample_name="${sample_name_key%%_fastq_ids}"
        print_sample_info "$sample_name" "${sample_info[$sample_name_key]}" \
            "${sample_info[${sample_name}_fastq_r1s]}" "${sample_info[${sample_name}_fastq_r2s]}" \
            "${sample_info[${sample_name}_in_fastq_folder]}"
    fi
done

echo "Grouped data has been written to $output_file"

# Disable debugging mode
set +x
