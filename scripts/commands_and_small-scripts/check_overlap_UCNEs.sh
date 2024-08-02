#!/bin/bash

# Read the start and end positions from the two files
file1="file1.csv"
file2="file2.csv"

# Sort the intervals by start position
sort -k1,1n -k2,2n $file1 > sorted_file1.csv
sort -k1,1n -k2,2n $file2 > sorted_file2.csv

# Initialize counters
overlap_count=0
exclusive_file1_count=0
exclusive_file2_count=0

# Function to check for overlaps
check_overlap() {
    local start1=$1
    local end1=$2
    local start2=$3
    local end2=$4

    if [[ $end1 -ge $start2 && $end2 -ge $start1 ]]; then
        return 0  # Overlap
    else
        return 1  # No overlap
    fi
}

# Read the sorted files and check for overlaps
while IFS=, read -r start1 end1; do
    overlap_found=false
    while IFS=, read -r start2 end2; do
        if check_overlap $start1 $end1 $start2 $end2; then
            overlap_found=true
            ((overlap_count++))
            break
        fi
    done < sorted_file2.csv
    if ! $overlap_found; then
        ((exclusive_file1_count++))
    fi
done < sorted_file1.csv

# Count exclusive segments in file2
while IFS=, read -r start2 end2; do
    overlap_found=false
    while IFS=, read -r start1 end1; do
        if check_overlap $start1 $end1 $start2 $end2; then
            overlap_found=true
            break
        fi
    done < sorted_file1.csv
    if ! $overlap_found; then
        ((exclusive_file2_count++))
    fi
done < sorted_file2.csv

# Output the results
echo "Number of overlapping segments: $overlap_count"
echo "Number of exclusive segments in file1: $exclusive_file1_count"
echo "Number of exclusive segments in file2: $exclusive_file2_count"