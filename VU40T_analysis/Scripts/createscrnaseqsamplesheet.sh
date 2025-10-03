#!/bin/bash

# Parse input arguments
while getopts "i:o:" opt; do
  case $opt in
    i) indir="$OPTARG" ;;
    o) outfile="$OPTARG" ;;
    \?) echo "Invalid option: -$OPTARG" >&2; exit 1 ;;
  esac
done

# Check required arguments
if [[ -z "${indir:-}" || -z "${outfile:-}" ]]; then
  echo "Usage: $0 -i <input_dir> -o <output_csv>"
  exit 1
fi

# Create header
echo "sample,fastq_1,fastq_2" > "$outfile"

# Populate samplesheet
for f1 in "$indir"/*R1_001.fastq.gz; do
  # Get the full path to fastq_2
  f2="${f1/R1_001.fastq.gz/R2_001.fastq.gz}"

  # Extract sample name from filename (customize as needed)
  sample=$(basename "$f1" | sed 's/_GB250724-MW_S[1-4]_L00[12]_R1_001.fastq.gz//')

  # Output CSV line
  echo "$sample,$f1,$f2"
done >> "$outfile"

