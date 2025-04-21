#!/bin/bash

# Script to modify FASTA headers for Kraken compatibility across multiple files

# Function to display the help menu
show_help() {
    echo "Usage: $0 -c <csv_file> -d <genome_dir> -o <output_dir>"
    echo
    echo "Options:"
    echo "  -c    Path to the CSV file (no headers, first column: file names, second column: tax IDs)"
    echo "  -d    Path to the directory containing the genome files"
    echo "  -o    Path to the output directory for modified FASTA files"
    echo "  -h    Show this help message and exit"
    echo
    echo "Description:"
    echo "  This script modifies the headers of multiple genome files based on a CSV file."
    echo "  It removes everything after the first whitespace in each header and appends '|kraken:taxid|<taxid>'."
    echo
    echo "Example:"
    echo "  $0 -c genomes.csv -d genomes_dir -o modified_genomes_dir"
}

# Parse command-line arguments
while getopts "c:d:o:h" opt; do
    case $opt in
        c) csv_file="$OPTARG" ;;
        d) genome_dir="$OPTARG" ;;
        o) output_dir="$OPTARG" ;;
        h) show_help; exit 0 ;;
        *) show_help; exit 1 ;;
    esac
done

# Check if all required arguments are provided
if [[ -z "$csv_file" || -z "$genome_dir" || -z "$output_dir" ]]; then
    echo "Error: Missing required arguments."
    show_help
    exit 1
fi

# Ensure the output directory exists
mkdir -p "$output_dir"

# Process each line in the CSV file
while IFS=',' read -r filename taxid; do
    # Find the genome file in the specified directory (with wildcard for extensions)
    genome_file=$(find "$genome_dir" -type f -name "${filename}.*" | head -n 1)

    # Check if the genome file exists
    if [[ -z "$genome_file" ]]; then
        echo "Warning: Genome file for '${filename}' not found in '${genome_dir}'. Skipping..."
        continue
    fi

    # Define the output file path
    output_file="${output_dir}/${filename}.fasta"

    # Process the genome file
    awk -v taxid="$taxid" '
        BEGIN { OFS = "" }
        /^>/ { 
            # Modify the header
            split($0, header_parts, " ")
            print header_parts[1], "|kraken:taxid|", taxid
        } 
        !/^>/ { 
            # Print sequence lines as-is
            print $0 
        }
    ' "$genome_file" > "$output_file"

    echo "Processed: $genome_file -> $output_file"
done < "$csv_file"

echo "All files processed. Modified FASTA files saved in '$output_dir'."