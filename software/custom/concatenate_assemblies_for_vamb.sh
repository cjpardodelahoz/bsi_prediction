#!/bin/bash

# Usage: ./concatenate_assemblies.sh output.fasta.gz sample1.fasta sample2.fasta ...
# Example: ./concatenate_assemblies.sh concatenated.fna.gz sample1/contigs.fasta sample2/contigs.fasta

# Check if at least two arguments are provided
if [[ $# -lt 2 ]]; then
    echo "Usage: $0 <output_fasta.gz> <assembly1> <assembly2> ..."
    exit 1
fi

# Output file
OUTPUT_FASTA=$1
shift  # Remove the first argument (output file) from the list

# Temporary file for concatenation
TEMP_FASTA=$(mktemp)

# Process each assembly
for ASSEMBLY in "$@"; do
    # Extract the sample name from the file path
    SAMPLE=$(basename $(dirname ${ASSEMBLY}))

    # Rename contigs and append to the temporary file
    sed "s|^>|>${SAMPLE}_C_|" ${ASSEMBLY} >> ${TEMP_FASTA}
done

# Gzip the concatenated FASTA file
gzip -c ${TEMP_FASTA} > ${OUTPUT_FASTA}

# Remove the temporary file
rm ${TEMP_FASTA}