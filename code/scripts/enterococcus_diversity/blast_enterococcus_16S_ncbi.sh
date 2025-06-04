#!/bin/bash

#SBATCH --mem=16G
#SBATCH -c 8
#SBATCH --time=96:00:00
#SBATCH --error=log/enterococcus_diversity/blast_meta_to_ncbi.err
#SBATCH --output=log/enterococcus_diversity/blast_meta_to_ncbi.out
#SBATCH --partition=cpu

# Activate conda
source $(conda info --base)/etc/profile.d/conda.sh
conda activate hmmerseqkit

# Set paths
QUERY_FASTA="analyses/enterococcus_diversity/16S/compiled/enterococcus_16S.fasta"
CHUNK_DIR="analyses/enterococcus_diversity/16S/compiled/chunks"
RESULTS_DIR="analyses/enterococcus_diversity/16S/compiled/chunk_results"
MERGED_RESULTS="analyses/enterococcus_diversity/16S/compiled/enterococcus_ncbi_blast.txt"

# Create directories for chunks and results
mkdir -p ${CHUNK_DIR}
mkdir -p ${RESULTS_DIR}

# Split the input FASTA file into chunks (e.g., 100 sequences per chunk)
seqkit split -p 108 -O ${CHUNK_DIR} ${QUERY_FASTA}

# Activate conda environment for BLAST
conda activate blast

# Process each chunk
for chunk in ${CHUNK_DIR}/*.fasta; do
    chunk_name=$(basename ${chunk%.*})  # Get the chunk name without extension
    output_file="${RESULTS_DIR}/${chunk_name}_blast.txt"

    # Run BLASTn search for the chunk
    blastn -query ${chunk} \
           -db nt \
           -out ${output_file} \
           -outfmt "6 qseqid sseqid staxid pident length nident mismatch gaps sacc" \
           -remote \
           -max_target_seqs 1

    # Add headers to the chunk result
    HEADER="qseqid\tsseqid\tsscinames\tpident\tlength\tnident\tmismatch\tgaps\tsacc"
    sed -i "1i${HEADER}" ${output_file}

    echo "Processed chunk: ${chunk_name}"
    sleep 30  # Add a sleep to avoid overwhelming the server
done

# Merge all chunk results into a single file
> ${MERGED_RESULTS}  # Initialize the merged results file
for result in ${RESULTS_DIR}/*.txt; do
    tail -n +2 ${result} >> ${MERGED_RESULTS}  # Skip headers from individual files
done

# Add a single header to the merged file
sed -i "1i${HEADER}" ${MERGED_RESULTS}