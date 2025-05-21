#!/bin/bash

#SBATCH --output=log/enterococcus_diversity/fastani.out
#SBATCH --error=log/enterococcus_diversity/fastani.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --partition=cpu

# Activate the conda environment with FastANI installed
source $(conda info --base)/etc/profile.d/conda.sh
conda activate fastani

# Define paths
GENOMES_DIR="analyses/enterococcus_diversity/genomes/mags"
OUTPUT_FILE="analyses/enterococcus_diversity/genomes/fastani/fastani_results.tsv"
mkdir -p $(dirname ${OUTPUT_FILE})

# Create a list of all genome files
GENOME_LIST="analyses/enterococcus_diversity/genomes/fastani/genome_list.txt"
ls ${GENOMES_DIR}/*.fna > ${GENOME_LIST}

# Run FastANI for all-against-all comparison
fastANI --ql ${GENOME_LIST} --rl ${GENOME_LIST} --threads 16 --output ${OUTPUT_FILE}