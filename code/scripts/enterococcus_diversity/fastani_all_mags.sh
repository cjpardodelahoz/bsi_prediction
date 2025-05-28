#!/bin/bash

#SBATCH --output=log/enterococcus_diversity/fastani_all_mags.out
#SBATCH --error=log/enterococcus_diversity/fastani_all_mags.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=512G
#SBATCH --partition=cpu_highmem

# Activate the conda environment with FastANI installed
source $(conda info --base)/etc/profile.d/conda.sh
conda activate fastani

# Define paths
GENOMES_DIR="analyses/yan_sd_2022/binning/compiled"
OUTPUT_FILE="analyses/yan_sd_2022/fastani/fastani_results.tsv"
mkdir -p $(dirname ${OUTPUT_FILE})

# Create a list of all genome files
GENOME_LIST="analyses/yan_sd_2022/fastani/genome_list.txt"
ls ${GENOMES_DIR}/*.fna > ${GENOME_LIST}

# Run FastANI for all-against-all comparison
fastANI --ql ${GENOME_LIST} --rl ${GENOME_LIST} --threads 32 --output ${OUTPUT_FILE}