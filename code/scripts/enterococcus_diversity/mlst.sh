#!/bin/bash

#SBATCH --output=log/enterococcus_diversity/mlst.out
#SBATCH --error=log/enterococcus_diversity/mlst.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --partition=cpu

# Activate the conda environment with mlst installed
source $(conda info --base)/etc/profile.d/conda.sh
conda activate mlst

# Define paths
GENOMES_DIR="analyses/enterococcus_diversity/genomes/mags"
OUTPUT_FILE="analyses/enterococcus_diversity/genomes/mlst/mlst_2002.csv"
mkdir -p $(dirname ${OUTPUT_FILE})

# Run mlst for all .fna files in the directory
mlst ${GENOMES_DIR}/*.fna --csv --threads 4 > ${OUTPUT_FILE}