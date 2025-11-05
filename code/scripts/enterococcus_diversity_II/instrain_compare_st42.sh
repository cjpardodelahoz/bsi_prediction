#!/bin/bash

#SBATCH --output=log/enterococcus_diversity_II/instrain_compare_st42.out
#SBATCH --error=log/enterococcus_diversity_II/instrain_compare_st42.err
#SBATCH --time=96:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=512G
#SBATCH --partition=cpu

# Load conda environment
source $(conda info --base)/etc/profile.d/conda.sh
conda activate instrain

# Output directory for inStrain comparison results
OUTPUT_DIR="analyses/enterococcus_diversity_II/strains/instrain/isabl1/st412/instrain_compare"

# inStrain profiles for ST42 genomes
PROFILES=$(find analyses/enterococcus_diversity_II/strains/instrain/isabl1/st412/output -mindepth 2 -maxdepth 2 -type d -name instrain_profile)

# Genome list
GENOME_LIST=$(ls analyses/yan_sd_2022/binning/drep98/dereplicated_genomes/*.fna)

# Run inStrain compare
inStrain compare -i ${PROFILES} \
    -o ${OUTPUT_DIR} \
    --stb $(echo ${GENOME_LIST}) \
    -p 16
