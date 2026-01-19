#!/bin/bash

#SBATCH --array=1-2053
#SBATCH --output=log/bacteroides_pul/filter_contigs_isabl1_%A_%a.out
#SBATCH --error=log/bacteroides_variation/filter_contigs_isabl1_%A_%a.err
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --partition=cpu

# Activate the hmmerseqkit environment
source $(conda info --base)/etc/profile.d/conda.sh
conda activate hmmerseqkit

# Get the sample name
sample_name=$(cat misc_files/bacteroides_pul/isabl1_succesful_assemblies.txt | sed -n ${SLURM_ARRAY_TASK_ID}p)

# Define paths
ASSEMBLY_DIR="analyses/bacteroides_pul/metagenomes/assembly/metaspades"

# Define input and output paths
INPUT_FILE="${ASSEMBLY_DIR}/${sample_name}/contigs.fasta"
OUTPUT_FILE="${ASSEMBLY_DIR}/${sample_name}/contigs_1000.fasta"

# Remove previous output file if it exists
if [[ -f "${OUTPUT_FILE}" ]]; then
    rm "${OUTPUT_FILE}"
fi

# Filter contigs >= 1000 bp
if [[ -f "${INPUT_FILE}" ]]; then
    seqkit seq -m 1000 ${INPUT_FILE} -o ${OUTPUT_FILE}
    echo "Filtered contigs for sample ${sample_name} saved to ${OUTPUT_FILE}"
else
    echo "Input file ${INPUT_FILE} does not exist for sample ${sample_name}"
fi