#!/bin/bash

#SBATCH --array=1-395%10  # Adjust the range (1-395) based on the number of samples
#SBATCH --output=log/bacteroides_variation/bacteroides_filter_contigs_%A_%a.out
#SBATCH --error=log/bacteroides_variation/bacteroides_filter_contigs_%A_%a.err
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --partition=cpu

# Activate the hmmerseqkit environment
source $(conda info --base)/etc/profile.d/conda.sh
conda activate hmmerseqkit

# Define paths
ASSEMBLY_DIR="analyses/bacteroides_variation/genomes/kraken_binning/bacteroides_metaspades"

# Get the sample name for the current task
sample_name=$(ls ${ASSEMBLY_DIR} | sed -n ${SLURM_ARRAY_TASK_ID}p)

# Define input and output paths
INPUT_FILE="${ASSEMBLY_DIR}/${sample_name}/contigs.fasta"
OUTPUT_FILE="${ASSEMBLY_DIR}/${sample_name}/contigs_1000.fasta"

# Filter contigs >= 1000 bp
if [[ -f "${INPUT_FILE}" ]]; then
    seqkit seq -m 1000 ${INPUT_FILE} -o ${OUTPUT_FILE}
    echo "Filtered contigs for sample ${sample_name} saved to ${OUTPUT_FILE}"
else
    echo "Input file ${INPUT_FILE} does not exist for sample ${sample_name}"
fi