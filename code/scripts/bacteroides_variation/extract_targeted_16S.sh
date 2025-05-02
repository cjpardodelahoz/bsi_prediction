#!/bin/bash

#SBATCH --array=1-395%20
#SBATCH --mem-per-cpu=8G
#SBATCH -c 2
#SBATCH --time=24:00:00
#SBATCH --error=log/bacteroides_variation/extract_targeted_16S_%A_%a.err
#SBATCH --output=log/bacteroides_variation/extract_targeted_16S_%A_%a.out
#SBATCH --partition=cpu

# Activate conda
source $(conda info --base)/etc/profile.d/conda.sh
conda activate hmmerseqkit

# Set paths
ASSEMBLY_DIR="analyses/bacteroides_variation/genomes/kraken_binning/bacteroides_metaspades"
OUTPUT_DIR="analyses/bacteroides_variation/16S/targeted"
HMM_PROFILE="hmms/16S/RF00177_16S.hmm"
EXTRACT_SCRIPT="software/custom/nhmmer_extract.sh"

# Sequence header prefix
HEADER_PREFIX="targeted"

# Get the sample name for the current task
READS_DIR="data/reads/yan_sd_2022"
sample_name=$(ls ${READS_DIR} | sed -n ${SLURM_ARRAY_TASK_ID}p)

# Create output directory for the sample
mkdir -p ${OUTPUT_DIR}/${sample_name}

# Input and output files
INPUT_FASTA="${ASSEMBLY_DIR}/${sample_name}/contigs.fasta"
OUTPUT_FASTA="${OUTPUT_DIR}/${sample_name}/16S_sequences.fasta"
OUTPUT_HITS="${OUTPUT_DIR}/${sample_name}/16S_hits.tbl"

# Run the 16S extraction script
sh ${EXTRACT_SCRIPT} -p ${HMM_PROFILE} \
                     -e 0 \
                     -i ${INPUT_FASTA} \
                     -s ${HEADER_PREFIX}_${sample_name} \
                     -o ${OUTPUT_FASTA} \
                     -t ${OUTPUT_HITS} \
                     -c 2