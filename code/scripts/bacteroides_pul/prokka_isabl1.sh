#!/bin/bash
#SBATCH --output=log/bacteroides_pul/prokka_%A_%a.out
#SBATCH --error=log/bacteroides_pul/prokka_%A_%a.err
#SBATCH --array=1-2 #1-2077
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH -p cpu
#SBATCH --time=48:00:00

# Activate conda
source $(conda info --base)/etc/profile.d/conda.sh
conda activate prokka

# Directory containing metaspades assemblies (symlinked)
ASSEMBLY_DIR=analyses/bacteroides_pul/metagenomes/assembly/metaspades
# Output directory for Prokka results
OUTDIR=analyses/bacteroides_pul/pul_prediction/prokka
mkdir -p "$OUTDIR"

# Get sample name (with brackets and quotes) for current task
SAMPLE_FILE="misc_files/enterococcus_diversity/all_target_samples_hostdepleted.txt"
sample_bracketed=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMPLE_FILE")
# Remove brackets and single quotes: e.g., ['275C'] -> 275C
SAMPLE=$(echo $sample_bracketed | sed "s/\[\('\([^']*\)'\)\]/\2/")

# Input contigs file
CONTIGS="$ASSEMBLY_DIR/$SAMPLE/contigs_1000.fasta"

# Output directory and prefix for Prokka
PROKKA_OUT="$OUTDIR/$SAMPLE"
PREFIX="$SAMPLE"
LOCUSTAG="$SAMPLE"

# Run Prokka
prokka --kingdom Bacteria --cpus 32 \
    --outdir "$PROKKA_OUT" \
    --prefix "$PREFIX" \
    --addgenes --addmrna --locustag "$LOCUSTAG" \
    "$CONTIGS"
