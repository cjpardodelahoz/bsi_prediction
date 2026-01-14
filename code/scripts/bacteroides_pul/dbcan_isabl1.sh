#!/bin/bash
#SBATCH --output=log/bacteroides_pul/dbcan_%A_%a.out
#SBATCH --error=log/bacteroides_pul/dbcan_%A_%a.err
#SBATCH --array=46-47 #1-2076
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH -p cpu
#SBATCH --time=24:00:00

# Activate conda
source $(conda info --base)/etc/profile.d/conda.sh
conda activate dbcan

# Get sample name for current task
SAMPLE_FILE="misc_files/bacteroides_pul/isabl1_samples.txt"
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMPLE_FILE")

# Directory containing metagenome assemblies
ASSEMBLY_DIR=analyses/bacteroides_pul/metagenomes/assembly/metaspades

# Input contigs file
CONTIGS="$ASSEMBLY_DIR/$SAMPLE/contigs_1000.fasta"

# Output directory for dbCAN results
OUTDIR=analyses/bacteroides_pul/pul_prediction/dbcan
mkdir -p "$OUTDIR"

# Output directory for this sample
SAMPLE_OUTDIR="$OUTDIR/$SAMPLE"
mkdir -p "$SAMPLE_OUTDIR"

# DBCAN database directory
DB_DIR="/data1/xavierj/carlos/dbs/dbcan"

# Run dbCAN using easy_substrate workflow
run_dbcan easy_substrate \
  --mode prok \
  --input_raw_data "$CONTIGS" \
  --input_gff "${SAMPLE_OUTDIR}/uniInput.gff" \
  --gff_type prodigal \
  --output_dir "$SAMPLE_OUTDIR" \
  --db_dir "$DB_DIR" \
  --threads 16
