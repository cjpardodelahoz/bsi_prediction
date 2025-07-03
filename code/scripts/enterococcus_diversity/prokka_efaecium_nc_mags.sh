#!/bin/bash
#SBATCH --job-name=prokka_efaecium
#SBATCH --output=log/enterococcus_diversity/prokka_efaecium_nc_mags_%A_%a.out
#SBATCH --error=log/enterococcus_diversity/prokka_efaecium_nc_mags_%A_%a.err
#SBATCH --time=024:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --array=1-168

# Conda environment for Prokka
source $(conda info --base)/etc/profile.d/conda.sh
source activate prokka

# Set paths
MAG_LIST_FILE="analyses/enterococcus_diversity/genomes/efaecium_nc.txt"
OUT_DIR="analyses/enterococcus_diversity/genomes/prokka"
REF_GBK="misc_files/entorococcus_diversity/efaecium_DO.gbk"

mkdir -p ${OUT_DIR}

# Get the MAG path
MAG=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $MAG_LIST_FILE)
BASENAME=$(basename "$MAG" .fna)

# Run Prokka
prokka \
  --outdir ${OUT_DIR}/${BASENAME} \
  --force \
  --prefix ${BASENAME} \
  --cpus 4 \
  --addgenes \
  --centre "contig" \
  --compliant \
  --genus Enterococcus \
  --species faecium \
  --usegenus \
  --metagenome \
  --proteins ${REF_GBK} \
  ${MAG}