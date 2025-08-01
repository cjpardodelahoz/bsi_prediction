#!/bin/bash
#SBATCH --output=log/enterococcus_diversity/prokka_ref_contigs.out
#SBATCH --error=log/enterococcus_diversity/prokka_ref_contigs.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=2G

# Conda environment for Prokka
source $(conda info --base)/etc/profile.d/conda.sh
source activate prokka

# Set paths
OUT_DIR="analyses/enterococcus_diversity/strains/ref_contigs"
REF_GBK="misc_files/entorococcus_diversity/efaecium_DO.gbk"
CONTIGS="analyses/enterococcus_diversity/genomes/drep/dereplicated_genomes/s_1044K_1_meta_vamb_single.fna"

# Run Prokka
prokka \
  --outdir ${OUT_DIR}/prokka \
  --force \
  --prefix efaecium_ref \
  --genus Enterococcus \
  --species faecium \
  --cpus 16 \
  --addgenes \
  --compliant \
  --usegenus \
  --metagenome \
  --proteins ${REF_GBK} \
  ${CONTIGS}