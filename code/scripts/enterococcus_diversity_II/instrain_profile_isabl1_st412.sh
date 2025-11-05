#!/bin/bash

#SBATCH --array=1-172
#SBATCH --output=log/enterococcus_diversity_II/instrain_st412_%A_%a.out
#SBATCH --error=log/enterococcus_diversity_II/instrain_st412_%A_%a.err
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --partition=cpu

# Activate conda
source $(conda info --base)/etc/profile.d/conda.sh

# Set directories and files
READS_PATHS="misc_files/enterococcus_diversity/all_target_samples_paths.txt"
SAMPLES_LIST="analyses/enterococcus_diversity_II/strains/instrain/isabl1/st412/misc/st42_samples.txt"
REF_FASTA="analyses/enterococcus_diversity/strains/ref_contigs/ref_contigs.fna"
OUTPUT_DIR="analyses/enterococcus_diversity_II/strains/instrain/isabl1/st412/output"

# Get sample name (plain) for current task
sample=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${SAMPLES_LIST})

# Add brackets and single quotes to match the format in the paths file
sample_bracketed="['${sample}']"

# Get first record for this sample from the paths table
record=$(awk -v s="$sample_bracketed" '$2==s{print; exit}' misc_files/enterococcus_diversity/all_target_samples_paths.txt)

# Extract R1 and R2 paths (columns 3 and 4)
R1=$(echo "$record" | cut -f3)
R2=$(echo "$record" | cut -f4)

# Replace storage path if needed
R1=$(echo "$R1" | sed 's|/data/brinkvd/|/data1/collab004/|')
R2=$(echo "$R2" | sed 's|/data/brinkvd/|/data1/collab004/|')

# Create output directory for this sample
SAMPLE_OUTDIR="${OUTPUT_DIR}/${sample}"
MAPPING_DIR="${SAMPLE_OUTDIR}/mapping"
mkdir -p ${MAPPING_DIR}

# Output BAM file
BAM_FILE="${MAPPING_DIR}/${sample}.sorted.bam"

# E. faecium reference contigs
REF_CONTIGS="analyses/enterococcus_diversity/strains/ref_contigs/ref_contig_list.txt"

# Genome list
GENOME_LIST=$(ls analyses/yan_sd_2022/binning/drep98/dereplicated_genomes/*.fna)

# Step 1: Index reference if needed
conda activate bowtie2
if [[ ! -f ${REF_FASTA}.1.bt2 ]]; then
    bowtie2-build ${REF_FASTA} ${REF_FASTA}
fi

# Step 2: Map reads and sort BAM
bowtie2 -x ${REF_FASTA} -1 ${R1} -2 ${R2} -p 16 | samtools view -bS - | samtools sort -o ${BAM_FILE}
samtools index ${BAM_FILE}
conda deactivate

# Step 3: Run inStrain profile
conda activate instrain
inStrain profile ${BAM_FILE} \
    ${REF_FASTA} \
    -o ${SAMPLE_OUTDIR}/instrain_profile \
    --stb $(echo ${GENOME_LIST}) \
    --scaffolds_to_profile ${REF_CONTIGS} \
    -p 16

# Step 4: Clean up BAM files if desired
rm ${BAM_FILE} ${BAM_FILE}.bai