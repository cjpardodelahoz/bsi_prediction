#!/bin/bash

#SBATCH --array=1-49%15 # Adjust the range (1-395) based on the number of samples
#SBATCH --output=log/bacteroides_variation/yan_sd_2022_vamb_multi_%A_%a.out
#SBATCH --error=log/bacteroides_variation/yan_sd_2022_vamb_multi_%A_%a.err
#SBATCH --time=96:00:00
#SBATCH --cpus-per-task=24
#SBATCH --mem=32G
#SBATCH --partition=cpu

# Activate the strobealign environment
source $(conda info --base)/etc/profile.d/conda.sh
conda activate strobealign

# Define paths
READS_DIR="data/reads/yan_sd_2022"
ASSEMBLY_DIR="analyses/yan_sd_2022/assembly/metaspades"
BINNING_DIR="analyses/yan_sd_2022/binning"
PATIENT_SAMPLES_DIR="${BINNING_DIR}/patient_samples"
VAMB_OUT_DIR="${BINNING_DIR}/vamb/multi"

# Get the patient ID for the current task
patient_id=$(ls ${PATIENT_SAMPLES_DIR} | sed -n ${SLURM_ARRAY_TASK_ID}p)

# Create output directories
mkdir -p ${VAMB_OUT_DIR}/${patient_id%%.txt}
AEMB_DIR="${VAMB_OUT_DIR}/${patient_id%%.txt}/aemb"
mkdir -p ${AEMB_DIR}

# Step 1: Concatenate assemblies (contigs >1000 bp) for all samples in the patient
PATIENT_SAMPLES=$(cat ${PATIENT_SAMPLES_DIR}/${patient_id})
CONCATENATED_FASTA="${VAMB_OUT_DIR}/${patient_id%%.txt}/contigs.fna.gz"
software/custom/concatenate_assemblies_for_vamb.sh \
    ${CONCATENATED_FASTA} $(for sample in ${PATIENT_SAMPLES}; do echo ${ASSEMBLY_DIR}/${sample}/contigs_1000.fasta; done)

# Step 2: Map reads from all samples to the concatenated assembly using strobealign
for sample in ${PATIENT_SAMPLES}; do
    R1="${READS_DIR}/${sample}/${sample}_R1_all.fastq.gz"
    R2="${READS_DIR}/${sample}/${sample}_R2_all.fastq.gz"
    strobealign -t 16 --aemb ${CONCATENATED_FASTA} ${R1} ${R2} > ${AEMB_DIR}/${sample}.tsv
done

# Step 3: Generate the abundance TSV file
conda activate python3
ABUNDANCE_TSV="${VAMB_OUT_DIR}/${patient_id%%.txt}/abundance.tsv"
python software/vamb/merge_aemb.py ${AEMB_DIR} ${ABUNDANCE_TSV}
conda deactivate

# Step 4: Run VAMB
vamb bin default \
    --outdir ${VAMB_OUT_DIR}/${patient_id%%.txt}/vambout \
    --fasta ${CONCATENATED_FASTA} \
    --abundance_tsv ${ABUNDANCE_TSV} \
    --minfasta 250000 \
    -m 1000 \
    -p 24 \
    -o _C_

