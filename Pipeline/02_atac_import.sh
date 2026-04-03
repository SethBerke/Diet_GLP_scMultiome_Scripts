#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=32GB
#SBATCH --time=24:00:00
#SBATCH --partition=pritykinlab,main
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=sb5023@princeton.edu
#SBATCH --job-name=atac_import
#SBATCH --output=atac_import_%j.out

source ~/miniforge3/etc/profile.d/conda.sh
conda activate multiome

echo "Job ID: ${SLURM_JOB_ID}"
echo "Job Name: ${SLURM_JOB_NAME}"
echo "Started: $(date)"
echo "Node: $(hostname)"
echo "-------------------------------------------"

BASE=/Genomics/pritykinlab/seth/Diet_WL_scMultiome
BARCODE_DIR=${BASE}/Diet_GLP_scMultiome_Scripts/gex_filtered/barcodes
ATLAS=/Genomics/pritykinlab/seth/ATACCompendium/results/all_atac_peaks.bed
TSS_BED=${BASE}/mm10_tss.bed
OUT_DIR=${BASE}/Diet_GLP_scMultiome_Scripts/atac_gex_filtered_pre_EA2
SCRIPT=${BASE}/Diet_GLP_scMultiome_Scripts/Pipeline/02_atac_import.py

mkdir -p ${OUT_DIR}

declare -A SAMPLES
SAMPLES["NK_1-SFD"]="NK_SFD"
SAMPLES["NK_2-HFD"]="NK_HFD"
SAMPLES["NK_3-HFD_GLP"]="NK_GLP"
SAMPLES["NK_4-HFD_CR"]="NK_CR"

for CR_NAME in "NK_1-SFD" "NK_2-HFD" "NK_3-HFD_GLP" "NK_4-HFD_CR"; do
    SHORT_NAME=${SAMPLES[$CR_NAME]}
    echo "-------------------------------------------"
    echo "Processing ${CR_NAME} -> ${SHORT_NAME}..."
    echo "Started: $(date)"

    python ${SCRIPT} \
        --fragment-file ${BASE}/cr_arc_outputs/${CR_NAME}/outs/atac_fragments.tsv.gz \
        --whitelist     ${BARCODE_DIR}/${SHORT_NAME}_barcodes.txt \
        --atlas         ${ATLAS} \
        --tss-bed       ${TSS_BED} \
        --output        ${OUT_DIR}/${SHORT_NAME}_atac.h5ad \
        --sample-name   ${SHORT_NAME} \
        --threads       8

    echo "Done ${CR_NAME}: $(date)"
done

echo "-------------------------------------------"
echo "All samples complete: $(date)"