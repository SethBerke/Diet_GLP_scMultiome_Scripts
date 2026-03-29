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
#SBATCH --job-name=cr_frip
#SBATCH --output=cr_frip_%j.out

source ~/miniforge3/etc/profile.d/conda.sh
conda activate multiome

echo "Job ID: ${SLURM_JOB_ID}"
echo "Job Name: ${SLURM_JOB_NAME}"
echo "Started: $(date)"
echo "Node: $(hostname)"
echo "-------------------------------------------"

BASE=/Genomics/pritykinlab/seth/Diet_WL_scMultiome
SCRIPT=${BASE}/Diet_GLP_scMultiome_Scripts/021_cr_frip_comparison.py
OUT=${BASE}/Diet_GLP_scMultiome_Scripts/atac_filtered/cr_frip_comparison.txt

python ${SCRIPT} > ${OUT}

echo "-------------------------------------------"
echo "Complete: $(date)"
echo "Output saved to: ${OUT}"