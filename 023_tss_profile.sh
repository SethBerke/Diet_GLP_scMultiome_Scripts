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
#SBATCH --job-name=tss_profile
#SBATCH --output=tss_profile_%j.out

source ~/miniforge3/etc/profile.d/conda.sh
conda activate multiome

echo "Job ID: ${SLURM_JOB_ID}"
echo "Started: $(date)"
echo "Node: $(hostname)"
echo "-------------------------------------------"

python /Genomics/pritykinlab/seth/Diet_WL_scMultiome/Diet_GLP_scMultiome_Scripts/023_tss_profile.py

echo "-------------------------------------------"
echo "Complete: $(date)"