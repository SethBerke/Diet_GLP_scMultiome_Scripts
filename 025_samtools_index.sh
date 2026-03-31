#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8GB
#SBATCH --time=2:00:00
#SBATCH --partition=pritykinlab,main
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=sb5023@princeton.edu
#SBATCH --job-name=samtools_index
#SBATCH --output=samtools_index_%j.out

echo "Job ID: ${SLURM_JOB_ID}"
echo "Started: $(date)"
echo "-------------------------------------------"

source /Genomics/argo/users/sb5023/miniforge3/etc/profile.d/conda.sh
conda activate multiome

cd /Genomics/pritykinlab/seth/ATACCompendium/results/atac_alignments

samtools index NK_GLP.bam &
samtools index NK_CR.bam &
wait

echo "-------------------------------------------"
echo "Complete: $(date)"
ls -lh NK_*.bai