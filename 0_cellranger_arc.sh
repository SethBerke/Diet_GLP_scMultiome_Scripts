#!/bin/bash
#SBATCH --mem=128GB
#SBATCH --time=5-00
#SBATCH --nodes=1
#SBATCH --job-name=cr_arc_%j
#SBATCH --mail-user=sb5023@princeton.edu
#SBATCH --mail-type=ALL
#SBATCH --output=cr_arc_%j.out
#SBATCH --array=0-7

# requires a sample_list.txt of sample names; $1 = out_path

# eval "$(conda shell.bash hook)"
# conda activate labenv_r

# set up
echo "Running cellranger-arc on multiome data" # Print status message to log
export PATH=/Genomics/pritykinlab/seth/cellranger-arc-2.1.0:$PATH # Add cellranger-arc 2.1.0 to PATH so it can be called by name
refdata="/Genomics/argo/users/skwalker/refdata-cellranger-arc-mm10-2020-A-2.0.0" # Path to mm10 reference genome

# ${1} = /Genomics/pritykinlab/seth/Diet_WL_scMultiome (passed at sbatch submission)
echo out_path: ${1} # Print the output path to log for confirmation
readarray -t sample_arr < ${1}/sample_list.txt # Print the output path to log for confirmation
echo "This is array task ${SLURM_ARRAY_TASK_ID}, the sample name is ${sample_arr[SLURM_ARRAY_TASK_ID]}" # Read sample_list.txt into an array (one sample name per element)

samplename=${sample_arr[SLURM_ARRAY_TASK_ID]} # Assign this task's sample name to a variable

# [SLURM_ARRAY_TASK_ID] becomes 0, 1, 2, 3, 4, 5, 6, or 7

cellranger-arc count --id=${samplename} \
    --reference=$refdata \
    --libraries=${1}/${samplename}/${samplename}_libraries.csv \
    --localcores=20 \
    --localmem=100 \
    --create-bam=true

echo $samplename DONE!