#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=16GB
#SBATCH --time=24:00:00
#SBATCH --partition=pritykinlab,main
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=sb5023@princeton.edu
#SBATCH --job-name=sinto_nk_subset
#SBATCH --output=sinto_nk_subset_%j.out

source ~/miniforge3/etc/profile.d/conda.sh
conda activate multiome

echo "Job ID: ${SLURM_JOB_ID}"
echo "Started: $(date)"
echo "Node: $(hostname)"
echo "-------------------------------------------"

BASE=/Genomics/pritykinlab/seth/Diet_WL_scMultiome
CR_BASE=${BASE}/cr_arc_outputs
BC_DIR=${BASE}/Diet_GLP_scMultiome_Scripts/nk_clean_barcodes_pre_EA2
OUT_DIR=/Genomics/pritykinlab/seth/ATACCompendium/results/atac_alignments

# NK_SFD
echo "Processing NK_SFD..."
sinto filterbarcodes \
    -b ${CR_BASE}/NK_1-SFD/outs/atac_possorted_bam.bam \
    -c ${BC_DIR}/NK_SFD_clean_barcodes.txt \
    --outdir ${OUT_DIR} \
    --barcodetag CB \
    -p 8
mv ${OUT_DIR}/atac_possorted_bam.bam ${OUT_DIR}/NK_SFD.bam
samtools index ${OUT_DIR}/NK_SFD.bam
echo "NK_SFD done: $(date)"

# NK_HFD
echo "Processing NK_HFD..."
sinto filterbarcodes \
    -b ${CR_BASE}/NK_2-HFD/outs/atac_possorted_bam.bam \
    -c ${BC_DIR}/NK_HFD_clean_barcodes.txt \
    --outdir ${OUT_DIR} \
    --barcodetag CB \
    -p 8
mv ${OUT_DIR}/atac_possorted_bam.bam ${OUT_DIR}/NK_HFD.bam
samtools index ${OUT_DIR}/NK_HFD.bam
echo "NK_HFD done: $(date)"

# NK_GLP
echo "Processing NK_GLP..."
sinto filterbarcodes \
    -b ${CR_BASE}/NK_3-HFD_GLP/outs/atac_possorted_bam.bam \
    -c ${BC_DIR}/NK_GLP_clean_barcodes.txt \
    --outdir ${OUT_DIR} \
    --barcodetag CB \
    -p 8
mv ${OUT_DIR}/atac_possorted_bam.bam ${OUT_DIR}/NK_HFD_GLP.bam
samtools index ${OUT_DIR}/NK_HFD_GLP.bam
echo "NK_GLP done: $(date)"

# NK_CR
echo "Processing NK_CR..."
sinto filterbarcodes \
    -b ${CR_BASE}/NK_4-HFD_CR/outs/atac_possorted_bam.bam \
    -c ${BC_DIR}/NK_CR_clean_barcodes.txt \
    --outdir ${OUT_DIR} \
    --barcodetag CB \
    -p 8
mv ${OUT_DIR}/atac_possorted_bam.bam ${OUT_DIR}/NK_HFD_CR.bam
samtools index ${OUT_DIR}/NK_HFD_CR.bam
echo "NK_CR done: $(date)"

echo "-------------------------------------------"
echo "Barcode sanity check:"
for name in NK_SFD NK_HFD NK_GLP NK_CR; do
    bc_file=${BC_DIR}/${name}_clean_barcodes.txt
    if [ -f "$bc_file" ]; then
        count=$(wc -l < "$bc_file")
        echo "  ${name}: ${count} barcodes"
    else
        echo "  ${name}: BARCODE FILE NOT FOUND at ${bc_file}"
    fi
done
echo "-------------------------------------------"
echo "All samples complete: $(date)"
echo "Output BAMs in: ${OUT_DIR}"