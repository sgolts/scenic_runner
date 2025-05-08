#!/bin/bash

#SBATCH --job-name=scenic_pipeline
#SBATCH --partition=gpu     # Adjusted to a common partition type
#SBATCH --account=lab_account
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=logs/scenic_%j.out
#SBATCH --error=logs/scenic_%j.err

mkdir -p ./logs 

# Activate conda environment
CONDA_BASE=$(conda info --base)
source "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate scenic_env

# Input files
INPUT_H5AD="./data/sample_input.h5ad"
TF_FILE="./resources/tfs/allTFs_hg38.txt"
MOTIF_ANNOTATIONS="./resources/motifs/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl"
DB_GLOB="./resources/motifs/*.feather"
MOTIF_DATABASES=( $(ls $DB_GLOB) )

# Output directory
OUTPUT_DIR="./results/scenic_run_$(date +%Y%m%d_%H%M%S)"

# Parameters
MIN_GENES=100
MIN_CELLS=5
NORM_TARGET_SUM="10000"
GRN_WORKERS=12
CTX_WORKERS=12
AUCELL_WORKERS=12
GRN_METHOD="grnboost2"
AUCELL_SEED=123

mkdir -p "${OUTPUT_DIR}"

PYTHON_SCRIPT_PATH="./scenic.py"

echo "Starting SCENIC pipeline run..."
echo "Input H5AD: ${INPUT_H5AD}"
echo "Output Directory: ${OUTPUT_DIR}"

"${PYTHON_SCRIPT_PATH}" \
    --input_h5ad "${INPUT_H5AD}" \
    --tf_file "${TF_FILE}" \
    --motif_annotations "${MOTIF_ANNOTATIONS}" \
    --motif_databases "${MOTIF_DATABASES[@]}" \
    --output_dir "${OUTPUT_DIR}" \
    --min_genes "${MIN_GENES}" \
    --min_cells "${MIN_CELLS}" \
    --norm_target_sum "${NORM_TARGET_SUM}" \
    --grn_method "${GRN_METHOD}" \
    --grn_workers "${GRN_WORKERS}" \
    --ctx_workers "${CTX_WORKERS}" \
    --aucell_workers "${AUCELL_WORKERS}" \
    --aucell_seed "${AUCELL_SEED}"

echo "SCENIC pipeline completed."

conda deactivate