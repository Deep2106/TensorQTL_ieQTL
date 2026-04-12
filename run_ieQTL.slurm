#!/bin/bash
#SBATCH --job-name=ieqtl_tensorqtl
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --time=12:00:00
#SBATCH --output=ieqtl_%j.log
#SBATCH -p compute

# Load your conda environment with tensorqtl
# source activate <TENSORQTL_ENV>

cd /path/to/LATEST_IEQTL_RUN_March2026

python main_ieqtl_430.py

echo "Job finished at $(date)"
