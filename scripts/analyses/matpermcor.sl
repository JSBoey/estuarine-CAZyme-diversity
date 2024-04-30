#!/bin/bash -e
#SBATCH --job-name=matpermcor
#SBATCH --account=uoa00348
#SBATCH --time=3-00:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=80
#SBATCH --partition=milan
#SBATCH --output=%x.%j.%A_%a.out
#SBATCH --error=%x.%j.%A_%a.err

module purge
module load GSL/2.7-GCC-12.3.0 R/4.3.2-foss-2023a

Rscript matpermcor.R rank_TPM.rds
