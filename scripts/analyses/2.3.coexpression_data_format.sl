#!/bin/bash -e
#SBATCH --job-name=format_change
#SBATCH --account=uoa00348
#SBATCH --time=30:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=32
#SBATCH --partition=milan
#SBATCH --output=slurm_out/%x.%j.%A_%a.out
#SBATCH --error=slurm_err/%x.%j.%A_%a.err
#SBATCH --array=0-1

module purge
module load \
    GSL/2.7-GCC-12.3.0 \
    R/4.3.2-foss-2023a

a=(*.spearman_edgelist.csv)
in_file=${a[${SLURM_ARRAY_TASK_ID}]}
b=${in_file%%.csv}

Rscript 2.3.coexpression_data_format.R \
    --in_file ${in_file} \
    --parquet_path ${b}.parquet \
    --cpus $SLURM_CPUS_PER_TASK
