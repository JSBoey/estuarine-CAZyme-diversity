#!/bin/bash -e
#SBATCH --job-name=coexp_format
#SBATCH --account=uoa00348
#SBATCH --time=1:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=32
#SBATCH --partition=milan
#SBATCH --array=0-1
#SBATCH --output=%x.%j.%A_%a.out
#SBATCH --error=%x.%j.%A_%a.err

module purge
module load \
    GSL/2.7-GCC-12.3.0 \
    R/4.3.2-foss-2023a

a=(filt_tpm.*.permcor.csv.gz)
in_file=${a[$SLURM_ARRAY_TASK_ID]}
base=${infile%%.csv.gz}
parquet_path=${base}.parquet

Rscript 2.2.coexpression_data_format.R \
    --in_file ${infile} \
    --parquet_path ${parquet_path} \
    --cpus $SLURM_CPUS_PER_TASK
    