#!/bin/bash -e
#SBATCH --job-name=generate_correlations
#SBATCH --account=uoa00348
#SBATCH --time=48:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=124
#SBATCH --partition=milan
#SBATCH --output=slurm_out/%x.%j.%A_%a.out
#SBATCH --error=slurm_err/%x.%j.%A_%a.err
#SBATCH --array=0-1
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jian.sheng.boey@auckland.ac.nz

module purge
module load \
    GSL/2.7-GCC-12.3.0 \
    R/4.3.2-foss-2023a

a=(filt_tpm.*.csv.gz)
in_file=${a[${SLURM_ARRAY_TASK_ID}]}
b=$(basename ${in_file} .csv.gz | sed -e 's/filt_tpm.//g')

Rscript 2.2a.coexpression_generate_correlations.R \
    --in_file ${in_file} \
    --out_file ${b}.spearman_edgelist.csv \
    --threads ${SLURM_CPUS_PER_TASK} \
    --permutations 2000 \
    --cutoff 0.8 \
    --tmp_dir ${SLURM_JOB_ID}.${b} \
    --keep_tmp_dir
