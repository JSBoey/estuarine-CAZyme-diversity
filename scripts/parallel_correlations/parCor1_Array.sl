#!/bin/bash -e
#SBATCH --job-name=parCor1_Array
#SBATCH --account=ga02676
#SBATCH --time=2:00:00
#SBATCH --mem=200G
#SBATCH --partition=milan
#SBATCH --cpus-per-task=64
#SBATCH --array=2
#SBATCH --output=slurm_out/%x.%j.out
#SBATCH --error=slurm_err/%x.%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jian.sheng.boey@auckland.ac.nz

module load R/4.3.1-gimkl-2022a

array=(data/*.one_third.tsv)
input=${array[$SLURM_ARRAY_TASK_ID]}

Rscript parCor_1.R ${input} results/parCor/
