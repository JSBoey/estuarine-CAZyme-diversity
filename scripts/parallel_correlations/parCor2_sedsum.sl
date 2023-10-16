#!/bin/bash -e
#SBATCH --job-name=parCor2_sedsum
#SBATCH --account=ga02676
#SBATCH --time=30:00
#SBATCH --mem=50G
#SBATCH --partition=milan
#SBATCH --cpus-per-task=64
#SBATCH --array=0
#SBATCH --output=slurm_out/%j.%x.%A_%a.out
#SBATCH --error=slurm_err/%j.%x.%A_%a.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jian.sheng.boey@auckland.ac.nz

module load R/4.3.1-gimkl-2022a

input1=data/sediment_sum.tpm.one_third.tsv
array=(tmp/parCor/sediment_sum.tpm.one_third/ffd81aa4c82c2991090a79c88de0c73e)
input2=${array[$SLURM_ARRAY_TASK_ID]}

Rscript parCor_2.R ${input1} ${input2}
