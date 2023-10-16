#!/bin/bash -e
#SBATCH --job-name=parCor2_water
#SBATCH --account=ga02676
#SBATCH --time=30:00
#SBATCH --mem=50G
#SBATCH --partition=milan
#SBATCH --cpus-per-task=64
#SBATCH --array=0-600
#SBATCH --output=slurm_out/%j.%x.%A_%a.out
#SBATCH --error=slurm_err/%j.%x.%A_%a.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jian.sheng.boey@auckland.ac.nz

module load R/4.3.1-gimkl-2022a

input1=data/water.tpm.one_third.tsv
array=($(ls -p tmp/parCor/water.tpm.one_third/ | grep -v /))
input2=${array[$SLURM_ARRAY_TASK_ID]}

printf "This is job %d of 1273\n" $(($SLURM_ARRAY_TASK_ID + 1))

Rscript parCor_2.R ${input1} tmp/parCor/water.tpm.one_third/${input2}
