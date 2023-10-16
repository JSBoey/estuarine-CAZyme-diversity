#!/bin/bash -e
#SBATCH --job-name=parCor3_sed
#SBATCH --account=ga02676
#SBATCH --time=1:00:00
#SBATCH --mem=50G
#SBATCH --partition=milan
#SBATCH --cpus-per-task=64
#SBATCH --array=0-1
#SBATCH --output=slurm_out/%j.%x.%A_%a.out
#SBATCH --error=slurm_err/%j.%x.%A_%a.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jian.sheng.boey@auckland.ac.nz

module load R/4.3.1-gimkl-2022a

input1=data/sediment.tpm.one_third.tsv
indir1=tmp/parCor/$(basename ${input1} .tsv)
array=($(ls -p ${indir1}/ | grep -v /))
input2=${array[$SLURM_ARRAY_TASK_ID]}

Rscript parCor_3.R ${input1} ${input2}
