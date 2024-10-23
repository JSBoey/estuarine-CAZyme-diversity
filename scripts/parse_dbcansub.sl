#!/bin/bash -e
#SBATCH --job-name=parse_dbcansub
#SBATCH --account=ga02676
#SBATCH --time=10:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=2
#SBATCH --output=slurm_out/%j.%A_%a.%x.out
#SBATCH --error=slurm_err/%j.%A_%a.%x.err
#SBATCH --array=0-35
#SBATCH --partition=milan

export PATH=$PATH:$(pwd)/bin

a=(tmp/dbcansubout/*.domtblout)
b=${a[$SLURM_ARRAY_TASK_ID]}
out=${b}_parsed

echo "Parsing ${b}"

hmmsearch_parser.sh ${b} > ${out}
