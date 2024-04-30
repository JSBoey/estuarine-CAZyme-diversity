#!/bin/bash -e
#SBATCH --job-name=prodigal
#SBATCH --account=uoa00348
#SBATCH --time=2:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=2
#SBATCH --output=slurm_out/%x.%j.%A_%a.out
#SBATCH --error=slurm_err/%x.%j.%A_%a.err
#SBATCH --array=0-250
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jian.sheng.boey@auckland.ac.nz

# Reference
printf "%s was called on %s\n" "$(basename $0)" "$(date)"

# Modules
module purge
module load prodigal/2.6.3-GCCcore-7.4.0

# Directories
DIRIN=data/1.bins
DIROUT=data/2.orf_prediction

mkdir -p $DIROUT

# Variables
ARR=($DIRIN/*.fna)
INPUT=${ARR[$SLURM_ARRAY_TASK_ID]}
NAME=$(basename $INPUT .fna)
OUTPUT=${DIROUT}/${NAME}_pred

# Run prodigal
prodigal \
  -i ${INPUT} \
  -f gff \
  -a ${OUTPUT}.faa \
  -d ${OUTPUT}.fna \
  -o ${OUTPUT}.gff \
  -p single