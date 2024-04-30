#!/bin/bash -e
#SBATCH --job-name=bam2cram
#SBATCH --account=uoa00348
#SBATCH --time=30:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=32
#SBATCH --output=slurm_out/%x.%j.%A_%a.out
#SBATCH --error=slurm_err/%x.%j.%A_%a.err
#SBATCH --array=0-65
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jian.sheng.boey@auckland.ac.nz
#SBATCH --partition=milan

# Modules
module purge
module load SAMtools/1.16.1-GCC-11.3.0

# Directories
BAMDIR=data/4.feature_count/bam
CRAMDIR=data/4.feature_count/cram

# Variables
ARR=($BAMDIR/*.bam)
INFILE=${ARR[$SLURM_ARRAY_TASK_ID]}
INBASE=$(basename ${INFILE} .bam)
OUTFILE=${CRAMDIR}/${INBASE}.cram

# Run
samtools view -C \
  -@ $SLURM_CPUS_PER_TASK \
  -T ${CRAMDIR}/allbins_scaffolds.fna \
  -o ${OUTFILE} \
  ${INFILE}
