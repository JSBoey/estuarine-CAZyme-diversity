#!/bin/bash -e
#SBATCH --job-name=ipr5
#SBATCH --account=uoa00348
#SBATCH --time=2:00:00
#SBATCH --mem=24G
#SBATCH --cpus-per-task=32
#SBATCH --output=slurm_out/%x.%j.%A_%a.out
#SBATCH --error=slurm_err/%x.%j.%A_%a.err
#SBATCH --array=0-9
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jian.sheng.boey@auckland.ac.nz

# Modules
module purge
module load Java/17

# Directories
INDIR=data/tmp
OUTDIR=data/3.annotation

mkdir -p ${OUTDIR}

# Variables
interproscan=bin/interproscan-5.61-93.0/interproscan.sh

ARR=(${INDIR}/allbins_pred.part*.faa)
INFILE=${ARR[$SLURM_ARRAY_TASK_ID]}
FILEPART=$(basename ${INFILE} .faa | sed -E 's/.*(part_[0-9]*)/\1/g')
OUTBASE=${OUTDIR}/allbins_pred.interpro.${FILEPART}
APPL=Pfam,TIGRFAM,CDD
FORMAT=xml,tsv
CPU=$(($SLURM_CPUS_PER_TASK - 2))


# Run InterProScan
$interproscan \
  --applications ${APPL} \
  --cpu ${CPU} \
  --formats ${FORMAT} \
  --input ${INFILE} \
  --output-file-base ${OUTBASE}
