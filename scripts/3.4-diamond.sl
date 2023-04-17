#!/bin/bash -e
#SBATCH --job-name=dmnd
#SBATCH --account=uoa00348
#SBATCH --time=3:00:00
#SBATCH --mem=24G
#SBATCH --cpus-per-task=32
#SBATCH --output=slurm_out/%x.%j.%A_%a.out
#SBATCH --error=slurm_err/%x.%j.%A_%a.err
#SBATCH --array=0-2
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jian.sheng.boey@auckland.ac.nz

# Modules
module purge
module load DIAMOND/2.1.1-GCC-11.3.0

# Directories
INDIR=data/2.orf_prediction
OUTDIR=data/3.annotation
DBDIR=data/0.db

# Variables
DBARRAY=(${DBDIR}/*.dmnd)
DB=${DBARRAY[$SLURM_ARRAY_TASK_ID]}
DBBASE=$(basename ${DB} | sed -E 's/([A-Za-z]+).*/\1/g')

INFILE=${INDIR}/allbins_pred.faa
INBASE=$(basename "${INFILE}" .faa)
OUTFILE=${OUTDIR}/${INBASE}.${DBBASE}.tsv

## Software variables
FORMAT='6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp scovhsp ppos'
MAXTARGET=0

# Run
printf "[%s] diamond blastp against %s\n" "$(date)" "${DBBASE}"

diamond blastp \
  --very-sensitive \
  --threads ${SLURM_CPUS_PER_TASK} \
  --db ${DB} \
  --query ${INFILE} \
  --max-target-seqs ${MAXTARGET} \
  --outfmt ${FORMAT} \
  --out ${OUTFILE}

