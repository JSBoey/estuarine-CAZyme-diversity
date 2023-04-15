#!/bin/bash -e
#SBATCH --job-name=kofam
#SBATCH --account=uoa00348
#SBATCH --time=24:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=32
#SBATCH --output=slurm_out/%x.%j.out
#SBATCH --error=slurm_err/%x.%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jian.sheng.boey@auckland.ac.nz

# Modules
module purge
module load \
  Ruby/3.0.1-GCC-11.3.0 \
  HMMER/3.3.2-GCC-11.3.0 \
  Parallel/20220922

# Directories
INDIR=data/2.orf_prediction
OUTDIR=data/3.annotation
SOFTDIR=/nesi/project/uoa02469/Software/kofam_scan_v1.3.0
TMPDIR=data/tmp

# Variables
INFILE=${INDIR}/allbins_pred.faa
INBASE=$(basename ${INFILE} .faa)
OUTFILE=${OUTDIR}/${INBASE}.kofam.tsv

## Software variables
kofamscan=${SOFTDIR}/bin/exec_annotation

KOLIST=${SOFTDIR}/db/ko_list
PROFILE=${SOFTDIR}/db/profiles
FORMAT=detail-tsv

# Run
$kofamscan \
  --format=$FORMAT \
  --profile=$PROFILE \
  --ko-list=$KOLIST \
  --cpu=${SLURM_CPUS_PER_TASK} \
  --tmp-dir=$TMPDIR \
  ${INFILE} \
  > ${OUTFILE}
