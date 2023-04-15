#!/bin/bash -e
#SBATCH --job-name=dbcan
#SBATCH --account=uoa00348
#SBATCH --time=24:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=32
#SBATCH --output=slurm_out/%x.%j.out
#SBATCH --error=slurm_err/%x.%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jian.sheng.boey@auckland.ac.nz

# Modules
module purge
module load HMMER/3.3.2-GCC-11.3.0

# Directories
INDIR=data/2.orf_prediction
OUTDIR=data/3.annotation
DBDIR=/nesi/project/uoa02469/Databases

# Variables 
INFILE=${INDIR}/allbins_pred.faa
INBASE=$(basename ${INFILE} .faa)

## Database for dbCAN2
DB1=${DBDIR}/dbCAN2_v11/dbCAN-HMMdb-V11
## Output for dbCAN2
OUTFILE1=${OUTDIR}/${INBASE}.dbcan.domtbl

## Database for dbCAN-sub
DB2=${DBDIR}/dbCAN-sub_20220811/dbCAN_sub.hmm
## Output for dbCAN-sub
OUTFILE2=${OUTDIR}/${INBASE}.dbcan-sub.domtbl

## dbCAN2 size 
SZDB1=$(hmmstat ${DB1} | tail -n 1 | cut -f 1 -d ' ')

## dbCAN-sub size
SZDB2=$(hmmstat ${DB2} | tail -n 1 | cut -f 1 -d ' ')

# Run
## dbCAN
printf "[%s]\tStart dbCAN annotation. Database size = %s\n" "$(date)" $SZDB1

hmmsearch \
  -Z $SZDB1 \
  --cpu ${SLURM_CPUS_PER_TASK} \
  --domtblout ${OUTFILE1} \
  -o /dev/null \
  $DB1 \
  $INFILE

printf "[%s]\tFinished dbCAN annotation. Domain table saved to %s\n" "$(date)" ${OUTFILE1}

## dbCAN-sub
printf "[%s]\tStart dbCAN-sub annotation. Database size = %s\n" "$(date)" $SZDB2

hmmsearch \
  -Z $SZDB2 \
  --cpu ${SLURM_CPUS_PER_TASK} \
  --domtblout ${OUTFILE2} \
  -o /dev/null \
  $DB2 \
  $INFILE

printf "[%s]\tFinished dbCAN-sub annotation. Domain table saved to %s\n" "$(date)" ${OUTFILE2}
