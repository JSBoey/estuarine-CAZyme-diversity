#!/bin/bash -e
#SBATCH --job-name=bt2map
#SBATCH --account=uoa00348
#SBATCH --time=02:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=36
#SBATCH --output=slurm_out/%x.%j.%A_%a.out
#SBATCH --error=slurm_err/%x.%j.%A_%a.err
#SBATCH --array=0-65
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jian.sheng.boey@auckland.ac.nz
#SBATCH --partition=milan

# Modules
module purge
module load \
  Bowtie2/2.4.5-GCC-11.3.0 \
  SAMtools/1.16.1-GCC-11.3.0

# Directories
WGSDIR=data/0.wgs
WTSDIR=data/0.wts
INDEXDIR=data/4.feature_count/index
BAMDIR=data/4.feature_count/bam

export TMPDIR=data/tmp/bt2map/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}
mkdir -p $TMPDIR

# Variables
ARR=(${WGSDIR}/*.R1.fastq.gz ${WTSDIR}/*.R1.fastq.gz)
INBASE=$(basename ${ARR[$SLURM_ARRAY_TASK_ID]} .R1.fastq.gz)
INDIR=$(dirname ${ARR[$SLURM_ARRAY_TASK_ID]})
INDEX=${INDEXDIR}/allbins_scaffolds
READ1=${INDIR}/${INBASE}.R1.fastq.gz
READ2=${INDIR}/${INBASE}.R2.fastq.gz
OUTBAM=${BAMDIR}/${INBASE}.bam

# Run
bowtie2 -p $SLURM_CPUS_PER_TASK -x $INDEX -1 $READ1 -2 $READ2 \
  | samtools view -@ $SLURM_CPUS_PER_TASK -bS - \
  | samtools sort -@ $SLURM_CPUS_PER_TASK -o $OUTBAM -

# Clean up
rm -rf $TMPDIR
