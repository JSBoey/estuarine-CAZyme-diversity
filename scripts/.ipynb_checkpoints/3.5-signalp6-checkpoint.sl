#!/bin/bash -e
#SBATCH --job-name=signalp6-gpu
#SBATCH --account=uoa00348
#SBATCH --time=10:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=12
#SBATCH --gpus-per-node=1
#SBATCH --output=slurm_out/%x.%j.out
#SBATCH --error=slurm_err/%x.%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jian.sheng.boey@auckland.ac.nz

# Modules
module purge
module load Python/3.9.9-gimkl-2020a
module load CUDA/12.0.0

# Environment
source /nesi/project/uoa02469/Software/signalp6_fast-gpu/bin/activate

# Directories
INDIR=data/2.orf_prediction
OUTDIR=data/3.annotation/signalp6

mkdir -p ${OUTDIR}

# Variables
INFILE=${INDIR}/allbins_pred.faa
FORMAT=none
ORGANISM=other
MODE=fast
BATCH_SIZE=500

# Run
signalp6 \
  --fastafile $INFILE \
  --output_dir $OUTDIR \
  --format $FORMAT \
  --organism $ORGANISM \
  --mode $MODE \
  --bsize $BATCH_SIZE \
  --write_procs $SLURM_CPUS_PER_TASK
