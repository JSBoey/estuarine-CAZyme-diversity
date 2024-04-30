#!/bin/bash -e
#SBATCH --job-name=signalp6-gpu
#SBATCH --account=uoa00348
#SBATCH --time=6:00:00
#SBATCH --mem=20G
#SBATCH --cpus-per-task=12
#SBATCH --gpus-per-node=1
#SBATCH --output=slurm_out/%x.%j.out
#SBATCH --error=slurm_err/%x.%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jian.sheng.boey@auckland.ac.nz

# Modules
module purge
module load \
  SeqKit/2.2.0 \
  Python/3.9.9-gimkl-2020a \
  CUDA/12.0.0

# Environment
source /nesi/project/uoa02469/Software/signalp6_fast-gpu/bin/activate

# Directories
INDIR=data/2.orf_prediction
OUTDIR=data/3.annotation/signalp6
TMPDIR=tmp/signalp6_$SLURM_JOB_ID

mkdir -p {$OUTDIR,$TMPDIR}

# Variables
PREDFILE=$INDIR/allbins_pred.faa
FORMAT=none
ORGANISM=other
MODE=fast
BATCH_SIZE=500

# Run
## Split input into multiple files
seqkit split2 -s 50000 $PREDFILE -O $TMPDIR

## Run SignalP6 on part files
for i in $TMPDIR/*.faa; do

  # Get part name
  part=$(basename $i .faa)

  # Create part result directory
  mkdir -p $TMPDIR/$part

  # Run
  signalp6 \
    --fastafile $i \
    --output_dir $TMPDIR/$part \
    --format $FORMAT \
    --organism $ORGANISM \
    --mode $MODE \
    --bsize $BATCH_SIZE \
    --write_procs $SLURM_CPUS_PER_TASK

  # Add part as prefix
  for j in $TMPDIR/$part/*; do

    # Get result basename
    filename=$(basename $j)

    # Move results to OUTDIR with part prefix
    mv $j $OUTDIR/${part}_${filename}

  done

done
