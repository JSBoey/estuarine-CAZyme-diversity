#!/bin/bash -e
#SBATCH --job-name=dbcansub
#SBATCH --account=ga02676
#SBATCH --time=4:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=36
#SBATCH --output=slurm_out/%j.%A_%a.%x.out
#SBATCH --error=slurm_err/%j.%A_%a.%x.err
#SBATCH --array=0-35
#SBATCH --partition=milan

# Modules
module purge
module load HMMER/3.3.2-GCC-12.3.0

# Array variables
a=(tmp/iprin/*)
seqin=${a[$SLURM_ARRAY_TASK_ID]}
inbase=$(basename ${seqin} .fasta)

db=db/dbCAN_sub.hmm
domtblout=tmp/dbcansubout/${inbase}.domtblout
out=${domtblout/domtbl/hmm}
z=$(hmmstat ${db} | grep -cv '^#')

if [[ ! -d $(dirname ${domtblout}) ]]; then
    mkdir -p $(dirname ${domtblout})
fi

echo "
This is ${SLURM_JOB_ID} searching ${seqin} against ${db_base} (Z = ${z})
"

hmmsearch --cpu $SLURM_CPUS_PER_TASK \
          --domtblout $domtblout \
          -o $out \
          -Z $z \
          $db \
          $seqin