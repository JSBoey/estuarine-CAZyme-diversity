#!/bin/bash -e
#SBATCH --job-name=tcdoms
#SBATCH --account=ga02676
#SBATCH --time=4:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=36
#SBATCH --output=slurm_out/%j.%x.out
#SBATCH --error=slurm_err/%j.%x.err
#SBATCH --partition=milan

# Modules
module purge
module load HMMER/3.3.2-GCC-12.3.0

# Variables
db=db/tcDoms/tcDomsGlobal/tcDoms.hmm
seqin=data/dastool_bins.faa
domtblout=data/annotations/dastool_bins.tcdoms.domtblout
out=${domtblout/domtbl/hmm}
z=$(hmmstat ${db} | grep -cv '^#')

if [[ ! -d $(dirname ${domtblout}) ]]; then
    mkdir -p $(dirname ${domtblout})
fi

echo "
This is ${SLURM_JOB_ID} searching ${seqin} against tcDoms.hmm (Z = ${z})
"

# Run
hmmsearch --cpu $SLURM_CPUS_PER_TASK \
          --domtblout $domtblout \
          -o $out \
          -Z $z \
          $db \
          $seqin
