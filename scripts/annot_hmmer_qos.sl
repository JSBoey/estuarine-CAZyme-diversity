#!/bin/bash -e
#SBATCH --job-name=annot_hmmer
#SBATCH --account=ga02676
#SBATCH --time=5:00
#SBATCH --mem=24G
#SBATCH --cpus-per-task=36
#SBATCH --output=slurm_out/%x.%j.out
#SBATCH --error=slurm_err/%x.%j.err
#SBATCH --array=1
#SBATCH --partition=milan
#SBATCH --qos=debug

# Modules
module purge
module load HMMER/3.3.2-GCC-12.3.0

# Array variables
a=(db/*.hmm)
db=${a[$SLURM_ARRAY_TASK_ID]}
db_base=$(basename ${db} .hmm)

seqin=data/test_query.faa
domtblout=data/annotations/dastool_bins.${db_base}.domtblout
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