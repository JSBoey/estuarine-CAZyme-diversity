#!/bin/bash -e
#SBATCH --job-name=ipr5
#SBATCH --account=ga02676
#SBATCH --time=12:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=36
#SBATCH --output=slurm_out/%x.%j.%A_%a.out
#SBATCH --error=slurm_err/%x.%j.%A_%a.err
#SBATCH --array=0-35
#SBATCH --partition=milan

# Modules
module load Java/11.0.4

# Add to path
export PATH=$PATH:$(pwd)/bin

# Array variables
a=(tmp/iprin/*)
input=${a[$SLURM_ARRAY_TASK_ID]}
inbase=$(basename ${input} .fasta)
outbase=tmp/iprout/${inbase}.ipr5
tempdir=./tmp/iprtmp/${SLURM_JOB_ID}

# Run variables
outfmt=xml,tsv
exclappl='ProSiteProfiles,ProSitePatterns' # Need PCRE2 library
cpu=$(( SLURM_CPUS_PER_TASK - 2 ))

if [[ ! -d ${tempdir} ]]; then
    mkdir -p ${tempdir}
fi

if [[ ! -d $(dirname ${outbase}) ]]; then
    mkdir -p $(dirname ${outbase})
fi

interproscan.sh \
    --input ${input} \
    --output-file-base ${outbase} \
    --formats ${outfmt} \
    --excl-applications ${exclappl} \
    --cpu ${cpu} \
    --tempdir ${tempdir} \
    --pathways --goterms
