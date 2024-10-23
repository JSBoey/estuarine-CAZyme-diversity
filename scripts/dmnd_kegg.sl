#!/bin/bash -e
#SBATCH --job-name=dmnd_kegg
#SBATCH --account=ga02676
#SBATCH --time=3:00:00
#SBATCH --mem=24G
#SBATCH --cpus-per-task=36
#SBATCH --output=slurm_out/%x.%j.out
#SBATCH --error=slurm_err/%x.%j.err
#SBATCH --partition=milan

# Modules
module purge
module load DIAMOND/2.1.9-GCC-11.3.0

db=db/kegg_2024.refprok.dmnd
db_base=$(basename ${db} .dmnd)

seqin=data/dastool_bins.faa
out=data/annotations/dastool_bins.${db_base}.b6o

evalue=1e-102
outfmt='6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send qcovhsp scovhsp bitscore corrected_bitscore evalue'

if [[ ! -d $(dirname ${out}) ]]; then
    mkdir -p $(dirname ${out})
fi

diamond blastp --threads $SLURM_CPUS_PER_TASK --db ${db} \
    --evalue ${evalue} --query ${seqin} \
    --out ${out} --outfmt ${outfmt}