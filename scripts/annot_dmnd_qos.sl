#!/bin/bash -e
#SBATCH --job-name=annot_dmnd
#SBATCH --account=ga02676
#SBATCH --time=5:00
#SBATCH --mem=24G
#SBATCH --cpus-per-task=36
#SBATCH --output=slurm_out/%x.%j.out
#SBATCH --error=slurm_err/%x.%j.err
#SBATCH --array=2
#SBATCH --partition=milan
#SBATCH --qos=debug

# Modules
module purge
module load DIAMOND/2.1.9-GCC-11.3.0

# Array variables
a=(db/*.faa)
db=${a[$SLURM_ARRAY_TASK_ID]}
db_base=$(basename ${db} .faa)

seqin=data/test_query.faa
out=data/annotations/dastool_bins.${db_base}.b6o
evalue=1e-10
outfmt='6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send qcovhsp scovhsp bitscore corrected_bitscore evalue'

if [[ ! -d $(dirname ${out}) ]]; then
    mkdir -p $(dirname ${out})
fi

# Make database if doesn't exist
if [[ ! -f ${db/faa/dmnd} ]]; then
    echo "Need to build database..."
    diamond makedb --threads $SLURM_CPUS_PER_TASK --db ${db/faa/dmnd} --in ${db}
    echo "Database built"
fi

# CAZymes need more stringent E-value
if echo ${db_base} | grep -vq "tcdb"; then
    evalue=1e-102
fi

echo "
This is ${SLURM_JOB_ID} searching ${seqin} against ${db_base} (E-value < ${evalue})
"

diamond blastp --threads $SLURM_CPUS_PER_TASK --db ${db/faa/dmnd} \
    --evalue ${evalue} --query ${seqin} --max-target-seqs 1 \
    --out ${out} --outfmt ${outfmt}
