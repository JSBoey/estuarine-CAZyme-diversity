#!/bin/bash -e
#SBATCH --job-name=bt2map
#SBATCH --account=uoa00348
#SBATCH --time=03:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=36
#SBATCH --output=slurm_out/%x.%j.%A_%a.out
#SBATCH --error=slurm_err/%x.%j.%A_%a.err
#SBATCH --array=0-65
#SBATCH --partition=milan

# Read mapping
module purge
module load Bowtie2/2.5.4-GCC-12.3.0 SAMtools/1.19-GCC-12.3.0

a=(data/w{g,t}s/*.R1.fastq.gz)
r1=${a[$SLURM_ARRAY_TASK_ID]}
r2=${r1/R1/R2}
base=$(basename ${r1} .R1.fastq.gz)
outdir=data/read_mapping/bam
index=$(dirname ${outdir})/index/drep_98ani_50comp_5cont
minin=0
maxin=500

# Modify insert size parameters
if grep -q "WGS" ${base}; then
    minin=200
    maxin=800
fi

bowtie2 --sensitive \
    -x ${index} \
    -p ${SLURM_CPUS_PER_TASK} \
    -I ${minin} -X ${maxin} \
    -1 ${r1} -2 ${r2} | \
    samtools sort -@ $SLURM_CPUS_PER_TASK \
        -o ${outdir}/${base}.bam -
