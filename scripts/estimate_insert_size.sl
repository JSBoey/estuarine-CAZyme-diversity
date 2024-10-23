#!/bin/bash -e
#SBATCH --job-name=insize
#SBATCH --account=uoa00348
#SBATCH --time=02:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=36
#SBATCH --output=slurm_out/%x.%j.%A_%a.out
#SBATCH --error=slurm_err/%x.%j.%A_%a.err
#SBATCH --array=0-65
#SBATCH --partition=milan

a=(data/w{g,t}s/*.R1.fastq.gz)
r1=${a[$SLURM_ARRAY_TASK_ID]}
r2=${r1/R1/R2}
base=$(basename ${r1} .R1.fastq.gz)
outdir=data/read_mapping/insert_size_estimates

# Subsample reads
module purge
module load seqtk/1.4-GCC-11.3.0

sub_n=1000000
seed=1647

seqtk sample -s ${seed} ${r1} ${sub_n} > ${outdir}/${base}.R1.fq
seqtk sample -s ${seed} ${r2} ${sub_n} > ${outdir}/${base}.R2.fq

# Map reads
module purge
module load Bowtie2/2.5.4-GCC-12.3.0 SAMtools/1.19-GCC-12.3.0

index=$(dirname ${outdir})/index/drep_98ani_50comp_5cont

bowtie2 --sensitive -x ${index} -p $SLURM_CPUS_PER_TASK -1 ${r1} -2 ${r2} | \
  samtools sort -@ $SLURM_CPUS_PER_TASK -o ${outdir}/${base}.bam

# Get metrics
samtools stats -@ $SLURM_CPUS_PER_TASK ${outdir}/${base}.bam > ${outdir}/${base}.samstats

module purge
module load R/4.3.2-foss-2023a picard/2.26.10-Java-11.0.4

picard CollectInsertSizeMetrics \
  I=${outdir}/${base}.bam \
  O=${outdir}/${base}.insert_size_metrics.txt \
  H=${outdir}/${base}.insert_size_histogram.pdf
