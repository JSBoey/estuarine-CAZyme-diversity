#!/bin/bash -e
#SBATCH --job-name=pul_blastp
#SBATCH --account=uoa00348
#SBATCH --time=02:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=36
#SBATCH --output=slurm_out/%x.%j.%x.out
#SBATCH --error=slurm_err/%x.%j.%x.err
#SBATCH --partition=milan

module purge && module load BLAST/2.16.0-GCC-12.3.0

db=db/PUL_20240105_blastdb
query=data/dastool_bins.faa
evalue=0.01
outfmt="6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovhsp"
out=data/annotations/dastool_bins.PUL_20240105.blastp.b6o

blastp -task blastp -query ${query} -max_hsps 1 -db ${db} -outfmt "${outfmt}" -out ${out} -evalue ${evalue} -num_threads ${SLURM_CPUS_PER_TASK}

header=$(echo ${outfmt} | sed -e 's/6 //g' -e 's/ /\t/g')
sed -i "1i${header}" ${out}

awk -v FS="\t" -v OFS="\t" '$3 >= 35 && $15 >= 70' ${out} > ${out}_parsed

module purge && module load pigz
pigz -p8 ${out}
pigz -p8 ${out}_parsed