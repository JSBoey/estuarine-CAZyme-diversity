#!/bin/bash -e
#SBATCH --job-name=graph_metrics_water
#SBATCH --account=uoa00348
#SBATCH --time=48:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --partition=milan
#SBATCH --output=slurm_out/%x.%j.%A_%a.out
#SBATCH --error=slurm_err/%x.%j.%A_%a.err

module purge
module load \
    GSL/2.7-GCC-12.3.0 \
    R/4.3.2-foss-2023a

Rscript 2.5.coexpression_graph_metrics.R \
    --graph water.rho8p05.graphs/global.graphml \
    --out_path water.global_graph.metrics \
    --cpus $SLURM_CPUS_PER_TASK