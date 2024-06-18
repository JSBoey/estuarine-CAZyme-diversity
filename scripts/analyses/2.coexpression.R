# Co-expression analyses

# 0. Set-up
## Graph based libraries are important here
library(igraph)
library(tidygraph)
library(ggraph)

# Need to subset the meta-transcriptome data for processing at HPC side
g <- ANNOTATION %$%
  (end - start + 1)
names(g) <- ANNOTATION$node

TPM_ns <- lapply(c("water" = "Filt", "sediment" = "Sed"), \(j) {
  m <- split_matrix(calc_TPM(COUNTS$WTS, g),
                    columns = j,
                    remove_zero_rowsum = TRUE)
  m[rowSums(ifelse(m > 0, 1, 0)) > 1, ]
})

for (i in seq_along(TPM_ns)) {
  fwrite(x = as.data.table(TPM_ns[[i]], keep.rownames = "node"),
         file = glue("sandbox/{names(TPM_ns)[i]}.tpm_nosingletons.csv"))
  rm(i)
}

# 1. Parallel correlations
# This is done in HPC with the following script
# Working directory: /nesi/nobackup/uoa00348/boey/large-correlations
# 
# #!/bin/bash -e
# #SBATCH --job-name=parallel_correlations
# #SBATCH --account=uoa00348
# #SBATCH --time=24:00:00
# #SBATCH --mem=64G
# #SBATCH --cpus-per-task=128
# #SBATCH --partition=milan
# #SBATCH --array=0-1
# #SBATCH --output=%x.%j.%A_%a.out
# #SBATCH --error=%x.%j.%A_%a.err
# 
# module purge
# module load \
# GSL/2.7-GCC-12.3.0 \
# R/4.3.2-foss-2023a
# 
# a=(*.tpm.csv)
# infile=${a[$SLURM_ARRAY_TASK_ID]}
# base=${infile%%.tpm.csv}
# outfile=${base}.permcor.csv
# 
# Rscript 0.parallel_correlations.R \
#   --in-file ${infile} \
#   --out-file ${outfile} \
#   --method spearman \
#   --threads $SLURM_CPUS_PER_TASK 
# module purge
# module pigz
#
# pigz -p $SLURM_CPUS_PER_TASK ${outfile}
