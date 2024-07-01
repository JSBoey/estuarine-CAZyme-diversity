# 2.1. Co-expression analyses - Data subset

# 1. Single- and doubletons are not useful for correlations, remove these 
#    transcripts from data (per environment)
# 2. For any pairs of genes, shared zeroes can deleteriously mask relationships, 
#    and no transformation can prevent that bias. Removal of shared zeroes needs 
#    to be performed on a per-pair basis.
# 3. Genes with a count of 5 or less are considered not active and converted to
#    zero.

# NOTE: This is performed locally.

cutoff <- 5
active_genes <- function(cutoff) {
  m <- COUNTS$WTS
  m[m < cutoff] <- 0
  m
}

senv <- c("water" = "Filt", "sediment" = "Sed")

tpm <- lapply(senv, \(ev) {
  calc_TPM(active_genes(cutoff), get_gene_length()) %>% 
    split_matrix(columns = ev, remove_zero_rowsum = TRUE)
})

# Filter out single- and doubletons
tpm <- lapply(tpm, \(m) {
  pa <- decostand(m, "pa")
  i <- rowSums(pa) > 2
  m[i, ]
})

# Export filtered TPM for calculating correlations on NeSI
walk2(tpm, names(tpm), \(x, nm) {
  fwrite(as.data.table(x, keep.rownames = "node"), 
         file = glue("sandbox/filt_tpm.{nm}.csv.gz"))
})