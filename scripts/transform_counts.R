# Transform count tables

# Preamble ----
# This script transforms count tables into Transcript Per Million (TPM) and 
# robust Centred Log-Ratios. In both methods, the counts will be standardised by
# gene length (per kilobase).
# Transformations here are relative to the count matrix alone. Therefore, a 
# per-taxon approach based on Klingenberg & Meinicke will also be applied in the
# interest of experimentation.

# Environment ----
library(vegan)
library(tidyverse)

# Import data ----
count_data <- list(
  "wgs" = "results/WGS_clean_count.tsv.gz",
  "wts" = "results/WTS_clean_count.tsv.gz"
) %>% 
  map(read_tsv)

annotation <- read_tsv("results/allbins_pred.annotation_table.tsv.gz")

# Transform functions ----
trans_RPK <- function(mat, gene_length) {
  sweep(mat, 1, gene_length/1e3, "/")
}

trans_TPM <- function(rpk) {
  sweep(rpk, 2, colSums(rpk)/1e6, "/")
}

trans_rCLR <- function(rpk) {
  decostand(rpk, method = "rclr", MARGIN = 2)
}

# Transform data ----
transform_count <- map(count_data, \(df) {
  # Create matrix
  MAT <- as.matrix(select(df, contains(c("Filt", "Sed"))))
  rownames(MAT) <- df$Geneid
  
  # Standardise counts by per kilobase gene length
  RPK <- trans_RPK(MAT, gene_length = df$Length)
  
  # Transform matrix
  list(
    "TPM" = trans_TPM(RPK),
    "rCLR" = trans_rCLR(RPK)
  )
}) %>% 
  list_flatten()

## Save
saveRDS(transform_count, file = "results/transformed_counts.rds")

# Per-taxon transform (NOT USED) ----
count_per_taxon <- map(count_data, \(df) {
  binned_counts <- left_join(
    select(df, Geneid, Length, starts_with(c("Filt", "Sed"))), 
    select(annotation, node, bin), 
    by = c("Geneid" = "node")
  ) %>% 
    group_by(bin) %>% 
    nest()
  
  binned_matrices <- binned_counts %>% 
    mutate(
      GENE_LENGTH = map(data, \(x) set_names(x$Length, x$Geneid)),
      MAT = map(data, \(x) {
        select(x, -Length) %>% 
          column_to_rownames("Geneid") %>% 
          as.matrix()
      })
    ) %>% 
    select(-data)
  
  return(binned_matrices)
})

transform_count_per_taxon <- map(count_per_taxon, \(tb) {
  transforms <- tb %>% 
    transmute(
      RPK = map2(MAT, GENE_LENGTH, \(x, y) trans_RPK(x, y)),
      TPM = map(RPK, \(x) {
        tpm <- trans_TPM(x)
        tpm[is.na(tpm)] <- 0
        
        tpm
      }),
      rCLR = map(RPK, \(x) {
        clog <- log(x)
        clog[is.infinite(clog)] <- NA
        mean_clog <- colMeans(clog, na.rm = TRUE)
        mean_clog[is.na(mean_clog)] <- 0
        xx <- log(x) - mean_clog
        xx[is.infinite(xx)] <- 0
        xx
      })
    )
  
  gc()
  
  list(
    "TPM" = do.call(rbind, transforms$TPM),
    "rCLR" = do.call(rbind, transforms$rCLR)
  )
}) %>% 
  list_flatten()

save(transform_count_per_taxon, file = "results/transform_count_per_taxon.RData")

