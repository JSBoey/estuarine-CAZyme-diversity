# Utility scripts

# Cleaning sample names from feature count outputs ----
fcNameClean <- function(x) {
  if (!requireNamespace("tidyverse", quietly = TRUE)) library(tidyverse)
  
  str_remove(x, ".*\\/") %>% 
    str_remove("W.S\\.") %>% 
    str_remove(".bam")
}

# Normalise counts to TPM based on TMM-normalised library sizes ----
count2TPM <- function(count, count_summary) {
  if (!requireNamespace("tidyverse", quietly = TRUE)) library(tidyverse)
  if (!requireNamespace("edgeR", quietly = TRUE)) library(edgeR)
  
  # Library size
  lib.size <- colSums(count_summary[, 2:ncol(count_summary)])/2
  
  # Convert count to DGEList
  data <- column_to_rownames(count, "Geneid")
  
  dge <- DGEList(
    counts = data, 
    annotation.columns = c("Chr", "Start", "End", "Strand", "Length"),
    lib.size = lib.size
  ) %>% 
    calcNormFactors()
  
  # Convert to TPM
  RPKM <- rpkm(dge)
  TPM <- t( t(RPKM) / colSums(RPKM) ) * 1e6
  
  # Clean up for tidyverse format
  as.data.frame(TPM) |> 
    rownames_to_column("Geneid")
}

# Count number ORFs with multiple annotation hits ----
countMultipleHits <- function(x, query_pattern) {
  if (!requireNamespace("tidyverse", quietly = TRUE)) library(tidyverse)
  
  # Column with ORF
  qcol <- x %>% 
    map_lgl(~ all(str_detect(.x, query_pattern))) %>% 
    which() %>% 
    names()
  # Find number of duplicates
  sum(table(x[[qcol]]) > 1)
}



