# Utility scripts

# Cleaning sample names from feature count outputs ----
fcNameClean <- function(x) {
  if (!requireNamespace("tidyverse", quietly = TRUE)) library(tidyverse)
  
  str_remove(x, ".*\\/") %>% 
    str_remove("W.S\\.") %>% 
    str_remove(".bam")
}

# Filter function for callback to subset data to CAZyme genes only ----
cazyFilter <- \(df, pos) {
  filter(df, Geneid %in% dbcan$query)
}

# Extract only count data ---- 
makeCountMatrix <- function(count_data) {
  if (!requireNamespace("tidyverse", quietly = TRUE)) library(tidyverse)
  
  count_data %>% 
    select(Geneid, contains(c("Filt", "Sed"))) %>% 
    column_to_rownames("Geneid") %>% 
    as.matrix()
}

# Normalize counts ----
# Methods are:
# - TPM (tpm)
# - TMM (tmm)
# - TPM with TMM normalised library size (tmm_tpm)
# - robust centered log ratio (rclr)
# - rCLR after gene length normalisation (trclr)
# Library size options:
# - All WTS reads (all): Obtain from mapping summary
# - Based on column sums of count data (csum)
# - Based on all mapped reads (mapped): Obtain from mapping summary
normaliseCounts <- function(count_data, mapping_summary, method, use.lib.size) {
  if (!requireNamespace("tidyverse", quietly = TRUE)) {library(tidyverse)}
  if (!requireNamespace("edgeR", quietly = TRUE)) {library(edgeR)}
  if (!requireNamespace("vegan", quietly = TRUE)) {library(vegan)}
  
  # Turn mapping summary into a matrix
  mapping_summary <- column_to_rownames(mapping_summary, "Status") %>% 
    as.matrix()
  
  # Library size
  if (use.lib.size == "all") {
    lib.size <- colSums(mapping_summary)
  } else if (use.lib.size == "mapped") {
    lib.size <- mapping_summary[1, ]
  } else if (use.lib.size == "csum") {
    lib.size <- colSums(makeCountMatrix(count_data))
  }
  
  # Convert to DGE list
  dge <- DGEList(
    counts = column_to_rownames(count_data, "Geneid"), 
    annotation.columns = c("Chr", "Start", "End", "Strand", "Length"),
    lib.size = lib.size
  )
  
  # If require TMM library size normalization
  if (method %in% c("tmm", "tmm_tpm")) {
    dge <- calcNormFactors(dge, method = "TMM")
  }
  
  # Return TPM if requested
  if (method %in% c("tpm", "tmm_tpm")) {
    RPKM <- rpkm(dge)
    norm_matrix <- t( t(RPKM) / colSums(RPKM) ) * 1e6
  }
  
  # Return TMM-normalised counts-per-million if requested
  if (method == "tmm") {
    norm_matrix <- cpm(dge, normalized.lib.sizes = TRUE, log = FALSE)
  }
  
  # If require gene length-normalisation (per kb) (i.e. normalised transcript)
  if (method == "trclr") {
    dge[["counts"]] <- sweep(getCounts(dge), 1, (dge[["genes"]]$Length / 1000), "/")
  }
  
  # Return robust CLR if requested
  if (method %in% c("rclr", "trclr")) {
    norm_matrix <- decostand(getCounts(dge), method = "rclr", MARGIN = 2)
  }
  
  return(norm_matrix)
  
}

# Matrix completion with OptSpace ----
completeMatrix <- function(m) {
  # Set zeroes to NA
  m[m == 0] <- NA
  # Run OptSpace
  op <- OptSpace(m)
  # Get completed matrix
  M <- with(op, X %*% S %*% t(Y))
  # Fill matrix
  m[is.na(m) == TRUE] <- M[is.na(m) == TRUE]
  
  return(m)
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


