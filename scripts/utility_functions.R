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
  if (!requireNamespace("tidyverse", quietly = TRUE)) library(tidyverse)
  if (!requireNamespace("edgeR", quietly = TRUE)) library(edgeR)
  if (!requireNamespace("vegan", quietly = TRUE)) library(vegan)
  
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

# Average frequency of a gene across a matrix of samples ----
aveGeneFreq <- function(count_matrix) {
  # Per sample frequency table 
  freq_table <- proportions(count_matrix, 2)
  # Number of samples
  N <- ncol(freq_table)
  # Average gene frequencies across samples
  p_i <- rowSums(freq_table) / N
  
  return(p_i)
}

# Specificity of a gene for a sample ----
geneSpecificity <- function(count_matrix, log.base = exp(1)) {
  # Per sample frequency table
  freq_table <- proportions(count_matrix, 2)
  # Number of samples
  N <- ncol(freq_table)
  # Per gene average frequency across samples
  p_i <- aveGeneFreq(count_matrix = count_matrix)
  # Calculate frequency ratios
  r <- sweep(freq_table, 1, p_i, "/")
  # Gene specificity
  S_i <- rowSums(r * log(r, base = log.base), na.rm = TRUE) / N
  
  return(S_i)
}

# Specificity of the transcript distribution of a sample ----
sampleSpecificity <- function(count_matrix, log.base = exp(1)) {
  # Per sample frequency table
  freq_table <- proportions(count_matrix, 2)
  # Gene specificity
  S_i <- geneSpecificity(count_matrix = count_matrix, log.base = log.base)
  # Average gene specificity
  aS_i <- sweep(freq_table, 1, S_i, "*")
  # Sample specificity
  d_j <- colSums(aS_i, na.rm = FALSE)
  
  return(d_j)
}

# Average transcript log frequencies of a sample ----
aveSampleFreq <- function(count_matrix, log.base = exp(1)) {
  # Per sample frequency table
  freq_table <- proportions(count_matrix, 2)
  # Log-scaled per gene average frequency across samples
  lp_i <- log(aveGeneFreq(count_matrix = count_matrix), base = log.base)
  # Per sample average transcript log frequencies
  H_Rj <- - colSums(sweep(freq_table, 1, lp_i, "*"), na.rm = T)
  
  return(H_Rj)
}

# Kullback-Leiber divergence of a sample ---- 
KLDivSample <- function(count_matrix, log.base = exp(1)) {
  H_Rj <- aveSampleFreq(count_matrix = count_matrix, log.base = log.base)
  H_j <- diversity(count_matrix, MARGIN = 2, index = "shannon")
  
  D_j <- H_Rj - H_j
  
  return(D_j)
}
