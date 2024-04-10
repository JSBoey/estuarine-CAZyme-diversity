# Normalise count data

# OUTPUTS FROM THIS SCRIPT SHOULD NOT BE USED FOR DIFFERENTIAL EXPRESSION 
# ANALYSES!

# All functions here require a count numeric matrix.

# Different analyses require different types of normalisation. Here, there are 2
# main types:
# (a) Within-sample: 
#     TPM is a popular choice for visualisation due to its ease of 
#     interpretation. It normalises counts by dividing counts by gene length 
#     (kbp) then by library size. It can highlight extremities in data within
#     sample to help with data filtering and/or visualisations. However, it does
#     not account for compositional biases.
# (b) Between-sample:
#     Two popular choices are TMM (edgeR) and RLE (DESeq2, also in edgeR). They 
#     account for sequencing depth and compositional biases. TMM requires one 
#     of the samples to act as reference, whereas RLE uses the ratio between 
#     counts and gene-wise geometric mean as a pseudo-reference. For RLE, there 
#     is a positive count only modification where the product of positive counts 
#     are averaged over all samples.
# (c) Combination:
#     Smid et al. (2018) introduced geTMM, where counts are
#     first transformed into RPK (reads per kilobase length) then normalised by
#     TMM. Based on the paper, it does seem that gene length correction might be
#     also valid for RLE. 
# (d) Per-taxon count:
#     This is a modification of (b) where genes are split into per-taxon tables
#     to be normalised individually where the library size is the sample sum per
#     table. See Klingenberg & Meinicke (2017) for implementation.
# (e) DNA-gene wise paired normalisation:
#     This is the RNA/DNA ratio when there is paired DNA gene-wise counts from 
#     the same sample. Zeroes will be a problem here and will require strict 
#     filtering. THIS IS ONLY VALID IF THE RNA IS CO-ELUTED WITH DNA! See Zhang
#     et al. (2021) for details.

# Notes:
# (a) For the purposes of network reconstruction using Pearson correlations, 
# between-sample normalisation is more important (Johnson & Krishnan, 2022).

# References:
# Robinson M & Oshlack A (2010) A scaling normalization method for differential
#   expression analysis of RNA-seq data doi:10.1186/gb-2010-11-3-r25
#
# Anders S & Huber W (2010) Differential expression analysis for sequence count 
#   data doi:10.1186/gb-2010-11-10-r106
# 
# Smid M, ven der Braak RRJC, van de Werken HJG..., Sieuwerts AM (2018) Gene 
#   length corrected trimmed mean of M-values (GeTMM) processing of RNA-seq data 
#   performs similarly in intersample analyses while improving intrasample
#   comparisons doi:10.1186/s12859-018-2246-7
# 
# Klingenberg H & Meinicke P (2021) How to normalize metatranscriptomic count 
#   data for differential expression analysis doi:10.7717/peerj.3859
# 
# Zhang Y, Thompson KN, Huttenhower C, Franzosa EA (2021) Statistical approaches 
#   for differential expression analysis in metatranscriptomics 
#   doi:10.1093/bioinformatics/btab327

# Environment ----
library(edgeR)
library(DESeq2)
library(tidyverse)

# Import data ----
# wgs <- read_tsv("results/WGS_clean_count.tsv.gz")
# wts <- read_tsv("results/WTS_clean_count.tsv.gz")
# wgs_summary <- read_tsv("results/WGS_count.tsv.summary")
# wts_summary <- read_tsv("results/WTS_count.tsv.summary")
# ann <- read_tsv("results/allbins_pred.annotation_table.tsv.gz")

# Clean data ----
rename_sample <- function(x) str_replace(x, ".*((Filt|Sed).S\\d_\\d).*", "\\1")
# 
# wgs_summary <- rename_with(wgs_summary, rename_sample)
# wts_summary <- rename_with(wts_summary, rename_sample)

# Preprocessing ----
## Create count matrices
# data <- list(
#   "wts_count" = wts,
#   "wgs_count" = wgs
# ) %>% 
#   map(\(x) {
#     select(x, Geneid, contains(c("Filt", "Sed"))) %>% 
#       column_to_rownames("Geneid") %>% 
#       as.matrix()
#   })

# Subset samples by those that were sequenced for the transcriptome
# wts_samples <- colnames(data$wts_count)
# data$wgs_count <- data$wgs_count[, colnames(data$wgs_count) %in% wts_samples]

# Gene length data
# data$gene_length <- ann %>% 
#   select(node, start, end) %>% 
#   mutate(gene_length = end - start + 1) %>% 
#   column_to_rownames("node") %>% 
#   as.matrix() %>% 
#   .[, "gene_length"]

# Mapping library sizes
# data$wts_lib_size <- colSums(select(wts_summary, contains(c("Filt", "Sed"))))[wts_samples]
# data$wgs_lib_size <- colSums(select(wgs_summary, contains(c("Filt", "Sed"))))[wts_samples]

# TPM ----
norm_TPM <- function(count_matrix, gene_length) {
  
  rpk <- sweep(count_matrix, 1, gene_length, "/")
  tpm <- sweep(rpk, 2, colSums(rpk)/1e6, "/")
  
  tpm
  
}

# TMM ----
# This will include the following method:
# 1. geTMM

norm_TMM <- function(count_matrix, method = "TMM", norm_gene_length = FALSE, 
                     ref_column = NULL, gene_length = NULL, logratioTrim = 0.3, 
                     Acutoff = -1e10, sumTrim = 0.05, out_type = "cpm", 
                     lib_size = NULL) {
  
  METHODS <- c("TMM", "TMMwsp")
  
  # Check if methods match
  if (!method %in% METHODS) {
    stop("Argument 'method' must be one of 'TMM' or 'TMMwsp'")
  }
  
  if (isTRUE(norm_gene_length) & is.null(gene_length)) {
    stop("Gene length must be provided if using normalising with gene length")
  }
  
  # Check if output type match
  OUT_TYPES <- c("cpm", "ctf")
  
  if (!out_type %in% OUT_TYPES) {
    stop("Argument 'out_type' must be one of 'cpm' or 'ctf'")
  }
  
  # If geTMM, calculate RPK
  if (isTRUE(norm_gene_length)) {
    count_matrix <- sweep(count_matrix, 1, gene_length, "/")
  }
  
  # Scaling factors
  scale_factors <- normLibSizes(count_matrix, method = method, 
                                refColumn = ref_column, sumTrim = sumTrim,
                                logratioTrim = logratioTrim, Acutoff = Acutoff)
  
  # Normalised library sizes
  if (is.null(lib_size)) {
    lib_size <- colSums(count_matrix)
  }
  
  normalised_library <- scale_factors * lib_size
  
  # Normalised counts
  norm_count <- switch(
    out_type,
    cpm = sweep(count_matrix, 2, normalised_library/1e6, "/"),
    ctf = sweep(count_matrix, 2, scale_factors, "/")
  )
  
}

# RLE ----
# This version of RLE only takes positive counts as in estimateSizeFactors 
# with type = "poscounts"
# This will also include the following method:
# 1. Gene length correction (as per Smid et al. (2018), where RLE scaling was
#    performed, then gene length were taken into account via a TPM calculation)

norm_RLE <- function(count_matrix, gene_length = NULL, 
                     norm_gene_length = FALSE) {
  
  if (isTRUE(norm_gene_length) & is.null(gene_length)) {
    stop("Gene length must be provided if using normalising with gene length")
  }
  
  # RPK if gene length normalisation is required
  if (isTRUE(norm_gene_length)) {
    count_matrix <- sweep(count_matrix, 1, gene_length, "/")
  }
  
  # Only use positive counts for geometric means
  size_factor <- estimateSizeFactorsForMatrix(count_matrix, type = "poscounts")
  
  # Normalise counts
  norm_count <- sweep(count_matrix, 2, size_factor, "/")
  
}

# Per-taxon ----
# This function only splits the data into individual bins and then maps any of 
# the other functions here to it.
# Béchade et al. (2023) used this with TPM on a per microbial symbiont basis.
norm_PerTaxon <- function(count_matrix, features, taxa, fun, ...) {
  
  # Split features by taxa
  features_per_taxon <- split(features, taxa)
  
  # Counts per taxon
  counts_per_taxon <- map(
    features_per_taxon, 
    \(x) count_matrix[rownames(count_matrix) %in% x, ]
  )
  
  # Apply function
  norm_counts_per_taxon <- map(
    counts_per_taxon,
    \(x) fun(x, ...)
  )
  
  # Recombine data
  do.call(rbind, norm_counts_per_taxon)
  
}

# Béchade B, Cabuslay C, Hu Y, Mendonca CM, Hassanpour B, Lin JY, et al. (2023) 
# Physiological and evolutionary contexts of a new symbiotic species from the 
# nitrogen-recycling gut community of turtle ants. ISME 17:1751-1764