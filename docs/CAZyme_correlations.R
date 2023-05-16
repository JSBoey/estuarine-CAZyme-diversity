# CAZyme correlations

# Environment ----
library(edgeR)
library(WGCNA)
library(vegan)
library(tidyverse)

# Utilities ----
# Feature count summaries tend to retain the file name as column names. This removes them.
fcNameClean <- function(x) {
  str_remove(x, ".*\\/") %>% 
    str_remove("W.S\\.") %>% 
    str_remove(".bam")
}

# Convert data into a DGEList
data2DGE <- function(count, count_summary) {
  # Get library size from count_summary
  lib.size <- colSums(count_summary[, 2:ncol(count_summary)])/2
  # Convert count to DGEList
  data <- column_to_rownames(count, "Geneid")
  DGEList(
    counts = data, 
    annotation.columns = c("Chr", "Start", "End", "Strand", "Length"),
    lib.size = lib.size
  ) %>% 
    calcNormFactors()
}

# Convert DGEList to TMM-normalised TPM
dge2tpm <- function(dge) {
  RPKM <- rpkm(dge)
  TPM <- t( t(RPKM) / colSums(RPKM) ) * 1e6
}

# Data ----
dbcan <- read_tsv("results/allbins_pred.dbcan_parsed.tsv.gz", col_names = F)
wts_count <- read_tsv("results/WTS_clean_count.tsv.gz")
wts_summary <- read_tsv("results/WTS_count.tsv.summary") |>  
  rename_with(.fn = fcNameClean, .cols = -Status)
wgs_count <- read_tsv("results/WGS_clean_count.tsv.gz")
wgs_summary <- read_tsv("results/WGS_count.tsv.summary") |> 
  rename_with(.fn = fcNameClean, .cols = -Status)

# Create DGEList ----
dge <- list(
  "WTS" = data2DGE(wts_count, wts_summary),
  "WGS" = data2DGE(wgs_count, wgs_summary)
)

# Obtain normalised counts ----
norm_counts <- map(dge, dge2tpm)

# Split data ----
# Water column and sediment are fundamentally different environments, they should not be put together
data <- map2(norm_counts, dge, ~ {
  # Get gene information
  gene_info <- .y$genes
  # Split count information
  water <- subset(.x, select = grepl("Filt", colnames(.x)))
  sediment <- subset(.x, select = grepl("Sed", colnames(.x)))
  # Return list
  list(
    "Water" = water,
    "Sediment" = sediment
  )
})

# Subset data ----
sub_data <- map_depth(data, 2, \(y) {
  subset(y, subset = rownames(y) %in% dbcan$X1)
})

# Correlate CAZyme gene and transcript count ----
# Question: Does transcript follow gene distribution?

corGeneTranscript <- function(gene, transcript, rclr = TRUE) {
  # Normalise data
  if (isTRUE(rclr)) {
    gene <- t(decostand(t(gene), method = "rclr"))
    transcript <- t(decostand(t(transcript), method = "rclr"))
    cat("Normalised matrix\n")
  }
  # Make sure the data is equivalent
  gene <- subset(gene, select = colnames(gene) %in% colnames(transcript))
  
  # Find genes with non-zero sums
  nzg <- names(which(rowSums(gene) > 0))
  nzt <- names(which(rowSums(transcript) > 0))
  
  # Subset matrices to shared row names
  sn <- base::intersect(nzg, nzt)
  sg <- subset(gene, subset = rownames(gene) %in% sn)
  st <- subset(transcript, subset = rownames(transcript) %in% sn)
  
  cat(paste("Found", length(sn), "non-zero genes and transcripts\n"))
  
  # Calculate correlation for same nodes
  nr <- length(sn)
  coef <- vector(length = nr)
  pval <- vector(length = nr)
  zg <- vector(length = nr)
  zt <- vector(length = nr)
  
  cat("Calculating correlations...\n")
  
  for (i in 1:nr) {
    r <- cor.test(sg[i, ], st[i, ], method = "spearman")
    coef[i] <- r$estimate
    pval[i] <- r$p.value
    zg[i] <- sum(sg[i, ] == 0)/length(sg[i, ])
    zt[i] <- sum(st[i, ] == 0)/length(st[i, ])
  }
  
  cat("Complete\n\n")
  
  # Adjust p-value
  adjp <- p.adjust(pval, method = "fdr")
  
  # Return data frame
  data.frame(
    "Geneid" = sn,
    "rho" = coef,
    "p" = pval,
    "adj.p" = adjp,
    "zero.gene" = round(zg, 3),
    "zero.transcript" = round(zt, 3)
  )
}

gt_cor <- map(
  transpose(sub_data), 
  ~ corGeneTranscript(gene = .x$WGS, transcript = .x$WTS)
)

## Plot correlations ----
par(mfrow = c(2, 2))

## Correlations
map2(gt_cor, names(gt_cor), ~ {
  x <- .x$rho[order(.x$adj.p)]
  y <- .x$adj.p[order(.x$adj.p)]
  col <- ifelse(y < 0.05, "red", "grey")
  main <- .y
  
  plot(x = x, y = y, col = col, pch = 20, 
       main = main, xlab = "Spearman's rho", ylab = "FDR")
})

## Zero structures
map2(gt_cor, names(gt_cor), ~ {
  x <- .x$rho[order(.x$zero.transcript)]
  y <- .x$zero.transcript[order(.x$zero.transcript)]
  col <- ifelse(y > 0.3, "grey", "green")
  main <- paste("Missingness in", .y)
  
  plot(x = x, y = y, col = col, pch = 20, 
       main = main, xlab = "Spearman's rho", ylab = "Transcript")
})

## Subset significant correlations ----
sig_gt_cor <- map(gt_cor, ~ {
  y <- select(dbcan, X1, X3)
  
  filter(.x, adj.p < 0.05) |> 
    left_join(select(dbcan, X1, X3), by = c("Geneid" = "X1"))
})
