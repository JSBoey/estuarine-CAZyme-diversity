# SparCC network

# Environment ----
library(NetCoMi)
library(vegan)
library(edgeR)
library(igraph)
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
  # Filter criteria
  fnode <- dbcan$X1[grepl("GH|PL", dbcan$X3)]
  data <- subset(y, subset = rownames(y) %in% fnode)
  
  # Remove nodes with less than 30% incidence
  imat <- ifelse(data > 0, 1, 0)
  fnode2 <- names(which(rowSums(imat) >= 0.3 * ncol(imat)))
  subset(data, subset = rownames(data) %in% fnode2)
})

# SPIEC-EASI ----
spieceasi_water <- spiec.easi(
  t(sub_data$WTS$Water),
  method = 'mb',
  lambda.min.ratio = 50e-3,
  nlambda = 20,
  pulsar.params = list(
    rep.num = 50
  )
)

spieceasi_sediment <- spiec.easi(
  t(sub_data$WTS$Sediment),
  method = 'mb',
  lambda.min.ratio = 50e-3,
  nlambda = 20,
  pulsar.params = list(
    rep.num = 50
  )
)

igMB <- list(
  "Water" = adj2igraph(getRefit(spieceasi_water)),
  "Sediment" = adj2igraph(getRefit(spieceasi_sediment))
)

igCoord <- map(igMB, layout.fruchterman.reingold)

pdf("results/CAZyme_SpiecEasi.pdf", width = 8, height = 4)
par(mfrow = c(1, 2))

pmap(list(igMB, igCoord, names(igCoord)), ~ {
  main <- paste("SPIEC-EASI (MB) network of", ..3)
  plot(..1, layout = ..2, vertex.size = 3, vertex.label = NA, main = main)
})

dev.off()

save.image("results/CAZyme_sparse_networks.RData")

# WGCNA ----
library(WGCNA)

powers <- c(c(1:20), seq(22, 40, by = 2))
sft <- map(sub_data$WTS, ~ {
  pickSoftThreshold(t(.x), powerVector = powers, verbose = 5)
})

par(mfrow = c(2, 2))
map2(sft, names(sft), function(sft, nm) {
  # Plot scale-free topology fit index as a function of soft-thresholding power
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
       main = paste("Scale independence for", nm))
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=0.9,col="red")
  # Mean connectivity as a function of soft-thresholding power
  plot(sft$fitIndices[,1], sft$fitIndices[,5],
       xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
       main = paste("Mean connectivity for", nm))
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=0.9,col="red")
})

net <- map2(sub_data$WTS, names(sub_data$WTS), ~ {
  blockwiseModules(t(.x), power = 7,
                   TOMType = "unsigned", minModuleSize = 30,
                   reassignThreshold = 0, mergeCutHeight = 0.25,
                   numericLabels = TRUE, pamRespectsDendro = FALSE,
                   saveTOMs = TRUE,
                   saveTOMFileBase = paste0("CAZymeTOM", .y),
                   verbose = 3)
})




