# Set up for parallelised correlation: TPM normalisation

# Environment ----
library(edgeR)
library(DESeq2)
library(tidyverse)

source("scripts/in_progress/normalise_counts.R")

# Import data ----
wgs <- read_tsv("results/WGS_clean_count.tsv.gz")
wts <- read_tsv("results/WTS_clean_count.tsv.gz")

# Clean data ----
gene_length <- wts$End - wts$Start + 1
names(gene_length) <- wts$Geneid

# Split tables ----
count_data <- list(
  "water" = select(wts, Geneid, starts_with("Filt")),
  "sediment" = select(wts, Geneid, starts_with("Sed"))
) %>% 
  map(\(x) as.matrix(column_to_rownames(x, "Geneid")))

# Create sediment sum table ----
sediment_patterns <- colnames(count_data$sediment) %>% 
  str_remove("_\\d+$") %>% 
  unique()

count_data$sediment_sum <- matrix(
  nrow = nrow(count_data$sediment), 
  ncol = length(sediment_patterns), 
  dimnames = list(rownames(count_data$sediment), sediment_patterns)
)

for (i in seq_along(sediment_patterns)) {
  
  pattern <- sediment_patterns[i]
  count_data$sediment_sum[, pattern] <- rowSums(
    count_data$sediment[, grep(pattern, colnames(count_data$sediment))]
  )
  
  rm(i, pattern)
  
}

# Normalise data ----
tpm_data <- map(count_data, \(x) norm_TPM(x, gene_length = gene_length))

# Data 1: Remove non-expressed ----
tpm_data <- map(tpm_data, \(x) x[rowSums(x) > 0, ])

# Data 2: Expressed in >= a 1/3 of samples ----
one_third_data <- map(tpm_data, \(x) x[rowSums(x > 0) >= ncol(x)/3, ])

## Save one_third_data as separate files for parallelised correlations
map2(one_third_data, names(one_third_data), \(x, nx) {
  filename <- paste0("data/", nx, ".tpm.one_third.tsv")
  data.table::fwrite(as.data.frame(x), filename, row.names = T)
})

# Data 3: Singletons ----
singletons <- map(tpm_data, \(x) x[rowSums(x > 0) == 1, ])

