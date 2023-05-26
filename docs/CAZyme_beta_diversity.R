# CAZyme beta diversity

# Environment ----
library(vegan)
library(edgeR)
library(compositions)
library(ROptSpace)

library(tidyverse)

source("scripts/utility_functions.R")

# Data ----
env_data <- read_tsv("data/sample_metadata.txt")

dbcan <- read_tsv("results/clean.dbcan.tsv.gz")

callback <- DataFrameCallback$new(cazyFilter)

count_data <- list(
  "wgs" = "results/WGS_clean_count.tsv.gz",
  "wts" = "results/WTS_clean_count.tsv.gz"
) %>% 
  map(read_tsv_chunked, callback = callback)

count_summary <- list(
  "wgs" = "results/WGS_count.tsv.summary",
  "wts" = "results/WTS_count.tsv.summary"
) %>% 
  map(read_tsv) %>% 
  map(rename_with, fcNameClean)

# Normalise count data ----
norm_method <- c("tmm", "tpm", "tmm_tpm", "rclr", "trclr") %>% 
  set_names(str_to_upper(.))

norm_count <- map2(count_data, count_summary, \(df, df_sum) {
  map(norm_method, \(method) {
    ncdf <- normaliseCounts(df, df_sum, method, "all")
    zncdf <- subset(ncdf, subset = rowSums(abs(ncdf)) > 0)
  })
})

# Matrix completion for rCLR transformed data ----
norm_count <- map(norm_count, \(l) {
  complete_rclr <- map(l[c("RCLR", "TRCLR")], completeMatrix)
  names(complete_rclr) <- paste0(names(complete_rclr), "_COMP")
  append(l, complete_rclr)
})

# Calculate distances ----
## For TPM, TMM, TMM_TPM, use Bray-Curtis
## For rCLR and trCLR, use Euclidean
D <- map(norm_count, \(x) {
  map2(x, names(x), \(m, nm) {
    if (nm %in% c("TMM", "TPM", "TMM_TPM")) {
      vegdist(t(m), method = "bray")
    } else {
      vegdist(t(m), method = "euclidean")
    }
  })
})

# Compare distances ----
dist_compare <- map(D, \(d) {
  method <- names(d)
  combination <- as.data.frame(t(combn(method, 2)))
  combination %>% 
    mutate(
      test = map2(V1, V2, ~ cor.test(d[[.x]], d[[.y]])),
      rho = map_dbl(test, "estimate"),
      p.value = map_dbl(test, "p.value")
    ) %>% 
    select(-test)
})

ggplot(dist_compare$wgs, aes(x = V1, y = V2, fill = rho)) +
  geom_tile()

# Ordinations ----


# Compare ordinations ----


# Plot ----


