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
  # Transformation and normalisation methods
  method <- names(d)
  
  # Combinations
  combination <- as.data.frame(t(combn(method, 2)))
  
  # Pairwise correlations
  combination %>% 
    mutate(
      test = map2(V1, V2, ~ cor.test(d[[.x]], d[[.y]])),
      rho = map_dbl(test, "estimate"),
      p.value = map_dbl(test, "p.value")
    ) %>% 
    select(-test)
  
})

# Ordinations ----
O <- map(norm_count, \(nc) {
  
  map2(nc, names(nc), \(normalised_count, method) {
    
    # Different ordination method depending on data type
    comp_count <- str_detect(method, "CLR")
    comm <- t(normalised_count)
    
    if (isFALSE(comp_count)) {
      
      metaMDS(
        comm, distance = "bray", try = 999, trymax = 999, autotransform = FALSE
      )
      
    } else {
      
      rda(comm)
      
    }
    
  })
  
})

# Compare ordinations ----
procrustes_self <- map(O, \(ordination) {
  
  # Transformation and normalisation methods
  method <- names(ordination)
  
  # Combinations
  combination <- as.data.frame(t(combn(method, 2)))
  
  # Pairwise Procrustes test
  combination %>%
    mutate(
      test = map2(V1, V2, ~ protest(ordination[[.x]], ordination[[.y]])),
      ss = map_dbl(test, "ss"),
      t0 = map_dbl(test, "t0"),
      signif = map_dbl(test, "signif"),
      across(c(ss, t0, signif), \(x) round(x, 3))
    )
  
}) %>% 
  bind_rows(.id = "data_type")

pdf(
  "results/Normalisation_transformation_comparisons_ordinations.pdf", 
  width = 8.3, height = 11.7
)
par(mfrow = c(3, 2), pty = "s")
for (i in 1:nrow(procrustes_self)) {
  plot_data <- procrustes_self$test[[i]]
  main <- paste(
    procrustes_self$data_type[[i]], procrustes_self$V1[[i]], procrustes_self$V2[[i]]
  )
  sub <- paste0(
    "SS: ", procrustes_self$ss[[i]], 
    " Correlation: ", procrustes_self$t0[[i]], 
    " P-value: ", procrustes_self$signif[[i]]
  )
  plot(plot_data, kind = 1, main = main, sub = sub)
  plot(plot_data, kind = 2)
}
dev.off()

# Compare interpretability of methods ----
interpretability <- map(norm_count, \(nc) {
  
  map2(nc, names(nc), \(normalised_count, method) {
    
    # Different tests depending on data type
    comp_count <- str_detect(method, "CLR")
    comm <- t(normalised_count)
    metadata <- filter(env_data, sample %in% rownames(comm))
    
    if (isFALSE(comp_count)) {
      
      dbrda(comm ~ type + salinity, metadata, sqrt.dist = F, distance = "bray")
      
    } else {
      
      rda(comm ~ type + salinity, metadata)
      
    }
    
  })
  
})

interpretability_test <- map_depth(interpretability, 2, \(x) {
  anova.cca(x, by = "margin", parallel = 4)
})

anova_df <- map_depth(interpretability_test, 2, \(x) {
  broom::tidy(x) %>% 
    rename_with(
      ~ str_replace(.x, "Variance|SumOfSqs", "Variation")
    )
}) %>% 
  map(bind_rows, .id = "transform") %>% 
  bind_rows(.id = "data_type") %>% 
  mutate(
    "Variation_type" = if_else(str_detect(transform, "CLR"), "Variance", "Sum of squares")
  )

write_tsv(anova_df, "results/Normalisation_transformation_PerMANOVA.tsv")

# Plot ----


