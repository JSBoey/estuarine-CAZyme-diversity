# CAZyme beta-diversity

# Preamble ----
# After some exploratory analysis, I've come to the conclusion that beta-diversity based on TPM and robust CLR is likely most helpful. Both account for library size variation in different ways: TPM re-casts counts on a per-million transcript scale, while robust CLR accounts for the compositional nature of the data. TPM also adjusts for read recruitment rates affected by gene length, whereas robust CLR does not. However, this may not be as important given that we are comparing across samples, and are not primarily interested in comparisons within samples.

# TPM:
# Ordinations: nMDS of Bray-Curtis dissimilarity. 
# Effects of constraints: dbRDA paired with tests of beta-dispersion based on sample type.

# Robust CLR:
# Ordinations: PCA of Euclidean distances.
# Effects of constraints: RDA with permutation.

# The overall beta-diversity tests are then followed up with comparisons between data types (WGS v. WTS) based on a complete pairwise samples. WTS ordinations will be fit into the ordination space of WGS.

# Environment ----
library(vegan)
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

# Normalise count ----
norm_method <- c("tpm", "rclr") %>% 
  set_names(str_to_upper(.))

norm_count <- map2(count_data, count_summary, \(df, df_sum) {
  map(norm_method, \(method) {
    ncdf <- normaliseCounts(df, df_sum, method, "all")
    zncdf <- subset(ncdf, subset = rowSums(abs(ncdf)) > 0)
  })
})

# Analysis 1: Beta-diversity ----
## Ordination ----
O <- map(norm_count, \(l) {
  map2(l, names(l), \(mat, nm) {
    if (nm == "TPM") {
      metaMDS(
        t(mat),
        distance = "bray",
        k = 2,
        try = 999,
        trymax = 999,
        autotransform = FALSE
      )
    } else {
      rda(t(mat))
    }
  })
})

## Constrained ordination ----
CO <- map(norm_count, \(l) {
  map2(l, names(l), \(mat, nm) {
    # Subset relevant rows in environmental data
    data <- filter(env_data, sample %in% colnames(mat))
    
    # Constrained ordinations
    if (nm == "TPM") {
      dbrda(
        t(mat) ~ type + salinity + type:salinity, data, 
        distance = "bray", 
        sqrt.dist = FALSE, 
        metaMDSdist = FALSE
      )
    } else {
      rda(t(mat) ~ type + salinity + type:salinity, data)
    }
  })
})

CO_summary <- map_depth(CO, 2, summary)

## Test of significance: Constrained ordination ----
TCO <- map_depth(CO, 2, \(x) anova.cca(x, by = "margin"))

## Beta-dispersion (variance) by sample type ----
BD <- map(norm_count, \(l) {
  map2(l, names(l), \(mat, nm) {
    # Subset relevant rows in environmental data
    data <- filter(env_data, sample %in% colnames(mat))
    
    # Calculate dissimilarities/distances
    if (nm == "TPM") {
      D <- vegdist(t(mat), method = "bray")
    } else {
      D <- vegdist(t(mat), method = "euclidean")
    }
    
    # Beta-dispersion
    betadisper(D, data$type, bias.adjust = TRUE)
  })
})

## Test of significance: Beta-dispersion ----
TBD <- map_depth(BD, 2, permutest)

# Plot Analysis 1 ----
## Constrained ordinations ----
### Extract scores

### Join tables

### Ordination plot

## Beta-dispersion ----
### Extract distance-to-centroids

### Boxplot


# Analysis 2: WGS v. WTS ----
## Subset relevent data ----
tNC <- transpose(norm_count) %>% 
  map(\(x) {
    # Subset data to have same samples
    x[["wgs"]] <- subset(
      x[["wgs"]], 
      select = colnames(x[["wgs"]]) %in% colnames(x[["wts"]])
    )
    
    return(x)
  })

## Generate ordinations
OPr <- map2(tNC, names(tNC), \(l, nm) {
  # Ordination
  if (nm == "TPM") {
    O <- map(l, \(mat) {
      metaMDS(
        t(mat),
        k = 2, 
        distance = "bray",
        autotransform = FALSE,
        try = 999,
        trymax = 999
      )
    })
  } else {
    O <- map(l, \(mat) {
      rda(t(mat))
    })
  }
})

## Procrustes projections ----
Pr <- map(OPr, \(ord) {
  procrustes(ord[["wgs"]], ord[["wts"]], symmetric = TRUE)
})

## Test of significance ----
TPr <- map(OPr, \(ord) {
  protest(ord[["wgs"]], ord[["wts"]])
})

## Plot Analysis 2 ----
### Extract points

### Extract rotaion plane

### Assemble tables

### Plot projections


