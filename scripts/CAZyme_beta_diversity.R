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
library(edgeR)
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
TCO <- map_depth(CO, 2, \(x) {
  model <- c("margin", "terms", "axis")
  names(model) <- model
  map(model, \(m) {
    anova.cca(x, by = m, permutations = 99)
  })
})

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
## Unconstrained ordinations ----
### Extract scores and join tables ----
O_score <- map_depth(O, 2, \(ord) {
  scores(ord, scaling = "sites", tidy = TRUE) %>% 
    left_join(env_data, by = c("label" = "sample"))
})

### Ordination plot ----
O_plot <- map_depth(O_score, 2, \(S) {
  # Filter score type to sites
  S <- filter(S, score == "sites")
  # Get axis names before standardisation
  axis_name <- names(S)[1:2]
  # Rename axes
  names(S)[c(1, 2)] <- c("x", "y")
  ggplot(
    data = S, 
    mapping = aes(x = x, y = y)
  ) +
    geom_hline(
      yintercept = 0, colour = "grey40", linetype = 2
    ) +
    geom_vline(
      xintercept = 0, colour = "grey40", linetype = 2
    ) +
    geom_point(
      mapping = aes(colour = salinity, shape = type)
    ) +
    labs(
      x = axis_name[[1]], 
      y = axis_name[[2]], 
      colour = "Salinity ‰", 
      shape = "Sample type"
    ) +
    scale_colour_viridis_c(
      limits = c(0, 35)
    ) +
    scale_shape_manual(
      values = c("sediment" = 15, "water" = 16),
      labels = str_to_title
    ) +
    theme_bw() +
    theme(
      aspect.ratio = 1,
      panel.grid = element_blank()
    )
})

## Constrained ordinations ----
### Extract scores and join tables 
CO_score <- map_depth(CO, 2, \(ord) {
  scores(ord, scaling = "sites", tidy = TRUE) %>% 
    left_join(env_data, by = c("label" = "sample"))
}) %>% 
  list_flatten()

CO_statistic <- list_flatten(TCO)

CO_arrowMultiplier <- map_depth(CO, 2, ~ {
  a <- ordiArrowMul(.x, display = "bp")
  
  if (is.infinite(a)) {
    a <- 1
  } else if (a < 2) {
    a <- 3
  }
  
  return(a)
}) %>% 
  list_flatten()

### Ordination plot
CO_plot <- pmap(
  list(CO_score, CO_statistic, CO_arrowMultiplier), \(scr, st, mul) {
    # Get axis names prior to standardisation
    axis_name <- names(scr)[1:2]
    # Standardise axis names for ggplot call
    names(scr)[1:2] <- c("x", "y")
    # Main plotting variables
    site_score <- filter(scr, score == "sites")
    centroid <- filter(scr, score == "centroids") %>% 
      mutate(
        label = str_to_title(str_remove(label, "type"))
      )
    biplot <- filter(scr, score == "biplot") %>% 
      mutate(
        across(c(x, y), \(axis) axis * mul)
      )
    
    # Axis variables
    axis_variation <- proportions(as.matrix(st[["axis"]][, 2]))[1:2]
    axis_label <- paste0(
      axis_name, " (", round(axis_variation * 100, 2), "%)"
    )
    
    ggplot() +
      geom_hline(
        yintercept = 0, colour = "grey40", linetype = 2
      ) +
      geom_vline(
        xintercept = 0, colour = "grey40", linetype = 2
      ) +
      geom_segment(
        mapping = aes(xend = x, yend = y, x = 0, y = 0),
        data = biplot, 
        arrow = arrow(length = unit(0.03, units = "npc"))
      ) +
      geom_point(
        mapping = aes(x = x, y = y, colour = salinity, shape = type),
        data = site_score
      ) +
      geom_text(
        mapping = aes(x = x, y = y, label = label),
        data = centroid
      ) +
      labs(
        x = axis_label[[1]], 
        y = axis_label[[2]], 
        colour = "Salinity ‰", 
        shape = "Sample type"
      ) +
      scale_colour_viridis_c(
        limits = c(0, 35)
      ) +
      scale_shape_manual(
        values = c("sediment" = 15, "water" = 16),
        labels = str_to_title
      ) +
      theme_bw() +
      theme(
        aspect.ratio = 1,
        panel.grid = element_blank()
      )
  }
)

patchwork::wrap_plots(CO_plot)
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

### Extract plane of rotation

### Assemble tables

### Plot projections


