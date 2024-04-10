# Ordinations

# Preamble ----
# This explores how community-wide prokaryotic carbohydrate degradation (via GH and PL) through their latent (genetic) and realised (transcriptomic) potentials. To what extent are both potentials coupled in each community?
# Note: This script requires output from transform_counts.R.

# Environment ----
library(vegan)
library(usedist)
library(tidyverse)
library(patchwork)

# Import data ----
load("results/transformed_counts.RData")
annotation <- read_tsv("results/allbins_pred.annotation_table.tsv.gz")
metadata <- read_tsv("data/sample_metadata.txt")
wgs_summary <- read_tsv("results/WGS_count.tsv.summary")
wts_summary <- read_tsv("results/WTS_count.tsv.summary")

# Filter data ----
# Retain rows with dbCAN assignments and are non-zero across all samples.

dbcan_nodes <- annotation$node[
  str_detect(annotation$dbcan_label, "GH|PL") & !is.na(annotation$dbcan_label)
]
transform_count <- map(transform_count, \(x) {
  x[rownames(x) %in% dbcan_nodes & rowSums(abs(x)) > 0, ]
})

# Distance matrices ----
D <- map2(
  transform_count, 
  names(transform_count), 
  \(x, nm) {
    d <- ifelse(grepl("TPM", nm), "bray", "euclidean")
    vegdist(t(x), method = d)
  }
)

D_Jac <- map(
  transform_count[c(1, 3)], \(x) {
    vegdist(t(x), method = "jaccard", binary = T)
  }
)

# Ordinations ----
O <- map2(
  D, names(D), 
  \(x, nm) {
    
    if (grepl("TPM", nm)) {
      metaMDS(x, try = 999, trymax = 4999, k = 2)
    } else {
      wcmdscale(x, k = 2)
    }
    
  }
)

O_Jac <- map(D_Jac, metaMDS, try = 999, trymax = 4999, k = 2) 

# Procrustes analyses ----
O_score <- map(O, scores)

protest_nonequal <- function(X, Y) {
  XX <- X[rownames(X) %in% rownames(Y), ]
  
  protest(XX, Y)
}

Pt <- list(
  "TPM" = protest_nonequal(O_score$wgs_TPM, O_score$wts_TPM), 
  "rCLR" = protest_nonequal(O_score$wgs_rCLR, O_score$wts_rCLR)
)

Pt_distance <- map(Pt, \(x) {
  r <- residuals(x)
  list(
    "residuals" = r, 
    "quantile" = quantile(r)
  ) %>% 
    map(as_tibble, rownames = "sample")
})

O_Jac_score <- map(O_Jac, scores)
O_Jac_score$wgs_TPM <- O_Jac_score$wgs_TPM[1:30, ]
Pt_Jac <- protest(O_Jac_score$wgs_TPM, O_Jac_score$wts_TPM)

par(mfrow = c(2, 2))
ordiplot(O_Jac$wgs_TPM, type = "t")
ordiplot(O_Jac$wts_TPM, type = "t")
plot(Pt_Jac, kind = 1)
plot(Pt_Jac, kind = 2)


# Procrustes residuals v salinity and a-diversity ----
Pt_d_TPM <- left_join(Pt_distance$TPM$residuals, 
                      metadata[, c("sample", "salinity", "type")], 
                      by = "sample") %>% 
  mutate(
    lambda_gene = diversity(transform_count$wgs_TPM, index = "invsimpson", MARGIN = 2)[1:30],
    lambda_transcript = diversity(transform_count$wts_TPM, index = "invsimpson", MARGIN = 2),
    lambda_ratio = lambda_transcript/lambda_gene,
    average_gene = matrixStats::colSds(
      ifelse(transform_count$wgs_TPM == 0, NA, transform_count$wgs_TPM),
      na.rm = TRUE
    )[1:30],
    average_transcript = matrixStats::colSds(
      ifelse(transform_count$wts_TPM == 0, NA, transform_count$wts_TPM),
      na.rm = TRUE
    ),
    average_ratio = average_transcript/average_gene,
    total_gene = colSums(transform_count$wgs_TPM)[1:30],
    total_transcript = colSums(transform_count$wts_TPM),
    total_ratio = total_transcript/total_gene
  )

ggplot(Pt_d_TPM, aes(x = salinity, y = value)) + 
  geom_point(aes(colour = total_ratio)) + 
  stat_smooth(method = "loess", se = F) + 
  labs(x = "Salinity", 
       y = "Procrustes residuals", 
       subtitle = "Trend line = LOESS"
       ) + 
  facet_wrap(~ type, ncol = 2, scales = "free_x")

ggplot(Pt_d_TPM, aes(x = lambda_gene, y = lambda_transcript)) + 
  geom_point(aes(colour = salinity)) +
  facet_wrap(~ type, ncol = 2, scales = "free_x") +
  scale_colour_viridis_c()

# Variation ----
V <- map(D, \(d) {
  d_samples <- attr(d, "Labels")
  M <- filter(metadata, sample %in% d_samples)
  adonis2(d ~ type + salinity + ammonia + nitrate + drp + sulfate + dnpoc, 
          data = M, na.action = na.exclude, by = "margin")
})

# Dispersion ----
Dsp <- map(D, \(x) {
  
  M <- filter(metadata, sample %in% attr(x, "Labels")) 
  dist_to_centroids(x, M$type) %>% 
    filter(
      str_detect(Item, "Filt") & CentroidGroup == "water" |
        str_detect(Item, "Sed") & CentroidGroup == "sediment"
    )
  
})

Dsp_test <- map(Dsp, \(x) {
  wilcox.test(CentroidDistance ~ CentroidGroup, x, paired = FALSE)
})

# Plot ----

## Parameters ----
theme_modification <- theme(
  panel.grid = element_blank()
)
  
## Ordinations ----
plot_O <- map2(O_score, names(O_score), \(x, nm) {
  
  x <- as_tibble(x, rownames = "sample") %>%
    left_join(metadata, by = "sample")
  names(x)[2:3] <- c("Axis1", "Axis2")
  
  ggplot(x, aes(x = Axis1, y = Axis2)) +
    geom_hline(yintercept = 0, colour = "grey60") +
    geom_vline(xintercept = 0, colour = "grey60") +
    geom_point(aes(colour = salinity, shape = type), 
               size = 3.5, alpha = 0.8) +
    labs(title = nm, x = "Axis 1", y = "Axis 2",
         colour = "Salinity (‰)", shape = "Sample type") +
    scale_colour_viridis_c(option = "D") +
    scale_shape_discrete(labels = str_to_title) +
    theme_bw() +
    theme_modification +
    theme(
      aspect.ratio = 1
    )

})

wrap_plots(plot_O[grepl('TPM', names(plot_O))], guides = 'collect')

## Procrustes ----
plot_Pt <- map2(Pt, names(Pt), \(x, nm) {
  
  M <- select(metadata, sample, type, salinity)
  
  X <- as_tibble(x$X, rownames = "sample")
  Yrot <- as_tibble(x$Yrot, rownames = "sample")
  
  names(X)[2:3] <- c("from_x", "from_y")
  names(Yrot)[2:3] <- c("to_x", "to_y")
  
  XYrot <- left_join(X, Yrot, by = "sample") %>% 
    left_join(M, by = "sample")
  
  rot <- x$rotation
  
  ggplot(XYrot) +
    geom_hline(yintercept = 0, linetype = 2, colour = "grey80") +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey80") +
    geom_abline(intercept = 0, slope = rot[1, 2]/rot[1, 1], colour = "grey80") +
    geom_abline(intercept = 0, slope = rot[2, 2]/rot[2, 1], colour = "grey80") +
    geom_segment(aes(x = from_x, xend = to_x, y = from_y, yend = to_y), 
                 arrow = arrow(length = unit(2, "mm")),
                 colour = "navy") +
    geom_point(aes(x = from_x, y = from_y, colour = salinity, shape = type),
               size = 3.5) +
    labs(title = nm, x = "Axis 1", y = "Axis 2",
         colour = "Salinity (‰)", shape = "Sample type") +
    scale_colour_viridis_c(option = "D") +
    scale_shape_discrete(labels = str_to_title) +
    theme_bw() +
    theme_modification
})

wrap_plots(plot_Pt, guides = 'collect')

plot_Pt_distance <- map2(Pt_distance, names(Pt_distance), \(x, nm) {
  
  M <- select(metadata, sample, type, salinity)
  r <- left_join(x$residuals, M, by = "sample")
  
  ggplot(r, aes(x = sample, y = value)) +
    geom_hline(yintercept = median(r$value)) + 
    geom_hline(yintercept = quantile(r$value)[c(2, 4)], linetype = 2) +
    geom_point(aes(colour = salinity, shape = type), size = 3.5) +
    labs(title = nm, x = "Sample", y = "Procrustes residuals",
         colour = "Salinity (‰)", shape = "Sample type") +
    scale_colour_viridis_c(option = "D") +
    scale_shape_discrete(labels = str_to_title) +
    theme_bw() +
    theme_modification +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    )
  
})

wrap_plots(plot_Pt_distance, guides = "collect")

## Dispersion ----
plot_Dsp <- map2(Dsp, names(Dsp), \(x, nm) {
  
  ggplot(x, aes(x = CentroidGroup, y = CentroidDistance)) +
    geom_boxplot() +
    labs(title = nm, x = "Sample type", y = "Distance to group centroid",
         colour = "Salinity (‰)", shape = "Sample type") +
    scale_x_discrete(labels = str_to_title) +
    scale_y_continuous(limits = c(0, 1)) +
    theme_bw() +
    theme_modification +
    theme(
      aspect.ratio = 2
    )
  
})

wrap_plots(plot_Dsp[grepl("TPM", names(plot_Dsp))], guides = 'collect')

# Figures ----
## Ordination

fig_ordination <- 
  (plot_O$wgs_TPM | plot_Dsp$wgs_TPM | plot_O$wts_TPM | plot_Dsp$wts_TPM) +
  plot_layout(nrow = 2, ncol = 2, guides = "collect") +
  theme(
    text = element_text(family = "sans")
  )

ggsave(
  filename = "results/Figures/ordination_with_centroid_distances.pdf",
  fig_ordination
)

ggsave(
  filename = "results/Figures/ordination_with_centroid_distances.jpg",
  fig_ordination
)

## Procrustes
fig_procrustes <- 
  (
    plot_Pt$TPM + 
      theme(aspect.ratio = 1)
  ) | 
  (
    plot_Pt_distance$TPM + 
      scale_y_continuous(limits = c(0, 0.2)) +
      theme(aspect.ratio = 1)
  ) +
    plot_layout(guides = "collect")

ggsave(
  filename = "results/Figures/Procrustes_with_residuals.pdf",
  fig_procrustes
)

ggsave(
  filename = "results/Figures/Procrustes_with_residuals.jpg",
  fig_procrustes
)
