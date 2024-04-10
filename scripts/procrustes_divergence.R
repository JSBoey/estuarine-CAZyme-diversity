# Procrustes divergence

# Preamble ----
# Investigating why there were larger divergence in freshwater and in non-saline and high salinity brackish sediments.

# Environment ----
library(vegan)
library(matrixStats)
library(tidyverse)
library(patchwork)

# Import data ----
load("results/transformed_counts.RData")
annotation <- read_tsv("results/allbins_pred.annotation_table.tsv.gz")
metadata <- read_tsv("data/sample_metadata.txt")

# Filter data ----
dbcan_nodes <- annotation$node[
  str_detect(annotation$dbcan_label, "GH|PL") & !is.na(annotation$dbcan_label)
]
transform_count <- map(transform_count, \(x) {
  x[rownames(x) %in% dbcan_nodes, ]
}) %>% 
  keep_at(at = c("wgs_TPM", "wts_TPM")) %>% 
  map_at("wgs_TPM", \(x) {
    keep <- colnames(transform_count$wts_TPM)
    x[, colnames(x) %in% keep]
  })

# Heatmap ----
heatmap_count <- map2(transform_count, names(transform_count), \(x, nm) {
  colnames(x) <- paste0(colnames(x), nm)
  x
}) %>% 
  do.call("cbind", .)
heatmap_count <- heatmap_count[, order(colnames(heatmap_count))]

heatmap(
  log2(heatmap_count + 0.5), 
  Rowv = NA, 
  Colv = NA, 
  labRow = NA, 
  scale = "none"
)

# Alpha diversity ----
alpha_diversity <- map(transform_count, \(x) {
  
  data.frame(
    "S" = specnumber(x, MARGIN = 2),
    "H1" = renyi(t(x), scales = 1, hill = T),
    "H2" = renyi(t(x), scales = 2, hill = T)
  ) %>% 
    as_tibble(rownames = "sample")
  
}) %>% 
  reduce(left_join, by = "sample", suffix = c("_wgs", "_wts")) %>% 
  left_join(select(metadata, sample, type, salinity))

ggplot(alpha_diversity, aes(x = salinity)) +
  geom_point(aes(y = S_wts), colour = "blue") +
  geom_point(aes(y = S_wgs), colour = "orange") +
  scale_y_log10() +
  facet_wrap(~ type)

# Binarise data ----
binary_count <- map(transform_count, \(x) ifelse(x > 0, 1, 0))

# Proportion expressed ----
csum <- map(binary_count, colSums)
prop_exp <- csum$wts_TPM / csum$wgs_TPM
plot(prop_exp)

# Many models ----
mlm_df <- map(transform_count, \(x) {
  x <- x[, 10:30]
  x <- log(x[rowSums(x) > 0, ] + 0.5)
  x <- as_tibble(as.data.frame(x), rownames = "node")
  pivot_longer(x, -node, names_to = "sample", values_to = "TPM")
}) %>% 
  .[c("wts_TPM", "wgs_TPM")] %>% 
  reduce(left_join, by = c("node", "sample"), suffix = c("_wgs", "_wts")) %>% 
  left_join(select(metadata, sample, type, salinity)) %>% 
  nest(data = c(sample, type, salinity, TPM_wts, TPM_wgs)) %>% 
  mutate(
    lmod = map(data, \(x) {
      try(lm(TPM_wts ~ salinity + TPM_wgs, x))
    }),
    aov = map(lmod, \(x) try(anova(x)))
  )

({
  i <- 10
  print(summary(mlm_df$lmod[[i]]))
  print(mlm_df$aov[[i]])
  par(mfrow = c(2, 2))
  plot(TPM_wts ~ TPM_wgs, mlm_df$data[[i]])
  plot(TPM_wts ~ salinity, mlm_df$data[[i]])
  plot(TPM_wgs ~ salinity, mlm_df$data[[i]])
})
