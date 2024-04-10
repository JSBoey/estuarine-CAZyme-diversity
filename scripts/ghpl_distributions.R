# GH/PL distribution

# Preamble ----
# Exploring how GH and PL encoding genes and trancripts are distributed at a family level across the environment.

# Environment ----
library(matrixStats)
library(VennDiagram)
library(ggvenn)
library(tidyverse)
library(data.table) 
library(patchwork)


# Import data ----
load("results/transformed_counts.RData")
annotation <- read_tsv("results/curated_annotation_table.tsv")
metadata <- read_tsv("data/sample_metadata.txt")
dbs_map <- read_tsv("https://bcb.unl.edu/dbCAN_sub/data/fam-substrate-mapping-08252022.tsv")

# Filter data ----
# Retain rows with dbCAN assignments and are non-zero across all samples.
dbcan_nodes <- annotation$node[
  str_detect(annotation$dbcan_label, "GH|PL") & !is.na(annotation$dbcan_label)
]
transform_count <- map(transform_count, \(x) {
  x[rownames(x) %in% dbcan_nodes & rowSums(abs(x)) > 0, ]
}) %>% 
  keep_at(c("wgs_TPM", "wts_TPM"))

# Aggregate data ----
dbcan_family <- annotation %>% 
  select(node, dbcan_label) %>% 
  filter(str_detect(dbcan_label, "GH|PL")) %>% 
  transmute(
    node,
    family = str_extract_all(dbcan_label, "(GH|PL)\\d+") %>% 
      map(unique)
  ) %>% 
  unnest(family)

aggregate_count <- map(transform_count, \(x) {
  x <- left_join(dbcan_family, as_tibble(x, rownames = "node")) %>% 
    select(-node)
  
  x <- as.data.table(x)
  
  x[, lapply(.SD, sum, na.rm = TRUE), by = family]
})

# Binarise data ----
binary_count <- map(aggregate_count, \(x) {
  mat <- as.matrix(x[, 2:ncol(x)])
  rownames(mat) <- x$family
  
  mat[mat > 0] <- 1
  mat
})

binary_count_eq <- map(binary_count, \(x) {
  x[, colnames(x) %in% colnames(transform_count$wts_TPM)]
})

# Venn ----
# Sum per sample type per family
venn_data <- map(binary_count_eq, \(x) {
  
  a <- c("Filt", "Sed") %>% 
    set_names(., nm = .)
  
  map(a, \(y) {
    # Split data by sample type
    Y <- x[, grep(y, colnames(x))]
    # Get families that are in each data subset
    rownames(Y)[rowSums(Y) > 0]
  })
  
}) %>% 
  list_flatten()

venn_data <- set_names(venn_data, str_remove(names(venn_data), "_TPM"))





venn_comparisons <- list(
  "water_sediment_capacity" = "wgs",
  "water_capacity_realised" = "Filt",
  "sed_capacity_realised" = "Sed"
) %>% 
  map(\(comp) venn_data[str_detect(names(venn_data), comp)])

venn_comparisons_fill <- list(
  c("#0081CCFF", "#E36200FF"),
  c("#0081CCFF", "#031A6BFF"),
  c("#E36200FF", "#8E3B46FF")
) %>% 
  set_names(nm = names(venn_comparisons))

venn_comparisons_diagram <- map2(
  venn_comparisons, venn_comparisons_fill, 
  \(dat, fill_color) {
    ggvenn(
      dat,
      show_percentage = F,
      fill_color = fill_color,
      stroke_color = "grey30",
      text_size = 5
    )  
  })


## Save plots
# filebase <- "results/Figures/venn_diagram.cazyme_family."
# extension <- c("jpg", "pdf")
# walk2(venn_comparisons_diagram, names(venn_comparisons_diagram), \(venn, nm) {
#   for (ext in extension) {
#     filename <- paste0(filebase, nm, ".", ext)
#     ggsave(filename, venn)
#   }
# })

## All sets
venn_fill <- c(
  "#0072B2FF",
  "#D55E00FF",
  "#56B4E9FF",
  "#E69F00FF"
)

venn_diagram <- ggvenn(
  venn_data, 
  show_percentage = F, 
  fill_color = venn_fill, 
  stroke_color = "grey30", text_size = 7
)

# walk(c("pdf", "jpg"), \(x) {
#   filename <- paste0("results/Figures/cazyme_family_venn_diagram.", x)
#   ggsave(filename = filename, plot = venn_diagram)
# })

## 3 in wgs_Sed
with(
  venn_data,
  setdiff(
    setdiff(wgs_Sed, wgs_Filt),
    wts_Sed
  )
)

## 4 in wts_Sed
with(
  venn_data,
  intersect(
    setdiff(wgs_Sed, wgs_Filt),
    wts_Sed
  )
)

# Boxplot of transcription ----
GH_PL_nodes <- c("GH" = "GH", "PL" = "PL")
GH_PL_TPM <- map(GH_PL_nodes, \(x) {
  nodes <- filter(annotation, str_detect(dbcan_label, x)) %>% 
    pull(node)
  x <- transform_count$wts_TPM[
    rownames(transform_count$wts_TPM) %in% nodes, 
  ]
  x <- as_tibble(x, rownames = "node")
}) %>% 
  bind_rows(.id = "family")

GH_PL_TPMlong <- pivot_longer(GH_PL_TPM, where(is.numeric), 
                              names_to = "sample", values_to = "TPM") %>% 
  filter(TPM > 0) %>% 
  select(-node) %>% 
  left_join(metadata, by = 'sample') %>% 
  mutate(
    type = factor(type, levels = c("water", "sediment"), labels = c("Water", "Sediment"))
  )

fig_boxplot <- ggplot(GH_PL_TPMlong, aes(x = sample, y = TPM)) +
  geom_boxplot(aes(fill = salinity)) +
  scale_y_continuous(trans = "log2") +
  scale_fill_viridis_c(option = "D", alpha = 0.75) +
  labs(fill = "Salinity (‰)", y = "TPM", x = "Sample") +
  facet_grid(rows = vars(family), cols = vars(type), 
             scales = "free_x", space = "free") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

# walk(c("pdf", "jpg"), \(x) {
#   filename <- paste0("results/Figures/GH_PL_distribution_boxplot.", x)
#   ggsave(filename, fig_boxplot)
# })

# Most transcribed CAZymes ----
## Aggregate by CAZy families
aggregate_perc <- map(aggregate_count, \(x) {
  mutate(x, across(where(is.numeric), \(y) proportions(y) * 100))
})

## Check distributions
par(mfrow = c(4, 2))

for (i in 2:30) {
  plot(
    density(log10(aggregate_perc$wts_TPM[[i]])), 
    main = names(aggregate_perc$wts_TPM)[i]
  )
}

## Most are unimodal, with some plateauing instead of perfectly unimodal.

## Most transcribed CAZymes in each environment
## Definition: CAZyme families which collective contribute to at least 5% of all 
##             GH and PL transcribed in each sample.

most_transcribed_threshold <- 5
most_transcribed_ghpl <- aggregate_perc$wts_TPM %>%
  filter(if_any(where(is.numeric), ~ . > most_transcribed_threshold))

## Substrates for most transcribed CAZymes
most_transcribed_substrate <- annotation %>% 
  select(dbcan_label, cazyme_substrate_higher) %>% 
  mutate(
    family = str_extract_all(dbcan_label, "(GH|PL)\\d+") %>% 
      map_chr(str_c, collapse = "; ")
  ) %>% 
  drop_na() %>% 
  filter(
    str_detect(
      family, 
      str_c(str_replace_all(
        most_transcribed_ghpl$family, "^|$", "\\\\b"
      ), collapse = "|")
    )
  ) %>% 
  distinct()

most_transcribed_substrate_curated <- tribble(
  ~family, ~substrate,
  "GH13" , "α-glucan",
  "GH3"  , "Arabinan, β-glucan, N-acetylhexose, polyphenol & xylan",
  "GH23" , "Chitin & peptidoglycan",
  "GH109", "Host glycan",
  "GH57" , "α- & β-glucan",
  "GH103", "Peptidoglycan",
  "GH20" , "Exo-polysaccharide & host glycan",
  "GH97" , "α-glucan",
  "GH73" , "Peptidoglycan",
  "GH16" , "Agarose, β-glucan, host glycan & porphyran",
  "GH77" , "α-glucan",
  "GH30" , "β-fucoside, β-glucan, host glycan & xylan",
  "GH42" , "Arabinogalactan protein",
  "PL3"  , "Pectate"
)

plotdata_most_transcribed_ghpl <- most_transcribed_ghpl %>%
  pivot_longer(cols = where(is.numeric),
               names_to = "sample",
               values_to = "relabd") %>%
  left_join(select(metadata, sample, salinity, type), by = "sample") %>% 
  left_join(most_transcribed_substrate_curated) %>% 
  mutate(
    type = factor(
      type,
      levels = c("water", "sediment"),
      labels = c("Water", "Sediment")
    ),
    family = factor(
      family,
      levels = c(
        # Peptidoglycan
        "GH23", "GH73", "GH103",
        # Alpha-glucan
        "GH13", "GH97", "GH77", "GH57",
        # Beta-glucan
        "GH3", "GH16", "GH30",
        # AGP, pectate, host glycan & LPS
        "GH42", "PL3", "GH109", "GH20"
      )
    )
  )

plot_labels_3 <- most_transcribed_substrate_curated %>% 
  mutate(
    label = str_c(family, " (", substrate, ")")
  ) %>% 
  pull(label, name = family)

plot_fill_3 <- c(
  "#EAC8CA",
  "#e69f00",
  "#009e73",
  "#f0e442",
  "#d55e00",
  "#cc79a7",
  "#005b8e",
  "#b87f00",
  "#8C271E",
  "#c0b634",
  "#6699CC",
  "#4A6C6F",
  "#77CBB9",
  "#33261D"
)

(
  plot_most_transcribed_ghpl <-
    ggplot(plotdata_most_transcribed_ghpl, aes(x = sample, y = relabd)) +
    geom_col(aes(fill = fct_rev(family))) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_fill_manual(
      values = plot_fill_3,
      labels = plot_labels_3
    ) +
    labs(x = "Sample", y = "% total GH and PL", fill = "CAZyme family") +
    facet_grid(cols = vars(type), scales = "free_x", space = "free_x") +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 90, vjust = .5)
    )
)

# walk(c("pdf", "jpg"), \(x) {
#   filename <- paste0("results/Figures/most_transcribed_GHPL.", x)
#   ggsave(filename, plot_most_transcribed_ghpl)
# })

# CAZyme present vs expressed per site ----
richness <- map(aggregate_count, \(x) {
  x <- as.matrix(x[, -1])
  colSums(x > 0)
})

richness_df <- data.frame(
  "status" = c(rep("Present", 36), rep("Transcribed", 30)),
  "sample" = c(names(richness$wgs_TPM), names(richness$wts_TPM)),
  "value" = c(richness$wgs_TPM, richness$wts_TPM)
) %>% 
  left_join(select(metadata, sample, type, salinity), by = "sample") %>% 
  mutate(
    type = factor(
      type,
      levels = c("water", "sediment"),
      labels = c("Water", "Sediment")
    )
  )

(
  plot_richness <- ggplot(richness_df, aes(x = salinity, y = value)) +
    geom_point(aes(shape = status), colour = "grey40") +
    geom_smooth(aes(linetype = status, group = status), 
                method = "loess", se = F) +
    labs(x = "Salinity (ppt)", y = "Number of CAZyme families") +
    facet_grid(cols = vars(type)) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      legend.title = element_blank()
    )
)

# walk(c("pdf", "jpg"), \(x) {
#   filename <- paste0("results/Figures/present_v_transcribed_richness.", x)
#   ggsave(filename, plot_richness)
# })

