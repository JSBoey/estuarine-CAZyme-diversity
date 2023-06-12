# CAZyme alpha diversity 

# Environment ----
library(vegan)
library(edgeR)
library(tidyverse)

source("scripts/utility_functions.R")
source("scripts/transcriptomes_information_theory.R")

# Utilities ----
# Filter function for callback to subset data to CAZyme genes only
cazyFilter <- \(df, pos) {
  filter(df, Geneid %in% dbcan$query)
}

# Summarise numeric data to various CAZyme classification levels.
sumFamily <- \(tpm, ann, lvl = "subfamily") {
  # Subset table and mutate family column
  A <- ann %>% 
    select(query, label) %>% 
    mutate(
      cazy_class = str_extract(label, "[A-Za-z]+"),
      cazy_family = str_extract(label, "[A-Za-z]+([0-9]+)?"),
      cazy_subfamily = label
    )
  # Join tables
  B <- left_join(A, tpm, by = c("query" = "Geneid"))
  # Group tables
  if (lvl == "class") {
    C <- group_by(B, cazy_class)
  } else if (lvl == "family") {
    C <- group_by(B, cazy_family)
  } else if (lvl == "subfamily") {
    C <- group_by(B, cazy_subfamily)
  }
  # Summarise tables
  summarise(C, across(where(is.numeric), \(x) sum(x)))
}

# Data ----
env_data <- read_tsv("data/sample_metadata.txt")
dbcan <- read_tsv("results/clean.dbcan.tsv.gz")

callback <- DataFrameCallback$new(cazyFilter)

wgs_count <- read_tsv_chunked("results/WGS_clean_count.tsv.gz", callback)
wts_count <- read_tsv_chunked("results/WTS_clean_count.tsv.gz", callback)

wgs_summary <- read_tsv("results/WGS_count.tsv.summary") %>% 
  rename_with(fcNameClean)
wts_summary <- read_tsv("results/WTS_count.tsv.summary") %>% 
  rename_with(fcNameClean)

# Subset environmental metadata ----
salinity_type <- env_data %>% 
  select(sample, type, salinity)

# CAZyme family-level summaries ----
wgs_sum <- sumFamily(
  tpm = select(wgs_count, Geneid, contains(c("Filt", "Sed"))), 
  ann = dbcan, lvl = "family"
)
wts_sum <- sumFamily(
  tpm = select(wts_count, Geneid, contains(c("Filt", "Sed"))), 
  ann = dbcan, lvl = "family"
)

# Present-transcribed (aka family level richness) ----
cazy_class_richness <- list(
  "wgs" = wgs_sum,
  "wts" = wts_sum
) %>% 
  map(\(df) {
    g <- str_extract(df$cazy_family, "[A-Za-z]+")
    m <- column_to_rownames(df, "cazy_family")
    map(split(m, g), specnumber, MARGIN = 2) %>% 
      bind_rows(.id = "cazy_class")
  }) %>% 
  bind_rows(.id = "data_type")

pd1 <- pivot_longer(
  cazy_class_richness, 
  where(is.numeric), 
  names_to = "sample",
  values_to = "S"
) %>% 
  left_join(
    select(env_data, sample, type, salinity)
  ) %>% 
  filter(
    cazy_class != "cohesin"
  ) %>% 
  mutate(
    cazy_class = factor(cazy_class, levels = c("GH", "GT", "CBM", "PL", "CE", "AA")),
    type = factor(type, levels = c("water", "sediment")),
    data_type = factor(data_type, levels = c("wgs", "wts"))
  )

(
  pt1 <- ggplot(pd1, aes(x = salinity, y = S, colour = cazy_class)) +
    geom_point() +
    geom_smooth(
      method = "lm",
      formula = y ~ splines2::naturalSpline(x, df = 3),
      se = FALSE
    ) +
    labs(
      x = "Salinity",
      y = "Number of CAZy families",
      colour = "CAZy class"
    ) +
    facet_grid(
      data_type ~ type,
      labeller = labeller(
        data_type = c(
          "wgs" = "Metagenome",
          "wts" = "Metatranscriptome"
        ),
        type = str_to_title
      )
    ) +
    theme_bw() +
    theme(
      aspect.ratio = 1,
      panel.grid = element_blank()
    )
)

# CAZyme gene and transcript Shannon index ----
tpm <- list(
  "wgs" = wgs_count,
  "wts" = wts_count
) %>% 
  map2(., list("wgs_summary" = wgs_summary, "wts_summary" = wts_summary), ~ {
    normaliseCounts(.x, .y, method = "tpm", use.lib.size = "all")
  })

aggregate_tpm <- tpm %>% 
  map(~ {
    tpm_df <- as.data.frame(.x) %>% 
      rownames_to_column("query")
    dbcan_tpm <- left_join(select(dbcan, query, label), tpm_df) %>% 
      mutate(
        label = str_remove(label, "_.*")
      ) %>% 
      group_by(label) %>% 
      summarise(across(where(is.numeric), ~ sum(.x))) %>% 
      column_to_rownames("label") %>% 
      as.matrix()
    
    return(dbcan_tpm)
  })

aggregate_tpm_sample_type <- with(
  tpm, 
  list(
    "wgs_water" = subset(wgs, select = str_detect(colnames(wgs), "Filt")),
    "wgs_sediment" = subset(wgs, select = str_detect(colnames(wgs), "Sed")),
    "wts_water" = subset(wts, select = str_detect(colnames(wts), "Filt")),
    "wts_sediment" = subset(wts, select = str_detect(colnames(wts), "Sed"))
  )
) %>% 
  map(~ subset(.x, subset = rowSums(.x) > 0))

cazy_shannon <- map(aggregate_tpm_sample_type, ~ {
  H_j <- diversity(.x, MARGIN = 2, base = 2) / log(nrow(.x), base = 2)
  delta_j <- sampleSpecificity(.x, log.base = 2) / log(ncol(.x), base = 2)
  D_j <- KLDivSample(.x, log.base = 2)
  
  data.frame(
    "H_j" = H_j,
    "delta_j" = delta_j,
    "D_j" = D_j
  ) %>% 
    rownames_to_column("sample")
}) %>% 
  bind_rows(.id = "data_type") %>% 
  left_join(select(env_data, sample, salinity, type)) %>% 
  mutate(
    data_type = str_remove(data_type, "_.*")
  )

ggplot(cazy_shannon, aes(x = H_j, y = delta_j, colour = salinity)) +
  geom_text(aes(label = sample)) +
  scale_colour_viridis_c() +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(0, 1)) +
  facet_grid(data_type ~ type) +
  theme_bw() +
  theme(
    panel.grid = element_blank()
  )

ggplot(cazy_shannon, aes(x = salinity, y = D_j, colour = H_j)) +
  geom_text(aes(label = sample)) +
  scale_colour_viridis_c() +
  facet_grid(data_type ~ type)

# CAZyme specificity ----
cazy_specifcity <- map(aggregate_tpm_sample_type, ~ {
  geneSpecificity(.x, log.base = 2) / log(ncol(.x), base = 2)
})

# CAZyme family specificity ----
family_count <- list(
  "wgs" = wgs_sum,
  "wts" = wts_sum
) %>% 
  map(\(x) as.matrix(column_to_rownames(x, "cazy_family")))

family_spec <- map(family_count, \(x) geneSpecificity(x) / log2(ncol(x))) %>% 
  map(~ .x[order(.x)])
family_pa <- map(family_count, \(x) rowSums(decostand(x, "pa")))
family_spec_pa <- map2(family_pa, family_spec, ~ {
  x <- data.frame(
    cf = names(.x),
    N = .x
  )
  y <- data.frame(
    cf = names(.y),
    S_i = .y
  )
  
  left_join(x, y) %>% 
    mutate(
      cc = str_extract(cf, "[A-Za-z]+")
    ) %>% 
    filter(N > 0)
})

ggplot(family_spec_pa$wts, aes(x=N, y=S_i)) +
  geom_point() +
  facet_wrap(~ cc)

family_Hj <- map(family_count, ~ diversity(.x, base = 2, MARGIN = 2) / log2(nrow(.x)))
family_dj <- map(family_count, ~ sampleSpecificity(.x) / log2(ncol(.x)))

family_Hj_dj <- map2(family_Hj, family_dj, ~ {
  x <- data.frame(
    Hj = .x,
    dj = .y
  ) %>% 
    rownames_to_column("sample") %>% 
    left_join(env_data)
})

ggplot(family_Hj_dj$wts, aes(x=Hj, y=dj, colour=salinity)) +
  geom_point() +
  facet_wrap(~ type) +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_colour_viridis_c()
