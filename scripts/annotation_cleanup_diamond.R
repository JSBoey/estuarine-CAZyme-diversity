# Clean up DIAMOND annotations ----

# This script aims to filter only to the best hits per query that also adheres
# to the following criteria:
# - CAZy E-value < 1e-102
# - General E-value < 1e-4
# - Percent positive >= percent ID
# - Query coverage >= 70%
# Best hits are matches with the highest bit scores


# Environment ----
library(tidyverse)
library(patchwork)

# Data ----
diamond_header <- c(
  "query",
  "target",
  "percent_id",
  "align_length",
  "mismatch",
  "gap_open",
  "query_start",
  "query_end",
  "target_start",
  "target_end",
  "evalue",
  "bitscore",
  "query_coverage",
  "target_coverage",
  "percent_positive"
)

# Utilities ----
callbackCazy <- function(df, pos) {
  df %>% 
    filter(
      evalue < 1e-102 &
        percent_positive >= percent_id & 
        query_coverage >= 70
    )
}

callbackDiamond <- function(df, pos) {
  df %>% 
    filter(
      evalue < 1e-4  &
        percent_positive >= percent_id & 
        query_coverage >= 70
    )
}

findBest <- function(df) {
  df %>% 
    group_by(query) %>% 
    filter(
      bitscore == max(bitscore)
    ) %>% 
    ungroup()
}

mergeHits <- function(df) {
  df %>% 
    mutate(
      query_range = paste(query_start, query_end, sep = "-"),
      target_range = paste(target_start, target_end, sep = "-")
    ) %>% 
    select(-ends_with(c("_start", "_end"))) %>% 
    group_by(query) %>% 
    reframe(
      across(everything(), list)
    ) %>% 
    mutate(
      across(
        where(is.list), ~ map_chr(.x, paste, collapse = ";")
      )
    )
}

checkClass <- function(df) {
  unlist(map(df, class))
}

plotChecks <- function(df) {
  data <- findBest(df)
  
  p1 <- ggplot(data, aes(x = query_coverage, y = target_coverage)) +
    geom_point(aes(colour = -log10(evalue), alpha = log2(align_length))) +
    labs(x = "Query coverage", y = "Target coverage",
         colour = "E-value", alpha = "Alignment length") +
    scale_colour_viridis_c() +
    theme_bw() +
    theme(
      panel.grid = element_blank()
    )
  p2 <- ggplot(data, aes(
    x = align_length, 
    y = percent_positive/percent_id
  )) +
    geom_point(aes(colour = -log10(evalue))) +
    labs(x = "Alignment length", y = "Positive:ID",
         colour = "E-value") +
    scale_colour_viridis_c() +
    theme_bw() +
    theme(
      panel.grid = element_blank()
    )
  
  p1 + p2
}

# Clean up: CAZy ----
cazy <- read_tsv_chunked(
  "results/allbins_pred.CAZyDB_filt.tsv.gz", 
  col_names = diamond_header, 
  callback = DataFrameCallback$new(callbackCazy),
  chunk_size = 50000
)

gc()

cazy_best <- findBest(cazy) %>% 
  mergeHits()

gc()

checkClass(cazy_best)

## Write out ----
write_tsv(cazy_best, "results/allbins_pred.CAZyDB_filt_best.tsv.gz")

## Clear data ----
rm(list = ls(pattern = "cazy"))


# Clean up: TC-DB ----
tcdb <- read_tsv_chunked(
  "results/allbins_pred.tcdb_filt.tsv.gz",
  col_names = diamond_header,
  callback = DataFrameCallback$new(callbackDiamond),
  chunk_size = 50000
)

gc()

tcdb_best <- findBest(tcdb) %>% 
  mergeHits()

gc()

## Write out ----
write_tsv(tcdb_best, "results/allbins_pred.tcdb_filt_best.tsv.gz")

## Clear data ----
rm(list = ls(pattern = "tcdb"))

# Clean up: Sulfatlas ----
sulfatlas <- read_tsv_chunked(
  "results/allbins_pred.sulfatlas_filt.tsv.gz",
  col_names = diamond_header,
  callback = DataFrameCallback$new(callbackDiamond),
  chunk_size = 50000
)

gc()

sulfatlas_best <- findBest(sulfatlas) %>% 
  mergeHits()

gc()

## Write out ----
write_tsv(sulfatlas_best, "results/allbins_pred.sulfatlas_filt_best.tsv.gz")

## Clear data ----
rm(list = ls(pattern = "sulfatlas"))



# Obsolete code ----
# cazy_best <- cazy %>% 
#   mutate(
#     query_range = paste(query_start, query_end, sep = "-"),
#     target_range = paste(target_start, target_end, sep = "-")
#   ) %>% 
#   select(-ends_with(c("_start", "_end"))) %>% 
#   group_by(query) %>% 
#   filter(
#     bitscore == max(bitscore)
#   )
# 
# cazy_merged <- cazy_best %>% 
#   group_by(query) %>% 
#   reframe(
#     across(everything(), list)
#   ) %>% 
#   mutate(
#     across(
#       where(is.list), ~ map_chr(.x, paste, collapse = ";")
#     )
#   )

