# Clean up annotations: Parsed dbCAN outputs 

# Environment ----
library(tidyverse)

# Data ----
dbcan_header <- c(
  "query",
  "query_length",
  "target",
  "target_length",
  "evalue",
  "target_start",
  "target_end",
  "query_start",
  "query_end",
  "target_coverage"
)
dbcan <- read_tsv("results/allbins_pred.dbcan_parsed.tsv.gz",
                  col_names = dbcan_header)
dbcansub <- read_tsv("results/allbins_pred.dbcan-sub_parsed.tsv.gz",
                     col_names = dbcan_header)

# Clean up: dbCAN ----
dbcan_merged <- dbcan %>% 
  arrange(query, query_start) %>% 
  mutate(
    target = str_remove(target, "\\.hmm"),
    target_range = paste(target_start, target_end, sep = "-"),
    query_range = paste(query_start, query_end, sep = "-"),
    target_coverage = round(target_coverage, 3)
  ) %>% 
  group_by(query, query_length) %>% 
  reframe(
    across(
      contains(c("target", "query_", "evalue")), list
    )
  ) %>% 
  mutate(
    across(
      -c(query, query_length),
      ~ map_chr(.x, paste, collapse = ";")
    )
  ) %>% 
  select(query, query_length, query_range, 
         target, target_range,
         evalue, target_coverage)

# Clean up: dbCAN-sub ----
dbcansub_merged <- dbcansub %>% 
  arrange(query, query_start) %>% 
  mutate(
    target = str_remove(target, "\\.hmm"),
    target_range = paste(target_start, target_end, sep = "-"),
    query_range = paste(query_start, query_end, sep = "-"),
    target_coverage = round(target_coverage, 3)
  ) %>% 
  group_by(query, query_length) %>% 
  reframe(
    across(
      contains(c("target", "query_", "evalue")), list
    )
  ) %>% 
  mutate(
    across(
      -c(query, query_length),
      ~ map_chr(.x, paste, collapse = ";")
    )
  ) %>% 
  select(query, query_length, query_range, 
         target, target_range,
         evalue, target_coverage)

# Write out ---- 
write_tsv(dbcan_merged, "results/allbins_pred.dbcan_parsed_clean.tsv.gz")
write_tsv(dbcansub_merged, "results/allbins_pred.dbcan-sub_parsed_clean.tsv.gz")