# Sum and average sediment per sample

# Sediment samples are highly heterogeneous. This obstructs correlation and
# construction of association metrics for network analyses. This script creates
# 2 objects: a sum and an averaged sediment gene/transcript count data per 
# biological sample across spatial replicates.

# Environment ----
library(tidyverse)

# Import data ----
wgs <- read_tsv("results/WGS_clean_count.tsv.gz")
wts <- read_tsv("results/WTS_clean_count.tsv.gz")

# Extract metadata ----
wgs_meta <- select(wgs, !starts_with(c("Filt", "Sed")))
wts_meta <- select(wts, !starts_with(c("Filt", "Sed")))

# Extract sediment data ----
wgs_sed <- select(wgs, Geneid, starts_with("Sed"))
wts_sed <- select(wts, Geneid, starts_with("Sed"))

# Means and sums based on column sample name ----
wgs_consolidated <- wgs_sed %>% 
  pivot_longer(-Geneid, names_to = "sample", values_to = "value") %>% 
  mutate(
    Geneid,
    biol_sample = str_remove(sample, "_\\d$"),
    value
  ) %>% 
  group_by(Geneid, biol_sample) %>% 
  summarise(
    mean = mean(value),
    sum = sum(value)
  )

wts_consolidated <- wts_sed %>% 
  pivot_longer(-Geneid, names_to = "sample", values_to = "value") %>% 
  mutate(
    Geneid,
    biol_sample = str_remove(sample, "_\\d$"),
    value
  ) %>% 
  group_by(Geneid, biol_sample) %>% 
  summarise(
    mean = mean(value),
    sum = sum(value)
  )

wgs_mean <- wgs_consolidated %>% 
  select(-sum) %>% 
  pivot_wider(names_from = "biol_sample", values_from = "mean") %>% 
  left_join(wgs_meta, ., by = "Geneid")

wgs_sum <- wgs_consolidated %>% 
  select(-mean) %>% 
  pivot_wider(names_from = "biol_sample", values_from = "sum") %>% 
  left_join(wgs_meta, ., by = "Geneid")

wts_mean <- wts_consolidated %>% 
  select(-sum) %>% 
  pivot_wider(names_from = "biol_sample", values_from = "mean") %>% 
  left_join(wgs_meta, ., by = "Geneid")

wts_sum <- wts_consolidated %>% 
  select(-mean) %>% 
  pivot_wider(names_from = "biol_sample", values_from = "sum") %>% 
  left_join(wgs_meta, ., by = "Geneid")

# Write out ----
write_tsv(wgs_mean, "results/WGS_sed_mean.tsv.gz")
write_tsv(wts_mean, "results/WTS_sed_mean.tsv.gz")
write_tsv(wgs_sum, "results/WGS_sed_sum.tsv.gz")
write_tsv(wts_sum, "results/WTS_sed_sum.tsv.gz")