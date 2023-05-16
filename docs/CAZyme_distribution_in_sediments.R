# CAZyme distribution in sediments

# Environment ----
library(tidyverse)
source("scripts/utility_functions.R")

# Data ----
dbcan <- read_tsv("results/allbins_pred.dbcan_parsed.tsv.gz", col_names = F)
wts_count <- read_tsv("results/WTS_clean_count.tsv.gz") |> 
  select(-contains("Filt")) |> 
  filter(Geneid %in% dbcan$X1)
wts_summary <- read_tsv("results/WTS_count.tsv.summary") |>  
  rename_with(.fn = fcNameClean, .cols = -Status) |> 
  select(-contains("Filt"))

# Subset to depolymerising CAZymes ----
depol_cazyme <- filter(dbcan, str_detect(X3, "GH|PL")) |> 
  select(X1, X3) |> 
  

# Normalise data ----
tpm <- count2TPM(wts_count, wts_summary)
## Remove zero-sums
tpm <- filter(tpm, rowSums(select(tpm, -Geneid)) > 0)

# Top 1% of depolymerising CAZymes ----
top_nodes <- tpm |> 
  filter(Geneid %in% depol_cazyme$X1) |> 
  pivot_longer(-Geneid, names_to = "sample", values_to = "tpm") |> 
  group_by(sample) |> 
  filter(tpm > quantile(tpm, 0.99))

top_cazymes <- left_join(top_nodes, depol_cazyme)

top_tpm <- filter(tpm, Geneid %in% unique(top_nodes$Geneid)) |> 
  left_join(dbcan[, c("X1", "X3")], by = c("Geneid" = "X1"))
