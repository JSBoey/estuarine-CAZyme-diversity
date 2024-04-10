# Clean annotations: CAZyDB

# Environment ----
library(tidyverse)

# Data ----
# - Add column name
# - Filter:
#   - E-value < 1e-4 (Global)
#   - E-value < 1e-102 (CAZy)
#   - Percent positive >= Percent identity
#   - Percent identity >= 30
#   - Query coverage >= 70
#   - Bitscore >= 50
# - Add 'label' column
filename <- "results/allbins_pred.CAZyDB_filt.tsv.gz"
dataname <- str_replace(filename, ".*\\.(.*)_.*", "\\1") %>% 
  str_to_lower()

header <- c(
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

callback <- \(x, pos) {
  x %>% 
    filter(
      evalue < 1e-102 & 
        bitscore >= 50 &
        percent_id >= 30 &
        query_coverage >= 70 &
        percent_positive >= percent_id
    )
}

data <- read_tsv_chunked(
  filename, DataFrameCallback$new(callback), col_names = header, 
  chunk_size = 50000, comment = "#"
)

gc()

# Clean up ----
label_pattern <- "(GH|GT|PL|CE|AA|CBM)\\d+(_\\d+)?"
ec_pattern <- "([0-9n-]+\\.){3}[0-9n-]+"

data_clean <- data %>% 
  mutate(
    label = str_replace(target, "^[^|]*\\|", "") %>% 
      str_extract_all(label_pattern),
    ec = str_replace(target, "^[^|]*\\|", "") %>% 
      str_extract_all(ec_pattern),
    across(c(label, ec), ~ map_chr(.x, paste, collapse = " "))
  )

gc()

# Best hits ----
## For best hit, need to keep accession
best_bitscore <- data_clean %>% 
  group_by(query) %>% 
  filter(
    bitscore == max(bitscore)
  ) %>% 
  ungroup() %>% 
  mutate(
    accession_best_hit = str_replace(target, "([^\\|])\\|.*", "\\1")
  ) %>% 
  select(-ec)

## For best hit with EC annotation, need to keep accession and EC
best_ec <- data_clean %>% 
  filter(
    ec != ""
  ) %>% 
  group_by(query) %>% 
  filter(
    bitscore == max(bitscore)
  ) %>% 
  ungroup() %>% 
  mutate(
    accession_best_hit_ec = str_replace(target, "([^\\|])\\|.*", "\\1")
  )

best_hits <- bind_rows(best_bitscore, best_ec) %>% 
  distinct() %>% 
  arrange(query)

# Write out ----
write_tsv(data_clean, paste0("results/clean.", dataname, ".tsv.gz"))
write_tsv(best_hits, paste0("results/best.", dataname, ".tsv.gz"))