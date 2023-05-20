# Clean annotations: TC-DB

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
filename <- "results/allbins_pred.tcdb_filt.tsv.gz"
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
      evalue < 1e-4 & 
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

# Clean up ----
label_pattern <- "[0-9]+\\.[A-Z]+\\.([0-9]+\\.){2}[0-9S]+"
family_pattern <- "[0-9]+\\.[A-Z]+\\.[0-9]+"

data_clean <- data %>% 
  mutate(
    label = str_extract_all(target, label_pattern),
    family = str_extract_all(target, family_pattern),
    across(c(label, family), ~ map_chr(.x, paste, collapse = " "))
  )

# Best hits ----
best_hits <- data_clean %>% 
  group_by(query) %>% 
  filter(
    bitscore == max(bitscore)
  ) %>% 
  ungroup() %>% 
  distinct()

# Write out ----
write_tsv(data_clean, paste0("results/clean.", dataname, ".tsv.gz"))
write_tsv(best_hits, paste0("results/best.", dataname, ".tsv.gz"))
