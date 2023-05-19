# Clean up KoFam annotations

# Environment ----
library(tidyverse)

# Data ----
kofam_header <- c(
  "significance",
  "query",
  "target",
  "score_threshold",
  "score",
  "evalue",
  "definition"
)

kofam <- read_tsv("results/allbins_pred.kofam_sig.tsv.gz",
                  col_names = kofam_header)

# Clean up ----
kofam_merged <- kofam %>% 
  select(-significance) %>% 
  mutate(
    ec_number = if_else(
      str_detect(definition, "EC:"),
      str_extract(definition, "\\[EC:.*\\]"),
      NA_character_
    ) %>% 
      str_remove_all("[\\[|\\]|EC:]")
  ) %>% 
  group_by(query) %>% 
  reframe(
    across(everything(), list)
  ) %>% 
  mutate(
    across(
      where(is.list),
      ~ map_chr(.x, paste, collapse = ";")
    )
  )

# Write out ----
write_tsv(kofam_merged, "results/allbins_pred.kofam_sig_clean.tsv.gz")