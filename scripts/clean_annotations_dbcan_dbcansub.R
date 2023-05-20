# Clean annotations: dbCAN and dbCAN-sub

# Environment ----
library(tidyverse)

# Data ----
filename <- list.files(path = "results/", pattern = "dbcan", full.names = T)
dataname <- str_replace(filename, ".*\\.(.*)_.*", "\\1") %>% 
  str_remove_all("-") %>% 
  str_to_lower()

header <- c(
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

data <- map(filename, ~ {
  read_tsv(.x, col_names = header, comment = "#")
}) %>% 
  set_names(dataname)
  
# Clean up ----
# - Add column names
# - Add 'label' column
# - Remove extra characters from GT2
data_clean <- map(data, ~ {
  .x %>% 
    mutate(
      label = str_remove(target, "\\.hmm") %>% 
        str_replace("(GT2)_[A-Z].*", "\\1")
    )
}) %>% 
  modify_at("dbcansub", ~ {
    .x %>% 
      mutate(
        label = str_replace(target, "([A-Z0-9]+_e[0-9]+)\\|.*", "\\1")
      )
})

# Write out ----
map2(
  data_clean,
  names(data_clean),
  \(x, y) write_tsv(x, paste0("results/clean.", y, ".tsv.gz"), col_names = T)
)
