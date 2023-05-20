# Clean annotations: dbCAN and dbCAN-sub

# Environment ----
library(tidyverse)

# Data ----
filename <- list.files(path = "results/", pattern = "kofam", full.names = T)
dataname <- str_replace(filename, ".*\\.(.*)_.*", "\\1") %>% 
  str_remove_all("-") %>% 
  str_to_lower()

header <- c(
  "significance",
  "query",
  "target",
  "score_threshold",
  "score",
  "evalue",
  "definition"
)

data <- map(filename, ~ {
  read_tsv(.x, col_names = header, comment = "#")
}) %>% 
  set_names(dataname)

# Clean up ----
# - Add column names
# - Add 'label' column
# - Add 'ec' column
# - Remove extra characters from GT2
data_clean <- map(data, ~ {
  .x %>% 
    select(-significance) %>% 
    mutate(
      ec = str_extract(definition, "\\[EC:(.*)\\]") %>% 
        str_replace("\\[EC:(.*)\\]", "\\1")
    )
})

# Write out ----
map2(
  data_clean,
  names(data_clean),
  \(x, y) write_tsv(x, paste0("results/clean.", y, ".tsv.gz"), col_names = T)
)
