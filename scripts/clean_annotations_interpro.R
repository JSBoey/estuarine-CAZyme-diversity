# Clean annotations: InterPro

# Environment ----
library(tidyverse)

# Data ----
filename <- list.files(path = "results/", pattern = "interpro.*.tsv.gz", full.names = T)
dataname <- str_replace(filename, ".*\\.(.*).tsv.gz", "\\1") %>% 
  str_remove_all("-") %>% 
  str_to_lower()

header <- c(
  "query",
  "md5",
  "query_length",
  "db",
  "target",
  "target_definition",
  "query_start",
  "query_end",
  "evalue",
  "status",
  "date",
  "interpro_accession",
  "interpro_definition"
)

data <- read_tsv(filename, col_names = header, comment = "#")

# Clean up ----
# - Add column names
# - Split by database
# - Write out per database

data_clean <- split(data, data$db) %>% 
  map(~ {
    .x %>% 
      select(-c(db, date, status))
  })
names(data_clean) <- str_to_lower(names(data_clean))

# Write out ----
map2(data_clean, names(data_clean), \(x, y) {
  write_tsv(x, paste0("results/clean.", y, ".tsv.gz"))
})
