# Clean annotations: SignalP

# Environment ----
library(tidyverse)

# Data ----
filename <- list.files(path = "results/", pattern = "signalp.*.txt.gz", full.names = T)
dataname <- "signalp"

header <- c(
  "query",
  "prediction",
  "other",
  "sp",
  "lipo",
  "tat",
  "tatlipo",
  "pilin",
  "cs_position"
)

data <- map(filename, ~ {
  read_tsv(.x, col_names = header, comment = "#")
}) %>% 
  set_names(dataname)

# Clean up ----
# No cleaning required
data_clean <- data

# Write out ----
map2(
  data_clean,
  names(data_clean),
  \(x, y) write_tsv(x, paste0("results/clean.", y, ".tsv.gz"), col_names = T)
)
