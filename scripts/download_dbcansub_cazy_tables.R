# Collect dbcansub to CAZy protein and substrate data

# Environment ----
library(tidyverse)
library(data.table)
library(furrr)

# Import data ----
# The hit table to limit the search
dbcansub <- read_tsv("results/clean.dbcansub.tsv.gz")

# The mapping file contains all the cluster names
dbcansub_map <- read_tsv("https://bcb.unl.edu/dbCAN_sub/data/combined23_6_22.txt")

# Search vector ----
clusters <- str_extract_all(dbcansub$target, ".*.hmm") %>% 
  unlist() %>% 
  unique() %>% 
  str_remove("\\.hmm") %>% 
  str_replace("_e", "_cluster_")

# Download ----
# URL prefix
url_prefix <- "https://bcb.unl.edu/dbCAN_sub/data/CAZy_with_ECs/"

# Download function for apply
f <- function(s) {
  # Generate URL
  u <- paste0(url_prefix, str_replace(s, "_", "/"), ".txt")
  # Generate message
  message(paste("Parsing", s))
  # Download the fasta file
  read_lines(u) %>% 
    # Extract header only
    str_subset("^>") %>%
    # Remove caret
    str_remove(">")
}

# Set futures
plan(multisession, workers = availableCores() - 4)

# Iterate
result <- future_map(set_names(clusters, clusters), \(x) try(f(x)))

# Stop futures
future:::ClusterRegistry("stop")

# Compile data
fresult <- discard(result, \(x) class(x) == "try-error" | length(x) == 0)

dbcansub_cazy <- tibble(
  "subfamily" = names(fresult),
  "description" = fresult
) %>% 
  unnest(description)

write_tsv(dbcansub_cazy, "data/dbcansub_cazy_map.tsv")

# Download for substrate ----
# URL prefix
url_prefix <- "https://bcb.unl.edu/dbCAN2/download/dbsub_data/subtable_v2/"

# Download function for apply
f <- function(s) {
  # Generate URL
  u <- paste0(url_prefix, str_replace(s, "_", "/"), ".txt")
  # Generate message
  message(paste("Parsing", s))
  # Download the fasta file
  read_tsv(u) %>% 
    rename_with(~ str_to_lower(str_replace_all(.x, " ", "_")))
}

# Set futures
plan(multisession, workers = availableCores() - 4)

# Iterate
result <- future_map(set_names(clusters, clusters), \(x) try(f(x)))

# Stop futures
future:::ClusterRegistry("stop")

# Compile data
fresult <- keep(result, ~ nrow(.x) > 0) %>% 
  map(\(df) mutate(df, across(everything(), as.character))) %>% 
  bind_rows(.id = "cluster")

write_tsv(fresult, "data/dbcansub_cazy_substrate_map.tsv")
