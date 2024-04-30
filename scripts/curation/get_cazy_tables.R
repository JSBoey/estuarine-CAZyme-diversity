# Download CAZy reference tables

# Environment ----
library(tidyverse)
library(furrr)
library(progressr)
source("scripts/download_cazy_tables.R")

# Data ----
## Consolidated annotations
# ann <- read_tsv("results/allbins_pred.annotation_table.tsv.gz")

## Update 30 Jan 2024:
## Better to use dbcan output straight away
dbcan <- read_tsv("results/clean.dbcan.tsv.gz")

# Subset families ----
## CBMs are not included in this list
# cazy_fam <- str_split(ann$dbcan_label, ";") %>%
#   unlist() %>% 
#   unique() %>% 
#   str_sort()

## Update 30 Jan 2024:
## Base cazy_fam on dbcan$label without subfamily
cazy_fam <- unique(str_remove(dbcan$label, "_\\d+")) %>% 
  str_subset("^CBM|cohesin", negate = TRUE) %>% 
  set_names(., .)

# Download tables ----
getTablesSafely <- possibly(getCAZyTables, otherwise = "Not found")
# cazy_table <- map(cazy_fam, \(x) {
#   # Print progress
#   cat(str_glue("{x}\tStarted\t"))
#   # Download tables
#   result <- getTablesSafely(x)
#   # Print progress
#   cat(str_glue("Complete"), "\n")
#   result
# })

## Update 30 Jan 2024:
## Use futures to parallelise process
plan(multisession, workers = availableCores() - 2)
with_progress({
  # Set up progress bar
  p <- progressor(steps = length(cazy_fam))
  # Download tables
  cazy_table <- future_map(cazy_fam, \(x) {
    # Update progressor
    p()
    getTablesSafely(x)
  })
})
plan(sequential)


# Clean up results ----
cazy_table_clean <- cazy_table %>% 
  keep(is.data.frame) %>% 
  bind_rows(.id = "cazy_family") %>% 
  ## Update 30 Jan 2024:
  ## Turn subfamilies into numeric
  mutate(
    subfamily = as.numeric(subfamily)
  )

# Write out environment ----
save.image("results/cazy_data.RData")
write_tsv(cazy_table_clean, "data/cazy_data.20240130.tsv.gz")
