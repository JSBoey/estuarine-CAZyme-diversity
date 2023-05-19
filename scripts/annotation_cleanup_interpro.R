# Clean up annotations: Parsed dbCAN outputs 

# Environment ----
library(tidyverse)

# Data ----
interpro_header <- c(
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
interpro <- read_tsv("results/allbins_pred.interpro.tsv.gz",
                     col_names = interpro_header)

# Clean up ----
# Note 1: Because the hits are from 3 different databases, need to merge outputs from each database before being joined to the main 'interpro_merged'.
# Note 2: The field separator for some of the definitions is a semicolon. This is not desirable and after checking, it can be safely converted into a comma.

# Start by applying Note 1 and splitting the data frame
interpro_split <- interpro %>% 
  mutate(
    target_definition = str_replace(target_definition, ";", ",") 
  ) %>% 
  split(., .$db)

interpro_merged <- map2(
  interpro_split, 
  names(interpro_split), 
  function(df, nm) {
    # Arrange and add query range; remove date, status, and MD5
    A <- df %>%
      arrange(query, query_start) %>% 
      mutate(
        query_range = paste(query_start, query_end, sep = "-")
      ) %>% 
      select(-c(query_start, query_end, date, md5, status))
    # Reframe multiple hits to list
    B <- A %>% 
      group_by(query, query_length) %>% 
      reframe(
        across(
          contains(c("query_range", "target", "interpro", "evalue")), list
        )
      ) %>% 
      mutate(
        across(
          -c(query, query_length),
          ~ map_chr(.x, paste, collapse = ";")
        )
      )
    # Add database to target and interpro column names
    C <- B %>% 
      rename_with(
        ~ paste0(str_to_lower(nm), "_", .x),
        contains(c("target", "interpro", "evalue"))
      )
    # Return results
    C
  }
) %>% 
  reduce(., full_join)

# Write out ----
write_tsv(interpro_merged, "results/allbins_pred.interpro_clean.tsv.gz")


