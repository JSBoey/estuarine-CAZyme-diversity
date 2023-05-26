# Consolidate annotations

# Environment ----
library(tidyverse)

# Data ----
dataname <- c(
  "dbcan",
  "dbcansub",
  "kofam",
  "pfam",
  "cdd",
  "tigrfam",
  "cazydb",
  "tcdb",
  "sulfatlas",
  "signalp"
)

metadata <- read_tsv("data/allbins_pred.metadata.tsv.gz")

data <- map(dataname, \(x) {
  if (x %in% c("cazydb", "sulfatlas", "tcdb")) {
    filename <- paste0("results/best.", x, ".tsv.gz")
  } else {
    filename <- paste0("results/clean.", x, ".tsv.gz")
  }
  read_tsv(filename)
}) %>%
  set_names(dataname)

# Collect only relevant columns ----
data_select <- map2(data, names(data), \(df, nm) {
  if (str_detect(nm, "dbcan|sulfatlas")) {
    columns <- c("label")
  } else if (str_detect(nm, "pfam|tigrfam|cdd")) {
    columns <- c("target", "target_definition", 
                 "interpro_accession", "interpro_definition")
  } else if (nm == "kofam") {
    columns <- c("target", "definition", "ec")
  } else if (nm == "cazydb") {
    columns <- c("label", "ec")
  } else if (nm == "tcdb") {
    columns <- c("label", "family")
  } else if (nm == "signalp") {
    columns <- c("prediction", "cs_position")
  }
  
  columns <- c("query", columns)
  
  if (nm %in% c("kofam", "signalp")) {
    a <- arrange(df, query)
  } else {
    a <- arrange(df, query, query_start, desc(query_end))
  }
  
  select(a, all_of(columns))
})

metadata_clean <- metadata %>% 
  select(bin, node, type, start, end, strand, partial, start_type, rbs_motif, 
         rbs_spacer, gc_cont, conf, score)


# Consolidate per query ----
data_consolidate <- map2(data_select, names(data_select), \(df, nm) {
  # Add database to column names
  a <- rename_with(df, \(x) paste(nm, x, sep = "_"), -query)
  
  # Process depend on annotation type
  dmnd <- c("tcdb", "sulfatlas")
  hmms <- c("dbcan", "dbcansub", "kofam", "pfam", "cdd", "tigrfam")
  sgnp <- "signalp"
  cazy <- "cazydb"
  
  # Distinct if necessary
  if (nm %in% c(dmnd, cazy)) {
    a <- distinct(a)
  }
  
  # Group and mutate
  if (nm %in% c(hmms, dmnd)) {
    results <- group_by(a, query) %>% 
      reframe(
        across(
          everything(), \(x) list(x)
        )
      ) %>% 
      mutate(
        across(
          -query, \(x) map_chr(x, paste, collapse = ";")
        )
      )
  } else if (nm %in% sgnp) {
    results <- a
  } else if (nm %in% cazy) {
    results <- group_by(a, query, cazydb_label) %>% 
      reframe(
        across(
          everything(), \(x) list(x)
        )
      ) %>% 
      mutate(
        across(
          -query, \(x) map_chr(x, ~ paste(.x[!is.na(.x)], collapse = ";"))
        )
      ) %>% 
      group_by(query) %>% 
      reframe(
        across(
          everything(), \(x) list(x)
        )
      ) %>% 
      mutate(
        across(
          -query, \(x) map_chr(x, ~ paste(.x[!is.na(.x)], collapse = ";"))
        )
      )
  }
  
  results
})

data_consolidate_clean <- data_consolidate %>% 
  modify_at("cazydb", \(df) {
    df %>% 
      mutate(
        cazydb_label = str_replace_all(cazydb_label, " ", ";"),
        cazydb_ec = str_replace_all(cazydb_ec, "^;+|;+$", "") %>% 
          str_replace_all(" ", ";"),
        across(
          starts_with("cazydb"), \(x) {
            str_split(x, ";") %>% 
              map_chr(~ paste(unique(.x), collapse = ";"))
          }
        )
      )
  }) %>% 
  modify_at("tcdb", \(df) {
    df %>% 
      mutate(
        tcdb_family = str_split(tcdb_family, ";") %>% 
          map_chr(~ paste(unique(.x), collapse = ";"))
      )
  })

# Join to form main table ----
main_table <- append(
  list("metadata" = metadata_clean),
  data_consolidate_clean
) %>% 
  reduce(left_join, by = c("node" = "query"))

# Write out ----
write_tsv(main_table, "results/allbins_pred.annotation_table.tsv.gz")

    