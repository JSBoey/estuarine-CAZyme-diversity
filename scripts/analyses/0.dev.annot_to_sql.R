library(RSQLite)
library(arrow)
library(data.table)
library(tidyverse)

setwd("C:/Users/jboe440/OneDrive/Current items/Projects/estuarine-CAZyme-diversity/")

con <- dbConnect(SQLite(), dbname = "sandbox/sql/raw_annotations.sqlite")
files <- list.files(path = "data/curated_annotations/", 
                    pattern = "allbins_pred.*.(tsv|txt).gz", 
                    full.names = TRUE)
names(files) <- c("CAZy", 
                  "dbCANsub", 
                  "dbCAN", 
                  "InterPro", 
                  "KOfam", 
                  "SignalP", 
                  "SulfAtlas", 
                  "TCDB")
headers <- list(
  "diamond" = c(
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
  ),
  "dbcan" = c(
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
  ),
  "ipr" = c(
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
  ),
  "kofamscan" = c(
    "significance",
    "query",
    "target",
    "score_threshold",
    "score",
    "evalue",
    "definition"
  ),
  "signalp" = c(
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
)

for (i in seq_along(files)) {
  
  db_name <- names(files)[i]
  
  hit_type <- case_when(
    db_name %in% c("CAZy", "SulfAtlas", "TCDB") ~ "diamond",
    db_name %in% c("dbCAN", "dbCANsub") ~ "dbcan",
    db_name == "KOfam" ~ "kofamscan",
    db_name == "InterPro" ~ "ipr",
    db_name == "SignalP" ~ "signalp"
  )
  
  cat(paste0("Reading ", db_name, "\n"))
  db <- fread(files[i], col.names = headers[[hit_type]])
  
  cat(paste0("Appending to database\n"))
  dbWriteTable(con, 
               name = db_name, 
               value = db)
}

tbl(con, "CAZy") %>% 
  filter(target_coverage > 80 & query_coverage > 80 & percent_id > 70 & evalue < 1e-10) %>% 
  group_by(target) %>% 
  tally() %>% 
  show_query()
  

