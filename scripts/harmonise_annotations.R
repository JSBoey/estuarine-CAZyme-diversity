# Harmonizing annotations

# Environment ----
library(tidyverse)

# Utilities ----
countMultipleHits <- function(x, query_pattern) {
  # Column with ORF
  qcol <- x %>% 
    map_lgl(~ all(str_detect(.x, query_pattern))) %>% 
    which() %>% 
    names()
  # Find number of duplicates
  sum(table(x[[qcol]]) > 1)
}

# Data ----
# Headers:
dbcan_header <- c(
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
diamond_header <- c(
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
)
kofam_header <- c(
  "significance",
  "query",
  "target",
  "score_threshold",
  "score",
  "evalue",
  "definition"
)
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
signalp_header <- c(
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
# Metadata from Prodigal predicted ORFs:
orf <- read_tsv("data/allbins_pred.metadata.tsv.gz")

# CAZyme annotations:
# Target = dbCAN v11 HMM profiles
dbcan <- read_tsv("results/allbins_pred.dbcan_parsed.tsv.gz",
                  col_names = dbcan_header)
dbcansub <- read_tsv("results/allbins_pred.dbcan-sub_parsed.tsv.gz",
                     col_names = dbcan_header)
# Target = CAZy sequences
# Not used due to overly large RAM requirements
# Only call by filtered chunks if required later
# cazydb <- read_tsv("results/allbins_pred.CAZyDB_filt.tsv.gz",
#                    col_names = diamond_header)

# KoFam annotations:
# Target = KEGG Orthology (KO)
kofam <- read_tsv("results/allbins_pred.kofam_sig.tsv.gz",
                  col_names = kofam_header)

# Sulfatase annotations:
# Target = Sulfatlas ID
sulfatlas <- read_tsv("results/allbins_pred.sulfatlas_filt.tsv.gz", 
                      col_names = diamond_header)

# Interpro annotations:
# Target: Database sepcific hits
interpro <- read_tsv("results/allbins_pred.interpro.tsv.gz",
                     col_names = interpro_header)

# Transporter annotations:
# Target: TC-DB ID
tcdb <- read_tsv("results/allbins_pred.tcdb_filt.tsv.gz",
                 col_names = diamond_header)

# Signal peptide annotations:
# Target: Prediction
signalp <- read_tsv("results/allbins_pred.signalp_prediction_results.txt.gz",
                    col_names = signalp_header, comment = "#")

# Summarise per ORF ----
## dbCAN ----
### Checks ----
# ORFs with multiple hits
countMultipleHits(dbcan, "NODE") * 100 / nrow(dbcan)

# Check distribution of filtering criteria
# E-value
plot(density(-log10(dbcan$evalue)), 
     main = "Distribution of dbCAN -log10(E-value)")
abline(v = c(-log10(1e-15), -log10(max(dbcan$evalue))))
summary(dbcan$evalue)
# HMM Coverage
plot(density(dbcan$target_coverage * 100), 
     main = "Distribution of dbCAN HMM coverage")
abline(v = c(30, mean(dbcan$target_coverage * 100)))
summary(dbcan$target_coverage)

### Merge hits ----
dbcan_merged <- dbcan %>% 
  arrange(query, query_start) %>% 
  mutate(
    target = str_remove(target, "\\.hmm"),
    target_range = paste(target_start, target_end, sep = "-"),
    query_range = paste(query_start, query_end, sep = "-"),
    target_coverage = round(target_coverage, 3)
  ) %>% 
  group_by(query, query_length) %>% 
  reframe(
    across(
      contains(c("target", "query_", "evalue")), list
    )
  ) %>% 
  mutate(
    across(
      -c(query, query_length),
      ~ map_chr(.x, paste, collapse = ";")
    )
  ) %>% 
  select(query, query_length, query_range, 
         target, target_range,
         evalue, target_coverage)


## dbCAN-Sub ----
### Checks ----
# ORFs with multiple hits
countMultipleHits(dbcansub, "NODE") * 100 / nrow(dbcansub)

# Check distribution of filtering criteria
par(mfrow = c(1, 2))
# E-value
plot(density(-log10(dbcansub$evalue)), 
     main = "Distribution of dbCAN-sub -log10(E-value)")
abline(v = c(-log10(1e-15), -log10(max(dbcansub$evalue))))
summary(dbcansub$evalue)
# HMM Coverage
plot(density(dbcansub$target_coverage * 100), 
     main = "Distribution of dbCAN-sub HMM coverage")
abline(v = c(30, mean(dbcansub$target_coverage * 100)))
summary(dbcansub$target_coverage)

### Merge hits ----
dbcansub_merged <- dbcansub %>% 
  arrange(query, query_start) %>% 
  mutate(
    target = str_remove(target, "\\.hmm"),
    target_range = paste(target_start, target_end, sep = "-"),
    query_range = paste(query_start, query_end, sep = "-"),
    target_coverage = round(target_coverage, 3)
  ) %>% 
  group_by(query, query_length) %>% 
  reframe(
    across(
      contains(c("target", "query_", "evalue")), list
    )
  ) %>% 
  mutate(
    across(
      -c(query, query_length),
      ~ map_chr(.x, paste, collapse = ";")
    )
  ) %>% 
  select(query, query_length, query_range, 
         target, target_range,
         evalue, target_coverage)


## InterPro ----
### Checks ----
# ORFs with multiple hits
countMultipleHits(interpro, "NODE") * 100 / nrow(interpro)

# Check distribution of filtering criteria
par(mfrow = c(1, 2))
# E-value
plot(density(-log10(interpro$evalue)), 
     main = "Distribution of InterPro -log10(E-value)")
abline(v = -log10(max(interpro$evalue)),
       col = "red")
summary(interpro$evalue)
# Note: It is okay if the E-value output for InterProScan is more than expected. For PFAM searches, the filter is the gathering score cut-off. The most reliable way to determine if the hit is appropriate is to check if the hit has an InterPro accession. If it does, it is likely to be a good match despite high E-value.

# Alignment length
plot(density(log(interpro$query_end - interpro$query_start + 1)), 
     main = "Distribution of InterPro log-scale alignment lengths")
abline(
  v = c(
    quantile(log(interpro$query_end - interpro$query_start + 1), .25),
    mean(log(interpro$query_end - interpro$query_start + 1)),
    median(log(interpro$query_end - interpro$query_start + 1)),
    quantile(log(interpro$query_end - interpro$query_start + 1), .75)
  ),
  col = c("black", "green", "chartreuse", "black")
)
summary(interpro$query_end - interpro$query_start + 1)

### Merge hits ----
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


## KoFam ----
### Checks ----
# ORFs with multiple hits
countMultipleHits(kofam, "NODE") * 100 / nrow(kofam)

# Check distribution of filtering criteria
par(mfrow = c(1, 2))
# E-value
plot(density(-log10(kofam$evalue)), 
     main = "Distribution of KoFam -log10(E-value)")
abline(v = -log10(max(kofam$evalue)),
       col = "red")
summary(kofam$evalue)
# Difference between threshold and scores
plot(
  density(
    kofam$score - kofam$score_threshold
  ),
  main = "Distribution of KoFam score - threshold"
)
abline(
  v = c(
    min(kofam$score - kofam$score_threshold),
    mean(kofam$score - kofam$score_threshold),
    median(kofam$score - kofam$score_threshold)
  ),
  col = c("red", "green", "blue")
)
summary(kofam$score - kofam$score_threshold)

### Merge hits ----
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

## SignalP ----
### Checks ----
# ORFs with multiple hits
countMultipleHits(signalp, "NODE") * 100 / nrow(signalp)
