# Diagnose and check CAZyme EC annotations

# Environment ----
library(tidyverse)
library(furrr)
library(data.table)
library(brendaDb)
library(progressr)

# Import data ----
# Filtered and clean hit tables
dbcan <- read_tsv("results/clean.dbcan.tsv.gz")
dbcansub <- read_tsv("results/clean.dbcansub.tsv.gz")
cazy <- read_tsv("results/clean.cazydb.tsv.gz")
kofam <- read_tsv("results/clean.kofam.tsv.gz")

# Mapping files
cazy_map <- read_tsv("data/cazy_data.20240130.tsv")
dbcansub_map <- read_tsv("https://bcb.unl.edu/dbCAN_sub/data/combined23_6_22.txt")
dbcansub_substrate <- read_tsv("https://bcb.unl.edu/dbCAN_sub/data/fam-substrate-mapping-08252022.tsv")
dbcansub_cazy <- read_tsv("data/dbcansub_cazy_map.tsv") %>% 
  rename(cluster = subfamily)
dbcansub_cazy_substrate <- read_tsv("data/dbcansub_cazy_substrate_map.tsv")
ec_map <- ReadBrenda(DownloadBrenda())

# Sanity checks ----
# Questions:
# (a) Can one target have more than one EC?
# (b) How much ambiguity is there in EC numbers? (i.e., 2 dashes? 1 dash?)

## CAZy BLAST results ----
# Each "target" has fields on information separated by pipes (|)
# <genbank_accession>|<(multi-)families>|<(multi-)ecs>
cazy %>% 
  filter(!is.na(ec)) %>% 
  mutate(nbar = str_count(target, "\\|")) %>% 
  pull(nbar) %>% 
  table()
# CAZy targets can have 4-8 fields
# Some targets have domains which do not contribute to EC assignment

cazy %>%
  drop_na() %>% 
  select(target, label, ec) %>% 
  filter(str_count(target, "\\|") == 5) %>% 
  distinct() %>% 
  View()
# For domains that have EC, it is not clear which EC is assigned to which domain.
# Remember that targets in this hit table are proteins, they can have multiple domains and multiple ECs.

str_split(cazy$target, "\\|") %>% 
  unlist() %>% 
  str_subset("^[1-7]\\.") %>% 
  str_subset("-") %>% 
  unique()

# Ambiguity = 1 dash

## dbCAN-sub HMM results ----
# Each "target" has fields on information separated by pipes (|) 
# <dbcan_cluster>|<cazy_family>:<n_seq_in_cluster>|<ec>:<n_seq_with_ec>
dbcansub %>% 
  select(target, label) %>% 
  mutate(nbar = str_count(target, "\\|")) %>% 
  pull(nbar) %>% 
  table()
NF <- 15
dbcansub %>%
  select(query, target, label) %>% 
  filter(str_count(target, "\\|") == NF) %>% 
  distinct() %>% 
  View()
dbcansub %>%
  filter(str_count(target, "\\|") == NF) %>% 
  pull(target) %>% 
  str_split("\\|") %>% 
  unlist() %>% 
  str_subset("^[1-7]\\.") %>% 
  str_remove(":\\d+") %>% 
  unique() %>% 
  str_sort(numeric = TRUE)
# dbCAN-sub targets can have 2-22 fields
# Those with 2 fields (i.e., 1 '|') do not have ECs assigned to them
# dbCAN-sub clusters can have multiple ECs, thus targets with more than 2 fields (i.e., > 2 '|') can have multiple ECs.

str_split(dbcansub$target, "\\|") %>% 
  unlist() %>% 
  str_subset("^[1-7]\\.") %>% 
  str_remove(":\\d+") %>% 
  unique() %>% 
  str_subset("-") %>% 
  str_sort(numeric = TRUE)

# Ambiguity = 2 dashes (transferases 2.4.-.-)

# WARNING:
# EC numbers from CBM dbCAN-sub clusters CANNOT be the sole decider of EC and/or substrate. Use CBM assignments as additional support for EC in active CAZyme families (i.e., AA, CE, GH, GT, PL). For example, if a node is assigned a CBM AND a GH, then the substrate and EC of the GH is narrowed based on the CBM.

## Kofam results ----
# WARNING:
# KoFam results should be the last resort if there is no overlapping EC between the dbCAN EC pool (ref: CAZy) and dbCAN-sub or CAZy BLAST.
unique(kofam$definition) %>% 
  str_sort() %>% 
  str_subset("[1-7]+\\.") %>% 
  length()

unique(kofam$definition) %>% 
  str_sort() %>% 
  str_subset("[1-7]+\\.\\d+\\.") %>% 
  length()

setdiff(
  str_subset(unique(kofam$definition), "[1-7]+\\."),
  str_subset(unique(kofam$definition), "[1-7]+\\.\\d+\\.")
) %>% 
  unique()

unique(kofam$definition) %>% 
  str_subset("\\[") %>% 
  length()

setdiff(
  str_subset(unique(kofam$definition), "\\["),
  str_subset(unique(kofam$definition), "[1-7]+\\.")
)

# EC in Kofam is stored in it's definition. This is usually in the form "[EC:<ec>]"
# For KOs with >1 EC, these are separated with a space.
# Ambiguity = 3 dashes

# Extract relevant data from hit tables ----
## dbCAN HMM
## Only keep CAZy (sub)family
DBCAN <- select(dbcan, query, target) %>% 
  transmute(query,
    subfamily = str_remove_all(target, "(.hmm|_[A-Za-z].+)"),
    family    = str_remove(subfamily, "_\\d+")
  )
rm(dbcan)
gc()

## CAZy BLAST
## Filter for best hit and best hit with EC and then keep target for list expansion
CAZY <- c("best", "best_ec") %>% 
  set_names(., nm = .) %>% 
  map(\(s) {
    # Group data
    r <- group_by(cazy, query)
    # Conditional filter for EC
    if (s == "best_ec") {
      r <- filter(r, str_detect(target, "\\|[1-7]\\."))
    }
    # Filter by minimum E-value and maximum bit score
    r <- slice_min(r, order_by = evalue) %>% 
      slice_max(order_by = bitscore) %>% 
      # Retain relevant columns only
      select(query, target) %>% 
      # Expand columns
      mutate(
        target = str_split(target, "\\|"),
        genbank_accession = map_chr(target, 1),
        family = map(target, \(l) str_subset(l, "^(AA|CBM|CE|GH|GT|PL)\\d+"))
      )
    # Add EC column if it's based on best hit with EC
    if (s == "best_ec") {
      r <- mutate(r, ec = map(target, \(l) str_subset(l, "^[1-7]\\.")))
    }
    # Remove target column
    select(r, -target)
  })
rm(cazy)
gc()

## dbCAN-sub HMM
DBCANSUB <- select(dbcansub, query, target) %>% 
  mutate(
    target = str_split(target, "\\|"),
    cluster = map_chr(target, \(l) str_remove(l[1], ".hmm")),
    family = str_remove(cluster, "_e\\d+"),
    ec = map(target, \(l) {
      r <- str_remove(str_subset(l, "^[1-7]\\."), ":\\d+")
      if (length(r) == 0) r <- NA_character_
      return(r)
    })
  ) %>% 
  select(-target)
rm(dbcansub)
gc()

## KOfam
KOFAM <- select(kofam, query, target, definition) %>% 
  mutate(
    ec = str_extract(definition, "\\[EC:.*\\]") %>% 
      str_remove_all("\\[|EC:|\\]") %>% 
      str_split(" ")
  )
rm(kofam)
gc()

# Extract relevant data from mapping files ----
## dbcansub_substrate is incomplete and requires supplementation
dbcansub_substrate <- dbcansub_substrate %>% 
  bind_rows(
    tribble(
      ~Substrate_high_level, ~Substrate_Simple, ~Family, ~Name, ~EC_Number,
      "peptidoglycan", "thiopeptidoglycan", "PL9", "thiopeptidoglycan lyase", "4.2.2.-",
      "alpha-glucan", "dextran", "GH31", "dextranase", "3.2.1.11",
      "beta-glucan", "beta-1,3-1,4-glucan (laminarin, lichenin and cereal D-glucans)", "GH3", "endo-1,3(4)-beta-glucanase", "3.2.1.6",
      "peptidoglycan", "peptidoglycan", "GH102", "peptidoglycan lytic transglycosylase", "4.2.2.-",
      "peptidoglycan", "peptidoglycan", "GH103", "peptidoglycan lytic transglycosylase", "4.2.2.-",
      "peptidoglycan", "peptidoglycan", "GH23", "peptidoglycan lytic transglycosylase", "4.2.2.-",
      "alpha-glucan", "kojibiose", "GH65", "kojibiose hydrolase", "3.2.1.216",
      "beta-mannan", "beta-mannan", "GH130", "beta-1,2-mannosidase", "3.2.1.197",
      "host glycan", "chondroitin", "PL35", "chondroitin lyase", "4.2.2.-",
      "ulvan", "ulvan", "PL40", "ulvan lyase", "4.2.2.-",
      "ulvan", "ulvan", "PL25", "ulvan lyase", "4.2.2.-",
      "ulvan", "ulvan", "PL24", "ulvan lyase", "4.2.2.-",
      "ulvan", "ulvan", "PL28", "ulvan lyase", "4.2.2.-",
      "alginate", "alginate", "PL39", "Poly-(MG)-lyase / MG-specific alginate lyase", "4.2.2.-",
      "alginate", "alginate", "PL39", "poly(beta-mannuronate) lyase / M-specific alginate lyase", "4.2.2.3",
      "alginate", "alginate", "PL39", "poly(alpha-L-guluronate) lyase / G-specific alginate lyase", "4.2.2.11"
    )
  ) %>% 
  distinct()

## Expand columns for dbcansub_cazy
dbcansub_cazy <- dbcansub_cazy %>% 
  mutate(
    description = str_split(description, "\\|"),
    genbank_accession = map_chr(description, 1),
    family = map(description, \(l) str_subset(l, "^(AA|CBM|CE|GH|GT|PL)\\d+")),
    ec = map(description, \(l) str_subset(l, "^[1-7]\\."))
  ) %>% 
  select(-description)

## Expand columns for cazy_map
cazy_map <- cazy_map %>%
  mutate(
    across(
      c(ec_number, genbank), \(s) str_split(s, " ")
    ),
    family = if_else(
      !is.na(subfamily) & !str_detect(cazy_family, "_"), 
      paste0(cazy_family, "_", subfamily),
      cazy_family
    )
  ) %>% 
  select(family, ec_number, protein_name, taxa_domain, organism, genbank) %>% 
  distinct()

## Clean and expand columns for dbcansub_cazy_substrate
## Supplement with additional substrate information
dbcansub_cazy_substrate <- dbcansub_cazy_substrate %>% 
  mutate(
    substrates = case_when(
      cluster %in% c("GH3_e227", "GH3_e0") & ecs == "3.2.1.6" ~ "beta-glucan",
      cluster %in% c("PL24_e0", "PL28_e0") ~ "ulvan",
      TRUE ~ substrates
    )
  ) %>% 
  mutate(
    across(-cluster, ~ str_replace(.x, "^-$", NA_character_)),
    cluster = str_replace(cluster, "_cluster_", "_e"),
    cazy_protein_id = str_split(cazy_protein_id, ","),
    substrates = str_split(substrates, ", ")
  ) %>% 
  select(-count)

# Matching data ----
unwrap <- function(list) {
  # Turns a list into a non-redundant vector
  unique(unlist(list))
}
poolEC <- function(SUBFAMILY) {
  # Takes the (sub)family and extracts valid EC from the CAZy mapping table
  # If it IS a subfamily, extract EC for subfamily, otherwise extract EC for
  #   entire family.
  pattern <- ifelse(
    grepl("_", SUBFAMILY),
    paste0(SUBFAMILY, '\\b'),
    paste0(SUBFAMILY, '(\\b|_)')
  )
  filter(cazy_map, grepl(pattern, family))$ec_number %>% 
    unwrap()
}
# getCAZY <- function(HITS, NODE, SUBFAMILY) {
#   # Subset CAZY hit tables
#   filter(HITS, query == NODE & map_lgl(family, ~ any(grepl(SUBFAMILY, .x))))
# }
getCAZY <- function(HITS, NODE, SUBFAMILY) {
  # Subset CAZY hit tables
  tb <- setDT(HITS)[
    query == NODE & sapply(family, function(x) any(x %in% SUBFAMILY))
  ]
  as_tibble(tb)
}
getHitStatus <- function(x) {
  # Tests a vector for hit and returns:
  y <- case_when(
    is.null(x) ~ 0, # (no hit)
    is.na(x)   ~ 1, # (hit with NA)
    !is.na(x)  ~ 2  # (valid hit)
  ) %>% 
    unique()
  y[length(y) == 0] <- 0
  y
}
overlapEC <- function(NODE, FAMILY, SUBFAMILY) {
  # Collects ECs from DBCANSUB and CAZY then finds agreement to assign EC
  
  # Pool of possible ECs
  ec_pool <- poolEC(SUBFAMILY)
  
  # Did the node and family have hits against hit tables?
  ec_dbcs <- filter(DBCANSUB, query == NODE & family == FAMILY)$ec %>% 
    unwrap()
  ec_cazy <- getCAZY(CAZY$best_ec, NODE, SUBFAMILY)$ec %>% 
    unwrap()
  
  # Hits must be part of EC pool
  ec_cazy <- intersect(ec_pool, ec_cazy)
  ec_dbcs <- intersect(ec_pool, ec_dbcs)
  
  # Hit conditions: 0 = no hits, 1 = hit w/out ECs, 2 = hit with ECs
  cond_dbcs <- getHitStatus(ec_dbcs)
  cond_cazy <- getHitStatus(ec_cazy)
  
  # No valid hits
  if (all(c(cond_dbcs, cond_cazy) < 2)) return(NA_character_)
  
  # Both have hits
  if (all(c(cond_dbcs, cond_cazy) == 2)) {
    # Find agreement between hits
    ec_hit <- intersect(ec_dbcs, ec_cazy)
    
    # If hits disagree, take from DBCANSUB
    if (length(ec_hit) == 0) ec_hit <- ec_dbcs
    
    return(ec_hit)
  }
  
  # Only one of them have hits
  if (between(cond_dbcs + cond_cazy, 2, 3)) {
    ec_hit <- c(ec_cazy, ec_dbcs)
    return(ec_hit)
  }
}
getKOEC <- function(NODE) {
  # Get EC numbers from KOfam
  ec <- filter(KOFAM, query == NODE)$ec %>% 
    unwrap()
  if (length(ec) == 0) ec <- NA_character_
  return(ec)
}
assignEC <- function(NODE, FAMILY, SUBFAMILY) {
  # Look for overlapping EC between CAZy and dbCAN based hits
  a <- overlapEC(NODE, FAMILY, SUBFAMILY)
  # If no overlapping EC, look for EC based on KO hits
  if (all(is.na(a))) {
    # Even EC from KO must be part of the EC pool from CAZy map
    a <- intersect(poolEC(SUBFAMILY), getKOEC(NODE))
  }
  a <- a[!is.na(a)]
  if (length(a) == 0) a <- NA_character_
  return(a)
}
## Functions for getting hit table annotations ----
getHitsDBCANSUB <- function(NODE, FAMILY) {
  # Get valid dbCAN-sub cluster hits from DBCANSUB hit table
  ht <- filter(DBCANSUB, query == NODE & family == FAMILY)$cluster %>% 
    unwrap()
  if (length(ht) == 0) ht <- NA_character_
  return(ht)
}
getHitsCAZY <- function(NODE, SUBFAMILY, HIT_TYPE) {
  # Get valid Genbank accession from CAZY hit table according to hit types
  ht <- getCAZY(CAZY[[HIT_TYPE]], NODE, SUBFAMILY)$genbank_accession %>% 
    unwrap()
  if (length(ht) == 0) ht <- NA_character_
  return(ht)
}

## Functions for getting substrates ----
getSubstrateFromCluster <- function(CLUSTER, EC) {
  sb <- filter(dbcansub_cazy_substrate, 
         cluster %in% CLUSTER & 
           ecs %in% EC)$substrates %>%
    unwrap()
  # Non-matches are NULL
  if (length(sb) == 0) sb <- NA_character_
  return(sb)
}
getSubstrateNoCluster <- function(FAMILY, EC) {
  sb <- filter(dbcansub_substrate,
               Family %in% FAMILY &
                 EC_Number %in% EC)$Substrate_high_level %>%
    unwrap()
  # Non-matches are character(0)
  if (length(sb) == 0) sb <- NA_character_
  return(sb)
}
getAmbigEnzName <- function(SUBFAMILY, EC, GENBANK) {
  # Find enzyme name for ambiguous EC numbers
  # EC and GENBANK can have lengths > 1
  # SUBFAMILY can be FAMILY, modify pattern accordingly
  pattern <- ifelse(
    grepl("_", SUBFAMILY),
    paste0(SUBFAMILY, '\\b'),
    paste0(SUBFAMILY, '(\\b|_)')
  )
  sb <- filter(cazy_map,
               grepl(pattern, family) &
                 map_lgl(ec_number, ~ any(.x %in% EC)) &
                 map_lgl(genbank, ~ any(.x %in% GENBANK)))$protein_name %>%
    unwrap()
  # Non-matches are character(0)
  if (length(sb) == 0) sb <- NA_character_
  return(sb)
}
getSubstrateAmbiguousEC <- function(SUBFAMILY, EC, GENBANK) {
  # Find substrate or enzyme name for ambiguous EC numbers
  # EC and GENBANK can have lengths > 1
  
  # GENBANK part of dbcansub_cazy_substrate?
  # Output can be character vector, NA, or NULL
  sb <- filter(dbcansub_cazy_substrate, 
    ecs %in% EC & 
      map_lgl(cazy_protein_id, ~ any(.x %in% GENBANK)))$substrates %>%
    unwrap()
  
  # GENBANK has enzyme name?
  # Output can be character vector or NA
  if (all(is.na(sb)) | length(sb) == 0) {
    sb <- getAmbigEnzName(SUBFAMILY, EC, GENBANK)
    if (all(!is.na(sb))) sb <- paste0("{{", sb, "}}")
  }
  
  return(sb)
}
getSubstrate <- function(FAMILY, SUBFAMILY, CLUSTER, GENBANK, EC) {
  # Family and subfamily are characters of length 1
  # Cluster, Genbank accession, and EC can be characters with multiple lengths 
  #   that contain NAs.
  
  # If no EC, forget it
  if (length(EC) == 0 | all(is.na(EC))) return(NA_character_)
  
  # Split analysis into those with full or ambiguous EC
  full_ec <- str_subset(EC, "-", negate = T)
  ambg_ec <- str_subset(EC, "-")
  # Non-matches return character(0)
  
  # Node has full EC assigned
  if (length(full_ec) > 0) {
    # Cluster assigned vs non-assigned
    a <- ifelse(
      all(!is.na(CLUSTER)) & CLUSTER %in% dbcansub_cazy_substrate$cluster,
      getSubstrateFromCluster(CLUSTER, EC),
      getSubstrateNoCluster(FAMILY, EC)
    )
    a <- a[!is.na(a)]
  } else {
    a <- NULL
  }
  
  # Node has ambiguous EC assigned
  if (length(ambg_ec) > 0) {
    # Cluster assigned vs non-assigned
    b <- ifelse(
      all(!is.na(CLUSTER)) & CLUSTER %in% dbcansub_cazy_substrate$cluster,
      getSubstrateFromCluster(CLUSTER, EC),
      getSubstrateAmbiguousEC(SUBFAMILY, EC, GENBANK)
    )
    b <- b[!is.na(b)]
  } else {
    b <- NULL
  }
  
  sb <- unwrap(c(a, b))
  
  if (length(sb) == 0) sb <- NA_character_
  return(sb)
  
}

## Run matches ----
augDBCAN <- vector("list", length = nrow(DBCAN))

# Set up cores
plan(multisession, workers = availableCores() - 2)

with_progress({
  # Set up progress bar
  p <- progressor(steps = length(augDBCAN))
  
  # Get hit matches
  augDBCAN <- future_map(
    seq_along(augDBCAN), \(i) {
      # Update progress bar
      p()
      
      # Variables
      NODE <- DBCAN$query[i]
      FAMILY <- DBCAN$family[i]
      SUBFAMILY <- DBCAN$subfamily[i]
      
      # Assign EC
      EC <- assignEC(NODE, FAMILY, SUBFAMILY)
      
      # Assign best Genbank hits
      GENBANK1 <- getHitsCAZY(NODE, SUBFAMILY, "best")
      
      # Assign best Genbank with EC
      GENBANK2 <- getHitsCAZY(NODE, SUBFAMILY, "best_ec")
      
      # Assign dbCAN-sub cluster
      CLUSTER <- getHitsDBCANSUB(NODE, FAMILY)
      
      # Assign substrates
      SUBSTRATE <- getSubstrate(FAMILY, SUBFAMILY, CLUSTER, GENBANK2, EC)
      
      # Compile results
      result <- list(
        'NODE'             = NODE,
        'FAMILY'           = FAMILY,
        'SUBFAMILY'        = SUBFAMILY,
        'DBCANSUB_CLUSTER' = CLUSTER,
        'CAZY_BEST'        = GENBANK1,
        'CAZY_BEST_EC'     = GENBANK2,
        'EC'               = EC,
        'SUBSTRATE'        = SUBSTRATE
      )
      
      return(result)
    }
  )
})

plan(sequential)

# Compile to data frame
CAZYMES <- augDBCAN %>% 
  purrr::transpose() %>% 
  as_tibble() %>% 
  mutate(
    across(c(NODE, FAMILY, SUBFAMILY), as.character)
  )

# Check and correct matches ----
CAZYMES <- CAZYMES %>% 
  mutate(
    # Updated EC for proteins in GENBANK records
    EC = case_when(
      map_lgl(CAZY_BEST_EC, ~ any("ADD62001.1" %in% .x)) ~ list("3.2.1.6"),
      map_lgl(CAZY_BEST_EC, ~ any("AAO77218.1" %in% .x)) ~ list("3.2.1.24"),
      TRUE ~ EC
    ),
    # Missing SUBSTRATE for:
    # - proteins above
    # - ulvan (PL24, PL28)
    # - peptidoglycan (GH23, GH102, GH103)
    # - peptidoglycan (PL9_4 thiopeptidoglycan lyase)
    SUBSTRATE = pmap(list(SUBSTRATE, FAMILY, EC), ~ {
      if (all(is.na(..1)) & all(!is.na(..3))) {
        a <- getSubstrateNoCluster(..2, ..3)
      } else {
        a <- ..1
      }
      
      if (length(a) == 0) a <- NA_character_
      return(a)
    })
  )

unmatched_substrate <- CAZYMES %>% 
  filter(!is.na(EC) & is.na(SUBSTRATE) & grepl("(PL|GH)", FAMILY))

if (nrow(unmatched_substrate) == 0) cat("IT'S DONE!")

# Curate substrates for enzyme names ----
substrate_pattern <- tribble(
  ~family, ~pattern, ~substrate,
  "GH2", "α-L-arabinopyranosidase", "pectin",
  "GH13", "maltopentaose-forming amylase A", "starch",
  "GH23", "peptidoglycan lytic transglycosylase", "peptidoglycan",
  "GH23", "soluble lytic transglycosylase", "peptidoglycan",
  "GH23", "lytic murein transglycosylase", "peptidoglycan",
  "GH38", "α-mannosidase / α-1,P-mannosidase", "alpha-mannan",
  "GH42", "β-galactosidase / α-L-arabinopyranosidase", "arabinogalactan protein",
  "GH42", "α-L-arabinopyranosidase", "arabinan",
  "GH43", "exo-α-1,5-L-arabinanase", "arabinan",
  "GH43", "α-1,5-L-arabinofuranosidase", "arabinan",
  "GH92", "α-1,4-mannosidase", "host glycan",
  "GH103", "exo-lytic transglycosylase", "peptidoglycan",
  "GH105", "d-4,5-unsaturated α-galacturonidase", "pectin",
  "GH110", "exo-α-galactosidase variant A", "host glycan",
  "GH143", "bimodular 2", "pectin",
  "GH163", "endo-β-N-acetylglucosaminidase cleaving GlcNAc-β-1,2-Man", "host glycan",
  "PL39", "alginate lyase", "alginate",
  "PL40", "ulvan lyase", "ulvan"
)

CAZYMES <- CAZYMES %>% 
  mutate(
    SUBSTRATE = map2(FAMILY, SUBSTRATE, \(f, s) {
      # SUBSTRATE (s) is actually an enzyme name
      # Subset substrate patterns
      sb <- setDT(substrate_pattern)[
        family == f & sapply(pattern, \(x) any(grepl(x, s)))
      ]$substrate
      # Append new substrate without removing enzyme names
      sb <- unwrap(append(s, sb))
      # Return NA if no matches
      if (length(sb) == 0) sb <- NA_character_
      return(sb)
    })
  )

## Manually modify EC and substrate for 
## GH125 [3.2.1.163]: alpha-mannan
## GH55 [3.2.1.58] AFZ68118.1: beta-glucan (laminarin)
CAZYMES <- CAZYMES %>% 
  mutate(
    EC = case_when(
      FAMILY == "GH125" ~ list("3.2.1.163"),
      FAMILY == "GH55" & !is.na(CAZY_BEST_EC) ~ list("3.2.1.58"),
      TRUE ~ EC
    ),
    SUBSTRATE = case_when(
      FAMILY == "GH125" ~ list("alpha-mannan"),
      FAMILY == "GH55" & !is.na(CAZY_BEST_EC) ~ list("beta-glucan"),
      TRUE ~ SUBSTRATE
    )
  )

# Save outputs ----
saveRDS(augDBCAN, file = "data/augmented_dbcan.rds")
saveRDS(CAZYMES, file = "data/curated_cazymes_substrates.rds")
saveRDS(dbcansub_substrate, file = "data/curated_dbcansub_substrate_map.rds")
save.image("rdata/curate_cazymes.RData")

