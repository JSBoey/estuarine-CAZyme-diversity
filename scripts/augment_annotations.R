# Augment annotations with metadata

# Description ----
# Appends additional information from web sources to the main annotation table.
# Adds the following:
# * dbCAN-sub EC
# * Consensus EC from all annotation methods
# * Labels for transporter super-families

# Environment ----
library(tidyverse)

# Data ----
## Consolidated annotations
ann <- read_tsv("results/allbins_pred.annotation_table.tsv.gz")
## CAZy tables
cazy <- read_tsv("data/cazy_data.tsv.gz")
## dbCAN-sub family-substrate mapping file
dbs_map <- read_tsv("https://bcb.unl.edu/dbCAN_sub/data/fam-substrate-mapping-08252022.tsv")
## dbCAN-sub cluster-EC mapping file
dbs_combined <- read_tsv("https://bcb.unl.edu/dbCAN_sub/data/combined23_6_22.txt")


# Annotate dbCAN-sub EC ----
## Modify dbs_combined for mapping
dbs_combined <- dbs_combined %>% 
  select(`sub-family`, hmm_ec_stats) %>% 
  mutate(
    `sub-family` = str_replace(`sub-family`, "cluster_", "e"),
    hmm_ec_stats = str_remove_all(hmm_ec_stats, "\\[|\\]") %>% 
      str_split(",") %>% 
      map(\(x) {
        unlist(x) %>% 
          str_trim() %>% 
          str_replace("([^:]):.*", "\\1")
      })
  ) %>% 
  arrange(`sub-family`) %>% 
  filter(
    map_int(hmm_ec_stats, \(x) length(x[x != ""])) > 0
  ) %>% 
  rename(
    "dbcansub_label" = "sub-family",
    "dbcansub_ec" = "hmm_ec_stats"
  )

ann_dbcansub <- ann %>% 
  ## Subset annotations to work with dbCAN-sub mapping
  filter(!is.na(dbcansub_label)) %>% 
  select(node, dbcansub_label) %>% 
  ## Modify to allow joining
  mutate(
    dbcansub_label = str_split(dbcansub_label, ";")
  ) %>% 
  unnest(dbcansub_label) %>% 
  ## Join with EC number
  left_join(dbs_combined) %>%
  ## Return to per-node summaries
  group_by(node) %>% 
  summarise(
    dbcansub_label = str_c(dbcansub_label, collapse = ";"),
    dbcansub_ec = list(dbcansub_ec) %>% 
      map_chr(\(x) {
        unlist(x) %>% 
          unique() %>% 
          str_c(collapse = ";")
      })
  ) %>% 
  ## Retain rows with ECs
  filter(dbcansub_ec != "")

## Join back to annotation table
ann <- ann %>% 
  left_join(ann_dbcansub) %>% 
  relocate(dbcansub_ec, .after = dbcansub_label)

# CAZyme consensus EC ----
## Rules:
## 1. If the gene/node contains a dbCAN assignment, it is assumed to be a 
##    CAZyme. Given that this is a metagenomic survey and sequences can be 
##    truncated, I am more averse to false negatives as opposed to false 
##    positives. This thinking is also applied to rule 3.
## 2. CAZy is the arbiter of truth when it comes to the possible EC assignments 
##    for each family.
## 3. Hierarchy of decisions:
##    3.1. Consensus EC = EC(dbCAN-sub) in EC(CAZy)
##         dbCAN-sub has been designed with substrate prediction in mind. 
##         Authors of the dbCAN tool kit divided each CAZy family into 
##         sub-families (eCAMI), then manually curated members of the sub-family
##         by considering sequence members with characterised EC and substrates,
##         and finally building HMM profiles from sub-family alignments. This
##         implies that significant hits to dbCAN-sub are likely to have hit the
##         active domain of the CAZyme, in addition to being internally
##         consistent with how "significance" is being interpreted across the 
##         dbCAN tools.
##    3.2. Consensus EC = EC(KOFam) in EC(CAZy)
##         This is a close second for the reasons: 
##           (a) it has curated EC annotations as part of the metadata 
##           (b) the method of homology is consistent with dbCAN (i.e., HMM) 
##           (c) matches in dbCAN and this means consistency across databases
##         Why is this not prioritised?
##           (a) the main premise of the study is CAZymes. While KO fam is 
##               is extensive, it is not comprehensive when it comes to newly
##               characterised CAZymes and can mis-annotate them as something
##               different.
##           (b) internal consistency within the dbCAN framework is more
##               relevant and important in this work.
##    3.3. Consensus EC = EC(CAZyDB) in EC(CAZy)
##         This is a last resort if no EC assignments are found in the above 
##         matches. Not only are the hits not directly comparable (HMM vs 
##         BLAST), the filter criteria is also less well defined (dbCAN has 
##         specific E-value and coverage criteria, KO fam has score thresholds).
##         The filter criteria used for BLAST-like searches are also set near
##         the twilight zone (30% identity, 70% query coverage). While less
##         reliable, it can sometimes be the only option left, especially when
##         considering less well characterised CAZymes with limited sequence
##         representation in the database.
##    3.4. Consensus EC = EC(CAZy) 
##         If EC(CAZy) = 1 and rule 3.2 is not satisfied, this rule overrides 
##         all other rules. I trust dbCAN's assignments, and if the CAZyme 
##         family is mono-specific, the annotation should reflect the record.

## Pattern definitions
ec_pattern <- "([0-9n-]+\\.){3}[0-9n-]+"
active_cazy_pattern <- "(AA|CE|GH|GT|PL)\\d+(_\\d+)?"

## Note:
## Tabular processing did not work. Trying list-based processing.

## Prepare CAZy for joining
cazy_ec <- cazy %>% 
  select(cazy_family, ec_number) %>% 
  group_by(cazy_family) %>% 
  summarise(
    family_ec = map(
      list(ec_number), unique
    )
  )

## Subset CAZymes from ann (nodes with dbCAN assignments) and select EC columns
ann_dbcan <- ann %>% 
  filter(!is.na(dbcan_label)) %>% 
  select(node, contains(c("kofam", "dbcan", "cazy")))

## Generate consensus EC
## Hierarchy:
## a. EC(dbCAN-sub) in EC(CAZy)
## b. EC(KOFam) in EC(CAZy)
## c. EC(CAZyDB) in EC(CAZy)
## d. If length(EC(CAZy)) == 1 & !b: EC == EC(CAZy)

extractEC <- function(x) {
  # Define regex
  pattern <- "([0-9n-]+\\.){3}[0-9n-]+"
  
  # Error if not character vector
  if (!is.character(x)) stop("Input must be a character vector.")
  
  # Extract
  str_extract_all(x, pattern)
}

extractActiveCAZyme <- function(x) {
  # Define regex
  pattern <- "(AA|CE|GH|GT|PL)\\d+(_\\d+)?"
  
  # Error if not character vector
  if (!is.character(x)) stop("Input must be a character vector.")
  
  # Extract
  str_extract_all(x, pattern)
}

collectObjects <- function(...) {
  dots <- list(...)
  set_names(dots, names(dots))
}

## Extract EC for each relevant column
ec_compare <- ann_dbcan %>% 
  transmute(
    node, dbcan_label,
    # Extract EC
    ec_list = pmap(
      # This list MUST BE IN THIS ORDER
      list(
        "dbcansub_ec" = dbcansub_ec, 
        "kofam_ec" = kofam_ec, 
        "cazydb_ec" = cazydb_ec
      ), 
      collectObjects
    ) %>% 
      map_depth(2, \(x) {
        extractEC(x) %>% 
          unlist() %>% 
          unique()
      }),
    # Get CAZyme EC based on CAZy
    cazy_ec = map(dbcan_label, \(x) {
      # Extract active CAZymes from string
      cazymes <- unlist(extractActiveCAZyme(x))
      # Filter CAZy based on cazymes
      sub_cazy <- filter(cazy, cazy_family %in% cazymes)
      # Extract EC numbers from table
      extractEC(sub_cazy$ec_number) %>% 
        unlist() %>% 
        unique()
    })
  )

## Get consensus EC based on dbCAN assignment and ECs in CAZy
dbcan_consensus_ec <- ec_compare %>% 
  mutate(
    # node, dbcan_label,
    consensus_ec = map2(ec_list, cazy_ec, \(x, ref) {
      # a. EC(dbCAN-sub) in EC(CAZy)
      # b. EC(KOFam) in EC(CAZy)
      # c. EC(CAZyDB) in EC(CAZy)
      # d. If length(EC(CAZy)) == 1 & !(a|b|c): EC == EC(CAZy)
      
      # Remove NA from list
      x <- discard(x, anyNA)
      # Iterative intersections
      y <- map(x, \(l) intersect(l, ref))
      # Remove character(0)
      y <- keep(y, \(x) length(x) > 0)
      
      # Consider condition d
      if (length(ref) == 1 & length(y) == 0) {
        ref
      } else {
        pluck(y, 1)
      }
    })
  )

## Split table into those with and without consensus EC
dbcan_with_consensus_ec <- dbcan_consensus_ec %>% 
  filter(map_lgl(consensus_ec, \(x) !is_null(x)))

dbcan_without_consensus_ec <- dbcan_consensus_ec %>% 
  filter(map_lgl(consensus_ec, \(x) is_null(x)))

## Manually curate consensus EC
## 1. Compare KOfam and CAZy EC
##    Some annotations are (1) functionally equivalent but the EC numbers are not 
##    harmonised; (2) share domain motif but are not functionally equivalent
curate_dbcan_kofam <- dbcan_without_consensus_ec %>% 
  # Add back KOfam information
  left_join(
    select(ann, node, contains("kofam"))
  ) %>% 
  filter(!is.na(kofam_target)) %>% 
  relocate(
    contains("kofam"), .after = ec_list
  )

## Write out unique kofam-dbcan pairs for manual curation
cazy_kofam <- curate_dbcan_kofam %>% 
  select(dbcan_label, contains("kofam"), cazy_ec) %>% 
  mutate(
    cazy_ec = map_chr(cazy_ec, str_c, collapse = " "),
    kofam_ec = str_replace_all(kofam_ec, ";", " ") %>% 
      str_remove_all("NA") %>% 
      str_trim()
  ) %>% 
  distinct() %>% 
  arrange(dbcan_label)

# write_csv(cazy_kofam, "results/cazy_kofam.csv")

## Read in curated list
cazy_kofam_curated <- readxl::read_excel("results/cazy_kofam_curated.xlsx")

## Join consensus EC to curate_kofam to get node IDs
dbcan_kofam_consensus <- list(
  "curate_dbcan_kofam" = curate_dbcan_kofam %>% 
    select(-consensus_ec),
  "cazy_kofam_curated" = cazy_kofam_curated %>% 
    select(dbcan_label, kofam_target, consensus_ec)
) %>% 
  reduce(left_join)

## Compile table of node and consensus EC
node_consensus_ec <- list(
  select(dbcan_with_consensus_ec, node, dbcan_label, consensus_ec),
  select(dbcan_kofam_consensus, node, dbcan_label, consensus_ec)
) %>% 
  map(\(x) mutate(x, consensus_ec = as.list(consensus_ec))) %>% 
  bind_rows() %>% 
  transmute(
    node, dbcan_label,
    cazyme_consensus_ec = map_chr(consensus_ec, \(x) {
      paste(x, collapse = ";")
    })
  ) %>% 
  filter(
    cazyme_consensus_ec != "NA"
  )

## Join consensus EC to annotations
ann <- left_join(ann, node_consensus_ec)

# CAZyme substrate ----
## There are CAZy family, single EC, and corresponding substrate and class in 
## dbs_map. Only assign a substrate if single EC and CAZy family matches. Nodes
## without cazyme_consensus_ec should NOT be assigned a substrate this way and 
## be subject to manual curation.

## `node_consensus_ec` has the required fields and observations, all joining 
## should be relative to this data frame. 

## 'List-ify' columns and strip CAZy sub-family (this is not accounted for in 
## the dbCAN-sub substrate table)
node_substrate <- node_consensus_ec %>% 
  mutate(
    across(-node, \(x) str_split(x, ";")),
    dbcan_label = map(dbcan_label, \(x) str_replace_all(x, "(.*)_.*", "\\1"))
  )

## Prepare dbs_map for list-wise matching
dbs_map_simple <- dbs_map %>% 
  rename_with(str_to_lower) %>% 
  select(-name) %>% 
  relocate(family, ec_number)

## Match tables
dbs_match <- function(cazy_family, consensus_ec, ref_table = dbs_map_simple) {
  # Subset reference table
  ref_table <- ref_table %>% 
    filter(family %in% cazy_family & ec_number %in% consensus_ec)
  
  # Extract substrates (if there are any)
  if (nrow(ref_table) > 0) {
    
    ssimple <- unique(ref_table$substrate_simple) %>% 
      str_c(collapse = "; ")
    
    shigher <- unique(ref_table$substrate_high_level) %>% 
      str_c(collapse = "; ")
    
    # Combine results to tibble
    tibble(
      "substrate_higher" = shigher,
      "substrate_simple" = ssimple
    )
  } else {
    # Return NA for null matches
     tibble(
       "substrate_higher" = NA,
       "substrate_simple" = NA
     )
  }

}

node_substrate <- node_substrate %>% 
  mutate(
    substrate = map2(dbcan_label, cazyme_consensus_ec, dbs_match)
  ) %>% 
  unnest(substrate)

## Curate genes without substrates
ec_without_substrate <- node_substrate %>% 
  filter(is.na(substrate_higher)) %>% 
  mutate(
    across(
      where(is.list), \(data) map_chr(data, \(x) str_c(x, collapse = "; "))
    )
  ) %>% 
  distinct(dbcan_label, cazyme_consensus_ec, substrate_higher, substrate_simple)

## Write out for curation (only do this once)
# write_csv(ec_without_substrate, "results/ec_without_substrate.csv")

## Rejoin curated EC without substrate with node_substrate
ec_without_substrate_curated <- readxl::read_xlsx("results/ec_without_substrate_curated.xlsx") %>% 
  select(-reasoning)

node_with_curated_substrate <- node_substrate %>% 
  filter(is.na(substrate_higher)) %>% 
  mutate(
    across(
      where(is.list), \(data) map_chr(data, \(x) str_c(x, collapse = "; "))
    )
  ) %>% 
  select(-contains("substrate")) %>% 
  left_join(ec_without_substrate_curated, by = c("dbcan_label", "cazyme_consensus_ec"))

node_with_matched_substrate <- node_substrate %>% 
  filter(!is.na(substrate_higher)) %>% 
  mutate(
    across(
      where(is.list), \(data) map_chr(data, \(x) str_c(x, collapse = "; "))
    )
  )

node_with_matched_curated_substrate <- bind_rows(
  node_with_curated_substrate,
  node_with_matched_substrate
)

node_substrate_for_joining <- node_with_matched_curated_substrate %>% 
  select(node, contains("substrate")) %>% 
  mutate(
    across(contains("substrate"), \(x) if_else(x == "NA", NA, x))
  ) %>% 
  rename_with(\(x) str_replace(x, "substrate", "cazyme_substrate"))

### Joining only non-GT CAZymes and renaming substrate to "cazyme_substrate"
ann <- left_join(
  ann, 
  node_substrate_for_joining,
  by = "node"
)

# Gene ontologies: Regulation of gene expression and signal transduction ----
## Requires the following packages on top of tidyverse
library(httr)
library(jsonlite)
library(xml2)

## InterPro2GO 
### Download
interpro2go <- read_delim(
  file = "https://ftp.ebi.ac.uk/pub/databases/interpro/current_release/interpro2go",   delim = " > ", 
  col_names = c("interpro", "go"), 
  comment = "!",
  trim_ws = TRUE
)

### Clean up
interpro2go_clean <- interpro2go %>% 
  transmute(
    interpro_accession = str_extract(interpro, "IPR\\d{6}"),
    interpro_definition = str_remove(interpro, "InterPro:IPR\\d{6} "),
    go_accession = str_extract(go, "GO:\\d{7}"),
    go_definition = str_replace(go, "GO:(.*) ; .*", "\\1")
  ) %>% 
  mutate(
    across(everything(), str_trim)
  )

## Function to request GO descendants
get_GO_descendants <- function(goterms, relations) {
  # Construct URL
  base_url <- "https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/"
  goterms <- str_replace(goterms, ":", "%3A")
  if (length(goterms) > 1) {
    goterms <- paste0(goterms, "%2C", collapse = "")
  }
  if (length(relations) > 1) {
    relations <- paste0(relations, "%2C", collapse = "")
  }
  requestURL <- paste0(
    base_url, goterms, "/descendants?relations=", relations
  )
  
  # Submit request to API
  r <- GET(requestURL, accept("application/json"))
  stop_for_status(r)

  json <- toJSON(content(r))
  fromJSON(json)$results
}

## Get descendants (is_a, part_of) for:
## - GO:0007165 Signal transduction 
## - GO:0010468 Regulation of gene expression
GOdesc_ST_rGE <- get_GO_descendants(
  goterms = c("GO:0007165", "GO:0010468"),
  relations = c("is_a", "part_of")
)

## Subset interpro2go to those of signal transduction and regulation of gene 
## of gene expression and their descendants
## Also, merge InterPro annotations with GO ancestors
main_go_accession <- c("GO:0007165", "GO:0010468")
children_go_accession <- list(
  "ST" = GOdesc_ST_rGE$children[[1]] %>% 
    filter(relation %in% c("is_a", "part_of")) %>% 
    pull(id),
  "rGE" =GOdesc_ST_rGE$children[[2]] %>% 
    filter(relation %in% c("is_a", "part_of")) %>% 
    pull(id)
)
descendants_go_accession <- list(
  "ST" = setdiff(GOdesc_ST_rGE$descendants[[1]], children_go_accession$ST),
  "rGE" = setdiff(GOdesc_ST_rGE$descendants[[2]], children_go_accession$rGE)
)

go_ids <- c(
  main_go_accession,
  unlist(children_go_accession),
  unlist(descendants_go_accession)
) %>% 
  unique()

interpro2go_clean_sub <- interpro2go_clean %>% 
  filter(go_accession %in% go_ids) %>% 
  mutate(
    go_heritage = case_when(
      go_accession %in% main_go_accession ~ "as_is",
      go_accession %in% unlist(children_go_accession) ~ "child_of",
      go_accession %in% unlist(descendants_go_accession) ~ "descendant_of"
    ),
    go_ancestor = case_when(
      go_accession %in% main_go_accession ~ go_accession,
      go_accession %in% c(GOdesc_ST_rGE$descendants[[1]]) ~ "GO:0007165",
      go_accession %in% c(GOdesc_ST_rGE$descendants[[2]]) ~ "GO:0010468",
    )
  )

### Save the relationship
write_tsv(interpro2go_clean_sub, file = "results/interpro2go_heritage.signal_transduction.regulation_of_gene_expression.tsv")

### Collect all InterPro annotations into one column
node_interpro <- ann %>% 
  select(node, contains("interpro_accession")) %>% 
  drop_na() %>% 
  unite(contains("interpro"), col = "interpro_accession", sep = ";") %>% 
  mutate(
    interpro_accession = str_split(interpro_accession, ";") %>% 
      map(\(ids) unique(str_subset(ids, "-", negate = TRUE)))
  )

### Match node with InterPro annotations with GO 
node_interpro_go <- node_interpro %>% 
  mutate(
    transcription_regulation = map(interpro_accession, \(ids) {
      x <- filter(interpro2go_clean_sub, interpro_accession %in% ids) %>% 
        pull(go_ancestor) %>% 
        unique()
      
      y <- str_replace(x, "GO:0010468", "regulation of gene expression") %>%
        str_replace("GO:0007165", "signal transduction") %>% 
        paste(collapse = "; ")

      if (length(x) == 0) {
        return(NA_character_)
      } else {
        return(y)
      }
      
    })
  )

## Merge the transcription_regulation into the main annotation table
node_interpro_go_tojoin <- node_interpro_go %>% 
  select(node, transcription_regulation) %>% 
  mutate(
    transcription_regulation = as.character(transcription_regulation)
  ) %>% 
  drop_na()

ann <- left_join(ann, node_interpro_go_tojoin, by = "node")

# Classify substrates at a higher level ----

## What substrates have been identified?
substrates <- ann %>% 
  select(dbcan_label, contains("cazyme")) %>% 
  filter(str_detect(dbcan_label, "GH|PL")) %>%
  distinct() %>% 
  drop_na()

## Write out for manual curation
# write_tsv(substrates, file = "manual_curation/substrates.tsv")

# Corrections ----


# Save the annotation table ----
save(ann, file = "results/curated_annotation_table.RData")
write_tsv(ann, "results/curated_annotation_table.tsv")

# Misc functions ----
search_annotations <- function(term, field) {
  ann %>% 
    filter(
      if_any(
        contains(field), \(x) str_detect(x, term)
      )
    ) %>% 
    View()
}




