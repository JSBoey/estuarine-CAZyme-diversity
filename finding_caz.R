#!/usr/bin/env Rscript

# Finding CAZymes

# Jian Sheng Boey
# 24 September 2024

# Packages ----
library(dplyr)
library(readr)
library(stringr)
library(tidyr)
library(purrr)
library(magrittr)
library(optparse)
library(furrr)
library(dtplyr)
library(KEGGREST)
library(digest)
library(progressr)

setwd("~/nobackup_uoa00348/boey/2024-estuarine-cazyme-diversity/")
load("finding_caz.RData")

# Functions ----
grep_EC <- function(string) {
  pattern <- "[0-9]\\.[0-9-]+\\.[0-9-]+\\.[0-9-A-Za-z]+"
  unlist(str_extract_all(string, pattern))
}
prefilter_blast <- function(x) {
  # Out of the Twilight Zone (Rost, 1999)
  filter(x, pident >= 35 & qcovhsp >= 70)
}
slice_blast <- function(x) {
  # Best match per query
  group_by(x, qseqid) %>% 
    slice_min(evalue, n = 1) %>% 
    slice_max(bitscore, n = 1) %>% 
    ungroup()
}
augment_kegg <- function(x) {
  
  ko_col <- which(sapply(x, function (j) any(str_detect(j, "K\\d{5}")))) %>% 
    names()
  ko <- unique(x[[ko_col]])
  
  klist <- future_map(
    split(ko, ceiling(seq_along(ko)/10)),
    function(k) {
      map(keggGet(k), function(i) i[c("ENTRY", "SYMBOL", "NAME")])
    }
  )
  
  kdf <- map_depth(klist, 2, as_tibble) %>% 
    future_map(list_rbind) %>% 
    list_rbind()
  
  left_join(x, kdf, by = join_by({{ko_col}} == "ENTRY")) %>% 
    rename_with(.fn = str_to_lower)
}
in_characterized_cazy <- function(cazyme_family, ec) {
  
  has_subfam <- str_detect(cazyme_family, "_")
  
  if (isTRUE(has_subfam)) {
    p1 <- str_extract(cazyme_family, "(AA|CBM|CE|GH|GT|PL)\\d+")
    p2 <- as.numeric(str_replace(cazyme_family, ".*_(\\d+)$", "\\1"))
    x <- characterized_cazy %>% 
      filter(family == p1 & subfamily == p2)
  } else {
    x <- characterized_cazy %>% 
      filter(family == cazyme_family)
  }
  
  x <- x %>% 
    pull(ec_number) %>% 
    grep_EC() %>% 
    unique()
  
  intersect(x, ec)
}
arrange_query <- function(x) {
  qcol <- str_subset(colnames(x), "q(uery|seqid)$")
  from <- str_subset(colnames(x), "q(.*)?(start|_from)")
  to <- str_subset(colnames(x), "q(.*)?(end|_to)")
  arrange(x,
          as.numeric(str_replace(.data[[qcol]], "Ww(\\d+).*", "\\1")),
          as.numeric(str_replace(.data[[qcol]], ".*NODE_(\\d+).*", "\\1")),
          as.numeric(str_extract(.data[[qcol]], "\\d+$")),
          .data[[from]],
          .data[[to]])
}
find_cgc <- function(intergenic_distance = 2, min_cazymes = 1) {
  
  cat(str_glue("[{Sys.time()}] Looking for CAZyme gene clusters with an intergenic distance of {intergenic_distance}"), "\n")
  
  cat(str_glue("[{Sys.time()}]\tStep 1: Gathering genes..."), "\n")
  
  # Setup
  q <- unique(c(dbcan$query, tcdom$query, tf$query, tcdb$qseqid, stp$query))
  
  # Find genes
  is_sig <- tibble(
    "query" = q,
    "is_caz" = q %in% dbcan$query,
    "is_tf" = q %in% tf$query,
    "is_tc" = q %in% tcdom$query | q %in% tcdb$qseqid,
    "is_stp" = q %in% stp$query
  )
  
  rm(q)
  
  # Combine tables
  sig_genes <- lazy_dt(is_sig) %>% 
    mutate(
      "scaffold" = str_remove(query, "_\\d+$"),
      "gene_order" = as.numeric(str_extract(query, "\\d+$"))
    ) %>% 
    select(query, scaffold, gene_order, starts_with("is")) %>% 
    as_tibble()
  
  cat(str_glue("[{Sys.time()}]\t\t>>> Gathered {nrow(sig_genes)} genes"), "\n")
  cat(str_glue("[{Sys.time()}]\tStep 2: Finding contiguous genes..."), "\n")
  
  # Find promising scaffolds
  good_scaffolds <- lazy_dt(sig_genes) %>% 
    group_by(scaffold) %>% 
    filter(n() > 1) %>% 
    summarise(
      n_genes = n(),
      across(.cols = starts_with("is_"), .fns = sum, .names = "n_{.col}")
    ) %>% 
    filter(n_is_caz >= !!min_cazymes & n_is_tf + n_is_tc + n_is_stp > 0) %>% 
    as_tibble()
  
  # Filter genes by scaffolds
  sig_genes <- filter(sig_genes, scaffold %in% good_scaffolds$scaffold)
  
  cat(str_glue("[{Sys.time()}]\t\t>>> Found {nrow(sig_genes)} genes across {nrow(good_scaffolds)} scaffolds"), "\n")
  cat(str_glue("[{Sys.time()}]\tStep 3: Forming gene clusters from scaffolds (this might take a while)..."), "\n")
  
  rm(good_scaffolds)
  
  # Consider intergenic distance
  split_scaffolds <- split(sig_genes, sig_genes$scaffold)
  
  clusters <- future_map(split_scaffolds, function(df) {
    
    # Clusters with intergenic distance
    genes <- setNames(df$gene_order, df$query)
    height <- intergenic_distance + 1
    
    hc <- hclust(dist(genes), method = 'single')
    cl <- cutree(hc, h = height)
    
    # Ignore clusters with single genes
    cl <- cl[cl %in% which(table(cl) > 1)]
    genes <- genes[names(genes) %in% names(cl)]
    
    if (length(genes) < 1) return(NULL)
    
    query <- names(genes)
    scaffold <- unique(str_remove(query, "_\\d+$"))
    cluster_id <- paste0(scaffold, "::", cl)
    cluster_md5 <- vapply(cluster_id, digest, character(1))
    
    tibble("query" = names(genes), "cluster_md5" = cluster_md5)
  })
  
  clusters <- list_rbind(compact(clusters))
  n_clusters <- length(unique(clusters$cluster_md5))
  
  cat(str_glue("[{Sys.time()}]\t\t>>> Formed {n_clusters} clusters"), "\n")
  cat(str_glue("[{Sys.time()}]\tStep 4: Filtering for CAZyme gene clusters..."), "\n")
  
  pre_cgc <- inner_join(sig_genes, clusters, by = 'query')
  rm(clusters, n_clusters)
  
  good_clusters <- pre_cgc %>% 
    group_by(cluster_md5) %>% 
    summarise(
      n_genes = n(),
      across(.cols = starts_with("is_"), .fns = sum, 
             .names = "n_{.col}")
    ) %>% 
    filter(n_is_caz >= !!min_cazymes & n_is_tf + n_is_tc + n_is_stp > 0) %>% 
    as_tibble()
  
  cat(str_glue("[{Sys.time()}]\t\t>>> Found {nrow(good_clusters)} CAZyme gene clusters"), "\n")
  
  pre_cgc <- filter(pre_cgc, cluster_md5 %in% good_clusters$cluster_md5)
  rm(good_clusters)
  
  # Padding with intermediate genes
  cat(str_glue("[{Sys.time()}]\tStep 5: Padding clusters to ensure contiguity..."), "\n")
  
  split_clusters <- split(pre_cgc, pre_cgc$cluster_md5)
  padded_clusters <- future_map(split_clusters, function(df) {
    
    x <- df$gene_order
    m <- min(x):max(x)
    p <- setdiff(m, x)
    
    if (length(p) < 1) return(df)
    
    nquery <- paste0(unique(df$scaffold), "_", p)
    ndf <- tibble("query" = nquery,
                  "scaffold" = unique(df$scaffold),
                  "gene_order" = p,
                  "is_caz" = FALSE,
                  "is_tf" = FALSE,
                  "is_tc" = FALSE,
                  "is_stp" = FALSE,
                  "cluster_md5" = unique(df$cluster_md5))
    
    bind_rows(df, ndf)
    
  })
  
  cat(str_glue("[{Sys.time()}]\t\t>>> Done"), "\n")
  cat(str_glue("[{Sys.time()}]\tReturning ordered tibble"), "\n")
  
  list_rbind(padded_clusters) %>%
    arrange(as.numeric(str_replace(query, "Ww(\\d+)_.*", "\\1")),
            as.numeric(str_replace(query, ".*_NODE_(\\d+)_.*", "\\1")),
            gene_order)
  
}

find_query <- function(x) {
  if (is.data.frame(x)) 
    x <- colnames(x)
  x[grepl("q(seqid|uery)", x)]
}

match_cgc2pul <- function(CGC, PUL) {
  
  # There is no best query per cluster filtering
  
  options(dplyr.summarise.inform = FALSE)
  on.exit(options(dplyr.summarise.inform = TRUE))
  
  cat(str_glue("[{Sys.time()}]\tPreparing CGC table"), "\n")
  
  # Prepare CGC table
  columns <- str_subset(names(CGC), "is_")
  replacement <- str_to_upper(str_remove_all(columns, "is_"))
  replacement <- str_replace(replacement, "CAZ", "CAZyme")
  names(replacement) <- columns
  for (j in columns) {
    CGC[[j]] <- ifelse(CGC[[j]], replacement[j], NA_character_)
  }
  CGC <- CGC %>% 
    mutate(
      cgc_gene_type = pmap(list(is_caz, is_tf, is_tc, is_stp),
                           function(...) {
                             s <- c(...)
                             paste(s[!is.na(s)], collapse = ";")
                           }),
      cgc_gene_type = map(str_split(cgc_gene_type, ";"), unlist)
    ) %>% 
    select(-starts_with("is_"))
  
  cat(str_glue("[{Sys.time()}]\tPreparing PUL table"), "\n")
  
  # Prepare PUL table
  PUL <- PUL %>% 
    mutate(
      pul_gene_type = map(
        str_split(sseqid, ":"), 
        function(s) unique(unlist(str_extract_all(s, "(CAZyme|TC|TF|STP)")))
      )
    )
  
  cat(str_glue("[{Sys.time()}]\tPreparing joined table"), "\n")
  
  genes <- inner_join(CGC, PUL, by = join_by("query" == "qseqid"))
  
  # Harmonise gene types
  genes <- genes %>% 
    mutate(
      harmonise_gene_type = map2(cgc_gene_type, pul_gene_type, intersect),
      harmonise_gene_type = map_chr(harmonise_gene_type, function(i) {
        ifelse(length(i) > 0, paste(i), NA_character_)
      })
    )
  
  cat(str_glue("[{Sys.time()}]\tScoring matches"), "\n")
  
  # Scoring per CGC-PUL pair
  # - Unique CGC genes
  # - Unique PUL genes
  # - CAZyme gene hits
  # - Total gene pairs
  # - Sum of bitscore
  # - Bitscore averaged across all unique PUL genes
  scoring <- genes %>% 
    group_by(cluster_md5, pul_id) %>% 
    filter(n() > 1) %>% # Don't waste resources on singletons
    summarise(
      n_uniq_cgc = length(unique(query)),
      n_uniq_pul = length(unique(sseqid)),
      cazyme_cazyme_pairs = sum(grepl("CAZyme", harmonise_gene_type)),
      total_pairs = n(),
      total_bitscore = sum(bitscore),
      average_bitscore = total_bitscore / n_uniq_pul
    )
  
  cat(str_glue("[{Sys.time()}]\tFiltering matches"), "\n")
  
  # Filter based on criteria
  # Unique CGC and PUL genes >= 2 
  # CAZyme-CAZyme pair >= 1 
  # Total gene pairs >= 2 
  # Average bitscore >= 50
  scoring <- scoring %>% 
    filter(
      n_uniq_cgc >= 2 & 
        n_uniq_pul >= 2 & 
        cazyme_cazyme_pairs >= 1 &
        total_pairs >= 2 &
        average_bitscore >= 50
    )
  
  # Best scoring CGC-PUL pair
  best_scores <- scoring %>% 
    slice_max(order_by = total_bitscore, n = 1)
  
  cat(str_glue("[{Sys.time()}]\tReturning best scoring clusters"), "\n")
  
  # Return filtered genes set
  inner_join(genes, best_scores[, c("cluster_md5", "pul_id")],
             by = join_by(cluster_md5, pul_id))
}

match_ec2substrate <- function() {
  
  on.exit(plan(sequential))
  
  # 0. Prepare mapping files
  map_table <- subfam_map %>% 
    mutate(family = str_extract(subfamily, "[A-Z]+\\d+"),
           substrate = str_split(substrate, ", ") %>% 
             map(function(i) str_replace_all(i, " ", "_")))
  
  ## Matching by dbCAN-sub families and EC
  dbsub_map <- distinct(map_table, subfamily, ec_number, substrate)
  
  ## Matching by CAZy sequences in dbCAN-sub families with EC
  cazy_map <- map_table %>% 
    mutate(genbank = str_split(genbank, ",")) %>% 
    unnest(genbank) %>% 
    distinct(family, ec_number, genbank, substrate)
  
  ## Matching by dbCAN families and EC obtained from KOfam and KEGG
  ec_map <- distinct(map_table, family, ec_number, substrate)
  
  # 1. dbCAN-sub
  q1 <- ec_final %>% 
    filter(ec_db == "dbCAN-sub") %>% 
    mutate(
      substrate = future_pmap(
        list(dbsub_family, ec), 
        function(x, y) {
          dbsub_map %>% 
            filter(subfamily == x & ec_number %in% unlist(y)) %>% 
            pull(substrate) %>% 
            unlist()
        }
      )
    )
  
  q2 <- ec_final %>% 
    filter(ec_db == "CAZyDB") %>% 
    mutate(
      substrate = future_pmap(
        list(cazy_accession, ec), 
        function(x, y) {
          cazy_map %>% 
            filter(genbank == x & ec_number %in% unlist(y)) %>% 
            pull(substrate) %>% 
            unlist()
        }
      )
    )
  
  q3 <- ec_final %>% 
    filter(ec_db %in% c("KEGG", "KOfam")) %>% 
    mutate(
      substrate = future_pmap(
        list(dbcan_family, ec), 
        function(x, y) {
          ec_map %>% 
            filter(family == str_extract(x, "[A-Z]+\\d+") & 
                     ec_number %in% unlist(y)) %>% 
            pull(substrate) %>% 
            unlist()
        }
      )
    )
  
  bind_rows(q1, q2, q3)
  
}

wrap <- function(i) unlist(as.list(i[!is.na(i)]))

# Import data ----
# TODO: Becomes opt
ann <- read_tsv("data/collated_dram_annotations.tsv")
cazy <- read_tsv("data/dastool_bins.CAZyDB.07142024.b6o_parsed.gz")
dbcan <- read_tsv("data/dastool_bins.dbCAN-HMMdb-V13.domtblout_parsed.gz")
dbsub <- read_tsv("data/dastool_bins.dbCAN_sub.domtblout_parsed.gz")
stp <- read_tsv("data/dastool_bins.stp.domtblout_parsed.gz")
tcdb <- read_tsv("data/dastool_bins.tcdb_20241918.b6o_parsed.gz")
tcdom <- read_tsv("data/dastool_bins.tcdoms.domtblout_parsed.gz")

tf <- bind_rows(
  read_tsv("data/dastool_bins.tf-1.domtblout_parsed.gz"),
  read_tsv("data/dastool_bins.tf-2.domtblout_parsed.gz")
)
pul <- read_tsv("data/dastool_bins.PUL_20240105.blastp.b6o_parsed.gz")
kegg <- read_tsv("data/dastool_bins.kegg_2024.refprok.b6o_parsed.gz")

characterized_cazy <- read_tsv("db/characterized_cazy.tsv")
subfam_map <- read_tsv("db/subfam_EC_mapping_20240105.tsv",
                       col_names = c("subfamily", 
                                     "ec_number",
                                     "nseq",
                                     "genbank",
                                     "substrate"))
kegg_dat <- read_tsv("db/kegg_2024.refprok.dat", 
                     col_names = c("accession",
                                   "ko",
                                   "num"))
puldb <- read_delim("db/pul_headers.txt", delim = ":",
                    col_names = c("puldb_gene_id",
                                  "puldb_id",
                                  "puldb_gene",
                                  "puldb_gene_locus",
                                  "puldb_accession",
                                  "puldb_gene_type",
                                  "puldb_annotation"))
puldb_meta <- readxl::read_xlsx("db/dbCAN-PUL_12-12-2023.xlsx")
tcdb_meta <- read_tsv("https://tcdb.org/cgi-bin/projectv/public/families.py",
                      col_names = c("transporter_family", 
                                    "transporter_name"))

# Pre-processing ----

## Filter ----
dbcan <- filter(dbcan, evalue < 1e-18 & domain_coverage >= 0.35)
dbsub <- filter(dbsub, evalue < 1e-18 & domain_coverage >= 0.35)
stp <- filter(stp, evalue < 1e-4 & domain_coverage >= 0.3)
tf <- filter(tf, evalue < 1e-4 & domain_coverage >= 0.35)
tcdom <- filter(tcdom, evalue < 1e-4 & domain_coverage >= 0.35)
pul <- prefilter_blast(pul)
tcdb <- prefilter_blast(tcdb) %>% 
  lazy_dt() %>% 
  slice_blast() %>% 
  as_tibble()
kegg <- lazy_dt(kegg) %>% 
  prefilter_blast() %>% 
  slice_blast() %>% 
  as_tibble()

cazy <- (
  function(x) {
    f0 <- function(i) prefilter_blast(i)
    f1 <- function(i) filter(f0(i), map_int(sseqid, \(a) length(grep_EC(a))) > 0)
    
    plan(multisession, workers = 2)
    
    Y <- future_map(list(f0, f1), function(f) {
      lazy_dt(x) %>% 
        f() %>% 
        slice_blast() %>% 
        as_tibble()
    })
    
    plan(sequential)
    
    list_rbind(Y) %>% 
      distinct() %>% 
      arrange(qseqid, qstart, qend) %>% 
      ungroup()
  }
)(cazy)


## Target sequence header expansion ----
dbcan_pattern <- "([A-Z]{2,3}\\d+(_\\d+)?)|([a-z]+)"
cazy_pattern <- "(AA|CBM|CE|GH|GT|PL)\\d+(_\\d+)?"
dbsub_pattern <- "[A-Z]{2,3}\\d+_e\\d+"
tcdb_pattern <- "\\d+\\.[A-Z]\\.(\\d+\\.)+\\d+"

plan(multisession, workers = availableCores() - 2)

dbcan <- dbcan %>% 
  mutate(dbcan_family = str_extract(target, dbcan_pattern))
dbsub <- dbsub %>% 
  mutate(dbsub_family = str_extract(target, dbsub_pattern),
         dbsub_ec = future_map(target, grep_EC))
tcdom <- tcdom %>% 
  mutate(tcdom_family = str_extract(target, "\\d+\\.[A-Z]\\.\\d+"))
cazy <- cazy %>% 
  mutate(cazy_family = str_remove(sseqid, "^[^|]+") %>% 
           str_extract_all(cazy_pattern),
         cazy_accession = str_extract(sseqid, "[^|]+"),
         cazy_ec = future_map(sseqid, grep_EC))
tcdb <- tcdb %>% 
  mutate(tcdb_protein = future_map_chr(strsplit(sseqid, "\\|"), 3),
         tcdb_accession = str_extract(sseqid, tcdb_pattern),
         tcdb_family = str_extract(tcdb_accession, "\\d+\\.[A-Z]\\.\\d+"))
pul <- pul %>% 
  mutate(pul_id = str_extract(sseqid, "PUL\\d+"),
         pul_gene_order = str_replace(sseqid, "PUL\\d+_(\\d+):.*", "\\1") %>% 
           as.numeric(),
         pul_query_contig = str_remove(qseqid, "_\\d+$"),
         pul_query_gene_order = as.numeric(str_extract(qseqid, "\\d+$")))
ann <- ann %>% 
  rename(query = ...1) %>% 
  mutate(kofam_ec = future_map(kegg_hit, grep_EC))
kegg <- left_join(kegg, kegg_dat[, c("accession", "ko")], 
                  by = join_by("sseqid" == "accession")) %>% 
  rename(kegg_ko = ko)

plan(sequential)

# Augment KEGG
tmp.ko <- kegg$kegg_ko
tmp.ko <- tmp.ko[!is.na(tmp.ko)]
tmp.ko <- unique(unlist(strsplit(tmp.ko, "/")))
tmp.kbat <- split(tmp.ko, ceiling(seq_along(tmp.ko)/10))
plan(multisession, workers = availableCores())
with_progress({
  p <- progressor(along = tmp.kbat)
  tmp.klist <- future_map(tmp.kbat, function(k) {
    p()
    map(keggGet(k), function(i) i[c("ENTRY", "SYMBOL", "NAME")])
  })
})
kegg_ko_def <- map_depth(tmp.klist, 2, as_tibble) %>% 
  future_map(bind_rows) %>% 
  bind_rows() %>% 
  rename_with(
    function(s) paste0("kegg_", str_to_lower(s))
  )
plan(sequential)

## Order tables ----
dbcan <- arrange_query(dbcan) %>% 
  select(query, starts_with("dbcan"), everything())
dbsub <- arrange_query(dbsub) %>% 
  select(query, starts_with("dbsub"), everything())
tcdom <- arrange_query(tcdom) %>% 
  select(query, starts_with("tcdom"), everything())
cazy <- arrange_query(cazy) %>% 
  select(qseqid, starts_with("cazy"), everything())
tcdb <- arrange_query(tcdb) %>% 
  select(qseqid, starts_with("tcdb"), everything())
pul <- arrange_query(pul) %>% 
  select(qseqid, starts_with("pul"), everything())
stp <- arrange_query(stp)
tf <- arrange_query(tf)

# CAZyme EC ----
## Assigning EC
## All valid hits to dbCAN are CAZymes
## EC is assigned based on annotations ranked in order of reliability:
## 1. dbCAN-sub with EC assignments are best
## 2. KoFamScan (via DRAM) overlap with characterized CAZymes for the family
## 3. CAZyDB assignments
## 4. KEGG (via DIAMOND blastp) overlap with characterized CAZymes for the family

plan(multisession, workers = availableCores() - 2)

### dbCAN-sub
ec_1 <- inner_join(select(dbcan, query, contains('dbcan')),
                   select(dbsub, query, contains('dbsub')),
                   by = 'query',
                   relationship = 'many-to-many') %>%
  filter(
    # dbCAN-sub and dbCAN CAZy families must match
    str_remove(dbsub_family, "_.*") == str_remove(dbcan_family, "_.*") &
      future_map_int(dbsub_ec, length) > 0
  )

### KoFamScan <> characterized CAZymes
ec_2 <- inner_join(select(dbcan, query, contains('dbcan')),
                   select(ann, query, ko_id)) %>% 
  filter(!is.na(ko_id)) %>% 
  augment_kegg() %>% 
  mutate(
    kofam_ec = future_map(name, grep_EC),
    kofam_ec_in_cazy = future_map2(dbcan_family, kofam_ec, 
                                   in_characterized_cazy)) %>% 
  rename(kofam_symbol = symbol,
         kofam_definition = name)

### CAZyDB
ec_3 <- inner_join(select(dbcan, query, contains('dbcan')),
                   select(cazy, qseqid, starts_with("cazy")),
                   by = join_by("query" == "qseqid"),
                   relationship = "many-to-many") %>% 
  filter(future_map_int(cazy_ec, length) > 0)

### KEGG <> characterized CAZymes
ec_4 <- inner_join(select(dbcan, query, contains('dbcan')),
                   select(kegg, qseqid, starts_with("kegg")),
                   by = join_by("query" == "qseqid"),
                   relationship = "many-to-many") %>% 
  filter(!is.na(kegg_ko)) %>% 
  mutate(kegg_ko = strsplit(kegg_ko, "/")) %>% 
  unnest(kegg_ko) %>% 
  augment_kegg() %>% 
  mutate(kegg_ec = future_map(name, grep_EC),
         kegg_ec_in_cazy = future_map2(dbcan_family, kegg_ec, 
                                       in_characterized_cazy)) %>% 
  rename(kegg_symbol = symbol,
         kegg_definition = name)

plan(sequential)

## Finalise EC by hierarchical rules
ec_joined <- reduce(list(ec_1, ec_2, ec_3, ec_4), full_join) %>% 
  distinct()
ec_final <- ec_joined %>% 
  mutate(
    ec = case_when(
      map_int(dbsub_ec, length) > 0 ~ dbsub_ec,
      map_int(kofam_ec_in_cazy, length) > 0 ~ kofam_ec_in_cazy,
      map_int(cazy_ec, length) > 0 ~ cazy_ec,
      map_int(kegg_ec_in_cazy, length) > 0 ~ kegg_ec_in_cazy
    ),
    ec_db = case_when(
      map_int(dbsub_ec, length) > 0 ~ "dbCAN-sub",
      map_int(kofam_ec_in_cazy, length) > 0 ~ "KOfam",
      map_int(cazy_ec, length) > 0 ~ "CAZyDB",
      map_int(kegg_ec_in_cazy, length) > 0 ~ "KEGG"
    )
  ) %>% 
  filter(map_int(ec, length) > 0) %>% 
  arrange(
    as.numeric(str_replace(query, "Ww(\\d+).*", "\\1")),
    as.numeric(str_replace(query, ".*_NODE_(\\d+)_.*", "\\1"))
  ) %>% 
  select(query, dbcan_family, ec, ec_db, everything())

ec_final %>% 
  mutate(across(where(is.list), 
                function(j) map_chr(j, paste, collapse = " "))) %>% 
  write_tsv(file = "data/ec_final.tsv")

rm(ec_joined)

# CGC prediction ----
# CGC = (2+ CAZymes) + (TF | TC | STP); d <= {2..10}
plan(multisession, workers = availableCores() - 2)
cgc <- find_cgc(intergenic_distance = 2, min_cazymes = 1)
plan(sequential)

# Match CGC to characterized PULs
cgc_pul <- match_cgc2pul(cgc, pul)
plan(sequential)

par(mfrow = c(1, 2))
hist(table(cgc$cluster_md5), 
     main = "Genes in CGC", 
     xlab = "Number of genes")
hist(table(cgc_pul$cluster_md5), 
     main = "Genes in CGC matching PUL", 
     xlab = "Number of genes")

# Substrate prediction ----
# 1. Matched to known PULs (best)
# 2. Majority voting in CGCs (better)
# 3. EC matching (meh)

## EC matching ----
### This is done first as (2) and (3) needs it
### Match
substrates_by_ec <- match_ec2substrate()

### Clean up
substrates_by_ec <- substrates_by_ec %>% 
  arrange(
    as.numeric(str_replace(query, "Ww(\\d+).*", "\\1")),
    as.numeric(str_replace(query, ".*_NODE_(\\d+)_.*", "\\1"))
  ) %>% 
  mutate(substrate = map(substrate, function(x) x[x != "-"])) %>% 
  select(query, dbcan_family, ec, ec_db, substrate, everything())

### Augment with manual curation of un-matched EC-substrate pairs

ec_no_substrate <- substrates_by_ec %>% 
  filter(map_int(substrate, length) == 0) %>% 
  mutate(across(where(is.list), function(j) map_chr(j, paste, collapse=" ")))

# GH
### 4.2.2.29 [GH23 GH102 GH103] - peptidoglycan
### 3.2.1.37 [GH39 GH5_22] - xylan
### 2.4.99.16 [GH13_3] - alpha-glucan synthesis from trehalose
### 3.2.1.23 [GH2 WLC21218.1] - possibly lactose (https://doi.org/10.1016/j.pep.2022.106074)
### 2.4.1.397 [GH189] - cyclic beta-1,2-glucan (important in bacterial infection of and symbiosis with hosts) synthesis (PMC10901819, PMC10834661)
### 3.2.1.67 [GH28] - pectin
### 3.2.1.39 [GH157 GH128] - beta-1,3-glucans (curdlan in GH128, an exopolysaccharide of soil bacteria) (https://doi.org/10.1186/s13068-024-02494-5)
### 2.4.1.97 [GH149] - beta-glucan (specifically beta-1,3-glucan) (https://doi.org/10.1016/j.chembiol.2019.03.017)
### 3.2.1.41 [GH57] - pullulan (non-starch alpha glucan) (https://doi.org/10.1007/s12010-010-9029-6)
### 3.2.1.153 [GH32] - inulin (fructan) (https://doi.org/10.1016/j.jbiosc.2014.08.003)
### 3.2.1.199 [GH188] - sulfoquinovosyl glycerol from α-d-sulfoquinovosyl diacylglyceride, a very abundant sulfolipid in plants (classified as polyphenol by dbCAN-sub) (https://doi.org/10.1021%2Fjacs.3c11126)
### 3.2.1.197 [GH130_3] - terminal beta-1,2-mannosidase in Candida alpha-mannan (https://doi.org/10.1074/jbc.m115.681460)
### 3.2.1.40
###     [GH106 ACM23671.1] alpha-rhamnosides of flavonones (e.g., naringin/hesperidin from citrus) as well as potentially other rhamnosylated oligos or glycoconjugates, but not polysaccharides (https://doi.org/10.1038%2Fs41598-019-52251-0)
###     [GH78 CDF79937.1] part of the ulvan degradation cascade (https://doi.org/10.1016/j.jbc.2021.101210)
### 3.2.1.216 [GH65] - Kojibiose, a disaccharide
### 3.2.1.18 [GH177] - sialic acid in host (https://doi.org/10.1099/jmm.0.05349-0)
### 3.2.1.81 [GH50] - agar (https://doi.org/10.3389%2Ffmicb.2018.00839)
### 3.2.1.32 [GH26] - xylan
### 3.2.1.214 2.4.1.391 [GH35] - beta-1,2-glucans (speculated) (https://doi.org/10.1016/j.jbc.2022.101606)
### 2.4.1.49 [GH94] - cellulose
### 3.2.1.205 [GH15] - cyclobis-(1->6)-alpha-nigerosyl (CNN), an intermediate product of starch degradation to glucose by Kribbella flavida (https://doi.org/10.1074%2Fjbc.M116.727305)
### 4.2.2.n1 [GH103] - peptidoglycan
### 3.2.1.53 [GH123_4] - various glycoconjugates (https://doi.org/10.1038/s41467-024-47653-2)
### 3.2.1.55 3.2.1.146 [GH159] - arabinogalactan and arabinoxylan-derived oligosaccharides (https://doi.org/10.3389%2Ffmolb.2022.907439)
### 3.2.1.1 [GH13_19] - starch
### 3.2.1.4 [GH26] - cellulose
### 3.2.1.28 [GH15] - trehalose
### 3.2.1.-
###   [GH16_13 QEP52094.1] - carrageenan
###   [GH172 EGK01279.1 EGK01275.1 EGK01287.1] - D-arabinan common in mycobacterial arabinogalactan (part of cell wall) (https://doi.org/10.1038/s41467-023-37839-5)
###   [GH178 VCT92674.1] - 3-O-methylmannose polysaccharide from mycobacteria (https://doi.org/10.1038%2Fs42003-023-04448-3)
###   [GH26 WP_077338539.1] - algal xylan of mixed beta-1,3/1,4-linkages (https://doi.org/10.1016/j.jbc.2023.105116)
###   [GH97 AMB20699.1] - dextran + change to 3.2.1.70 (https://doi.org/10.1016/j.jsb.2016.09.009)
### 2.4.1.- [GH31_12] - biosynthesis of cyclic tetrasaccharide, a bacterial alpha-glucan

# PL
### 4.2.2.2 [PL1_2] - pectin
### 4.2.2.3 [PL31] - alginate
### 4.2.2.3 4.2.2.11 [PL39] - alginate
### 4.2.2.-
###   [PL9_4] - bacterial sheath polysaccharides (https://doi.org/10.1128/aem.66.11.4998-5004.2000)
###   [PL24] - ulvan
###   [PL28] - ulvan
### 4.2.2.- 4.2.2.- [PL25 KEZ94292.1] - ulvan

# CE
### 3.1.1.72 3.1.1.73 [CE1] - xylan and hemicellulosic components
### 3.2.1.8 [CE1] - xylan (missing the xylanase portion of the protein?)
### 3.2.1.73 [CE1] - xylan (speculative aid in removing feruloyl residues from xylans) (https://doi.org/10.1016%2Fj.jbc.2021.100500)
### 3.5.1.104 [CE4] - peptidoglycan
### 4.2.2.- [CE8] - pectin + should be 3.1.1.11?
### 3.1.1.86 [CE12] - pectin (rhamnogalacturonan) (https://doi.org/10.1074/jbc.270.45.27172)
### 3.1.1.- 3.1.1.72 [CE15] - lignocellulose (https://doi.org/10.1007/s10930-024-10221-0, https://doi.org/10.1186/s13068-018-1213-x)
### 3.1.1.- [CE15] - lignocellulose
### 3.2.1.40 [CE15] - xylan? (https://doi.org/10.1016/j.jbiotec.2014.04.024)
### 3.1.1.117 [CE15] - lignocellulose
### 3.1.1.53 [CE20] - sialic acid

manually_curated_substrate_from_ec <- read_csv("manually_curated_substrates_from_ec.csv")

ec_substrate_curated <- ec_no_substrate

for (i in 1:nrow(manually_curated_substrate_from_ec)) {
  a <- manually_curated_substrate_from_ec$dbcan_family[i]
  b <- manually_curated_substrate_from_ec$ec[i]
  cz <- manually_curated_substrate_from_ec$cazy_accession[i]
  sb <- manually_curated_substrate_from_ec$substrate[i]
  
  ex <- "dbcan_family == a & ec == b"
  
  if (!is.na(cz))
    ex <- paste0(ex, " & cazy_accession == cz")
  
  nr <- which(eval(parse(text = ex), ec_substrate_curated))
  ec_substrate_curated[nr, "substrate"] <- sb
  
  rm(i, a, b, cz, sb, ex, nr)
}

tmp.is_list <- names(substrates_by_ec)[sapply(substrates_by_ec, is.list)]
ec_substrate_curated <- ec_substrate_curated %>%
  mutate(across(any_of(tmp.is_list), function(j) str_split(j, " ")))
rm(list = ls(pattern = "tmp"))

substrates_by_ec_final <- bind_rows(
  filter(substrates_by_ec, map_int(substrate, length) > 0),
  ec_substrate_curated
)

if (all(dim(substrates_by_ec_final) == dim(substrates_by_ec))) {
  substrates_by_ec_final <- substrates_by_ec_final %>%
    filter(map_int(substrate, function(i) length(str_subset(i, "[a-z]+"))) > 0) %>% 
    arrange(
      as.numeric(str_replace(query, "Ww(\\d+)_.*", "\\1")),
      as.numeric(str_replace(query, ".*NODE_(\\d+)_.*", "\\1")),
      as.numeric(str_extract(query, "\\d+$"))
    )
  rm(ec_no_substrate, ec_substrate_curated)
}

substrates_by_ec_final %>% 
  mutate(across(where(is.list), 
                function(j) map_chr(j, paste, collapse = " "))) %>% 
  write_tsv("data/substrate_by_ec_final.tsv")

## CGC matched to PUL ----
substrates_by_pul <- match_cgc2pul(cgc, pul)
substrates_by_pul <- left_join(substrates_by_pul,
                               puldb_meta[, c("ID",
                                              "substrate_final",
                                              "degradation_biosynthesis",
                                              "organism_name")],
                               by = join_by("pul_id" == "ID")) %>%
  distinct(cluster_md5,
           pul_id,
           substrate_final,
           degradation_biosynthesis,
           organism_name)

## CGC substrates by majority voting ----
# The original dbCAN3 code had weights for substrates depending on the order of motifs found in the gene. I'm not entirely certain of the biological reason for those weightings. Furthermore, I have included additional substrates predicted from ECs assigned via other methods such as CAZyDB, KOfam, and KEGG. Therefore, I'm inclined to use equal weighting for each substrate per gene.

substrates_by_vote <- inner_join(
  cgc[, c("cluster_md5", "query")],
  substrates_by_ec_final[, c("query", "substrate")],
  by = join_by("query")
) %>%
  group_by(cluster_md5, query) %>% 
  # Substrates per gene
  summarise(
    substrate = map(list(substrate), function(i) unique(reduce(i, union)))
  ) %>% 
  # Substrates per CGC
  group_by(cluster_md5) %>% 
  summarise(
    substrate = map(list(substrate), function(i) {
      votes <- table(unlist(i))
      names(votes)[votes == max(votes)]
    })
  )

## Hierarchical substrate designation ----
tmp.substrate_1 <- inner_join(
  cgc[, c("cluster_md5", "query")],
  substrates_by_pul[, c("cluster_md5", "pul_id", "substrate_final")],
  by = join_by("cluster_md5"),
  relationship = "many-to-many"
) %>% 
  mutate(prediction_method = "CGC-PUL match") %>% 
  distinct() %>% 
  rename("substrate" = "substrate_final")

tmp.substrate_2 <- inner_join(
  cgc[, c("cluster_md5", "query")], 
  substrates_by_vote,
  by = join_by("cluster_md5"),
  relationship = "many-to-many"
) %>%
  filter(!cluster_md5 %in% tmp.substrate_1$cluster_md5) %>% 
  mutate(prediction_method = "CGC majority voting") %>% 
  distinct()

tmp.substrate_3 <- substrates_by_ec_final %>% 
  filter(!query %in% unique(c(tmp.substrate_1$query, tmp.substrate_2$query))) %>% 
  select(query, substrate) %>% 
  mutate(prediction_method = "EC-substrate match") %>% 
  distinct()

# Check for collisions
length(intersect(tmp.substrate_1$query, tmp.substrate_2$query))
length(intersect(tmp.substrate_1$query, tmp.substrate_3$query))
length(intersect(tmp.substrate_2$query, tmp.substrate_3$query))

# Combine
substrates_final <- list(tmp.substrate_1, 
                         tmp.substrate_2, 
                         tmp.substrate_3) %>% 
  map(function(df) {
    if (class(df$substrate) != "list")
      df$substrate <- map(df$substrate, function(i) unlist(as.list(i)))
    df
  }) %>% 
  bind_rows() %>% 
  arrange(
    as.numeric(str_replace(query, "Ww(\\d+)_.*", "\\1")),
    as.numeric(str_replace(query, ".*_NODE_(\\d+)_.*", "\\1")),
    as.numeric(str_extract(query, "\\d+$"))
  ) %>% 
  mutate(substrate = map(substrate, 
                         function(i) unique(str_replace_all(i, " ", "_"))))

rm(list = ls(pattern = "tmp."))

# Deduplicate ----
dedup_substrates_final <- substrates_final %>% 
  group_by(cluster_md5, query, prediction_method) %>% 
  summarise(
    pul_id = map(list(pul_id), wrap),
    putative_substrate = map(list(substrate), wrap)
  )

dedup_ec_final <- ec_final %>% 
  group_by(query) %>% 
  summarise(
    cazyme_ec = map(list(ec), wrap),
    cazyme_ec_db = map(list(ec_db), wrap)
  )

dedup_dbcan <- dbcan %>% 
  group_by(query) %>% 
  summarise(
    dbcan_family = map(list(dbcan_family), wrap),
    dbcan_query_from = map(list(query_from), wrap),
    dbcan_query_to = map(list(query_to), wrap)
  )

dedup_dbsub <- dbsub %>% 
  group_by(query) %>% 
  summarise(
    dbsub_subfamily = map(list(dbsub_family), wrap),
    dbsub_ec = map(list(dbsub_ec), wrap),
    dbsub_query_from = map(list(query_from), wrap),
    dbsub_query_to = map(list(query_to), wrap)
  )

dedup_cazy <- cazy %>% 
  group_by(qseqid) %>% 
  summarise(
    cazy_family = map(list(cazy_family), wrap),
    cazy_accession = map(list(cazy_accession), wrap),
    cazy_ec = map(list(cazy_ec), wrap),
    cazy_subject = map(list(sseqid), wrap),
    cazy_query_from = map(list(qstart), wrap),
    cazy_query_to = map(list(qend), wrap)
  ) %>% 
  rename(query = qseqid)

tmp.tcdom <- tcdom[, c("query", "tcdom_family", "query_from", "query_to")]
tmp.tcdb <- tcdb[, c("qseqid", "tcdb_family", "qstart", "qend")]
tmp.tcdb <- anti_join(tmp.tcdb, tmp.tcdom, 
                      by = join_by("qseqid" == "query", 
                                   "tcdb_family" == "tcdom_family"))
dedup_transporter <- list(
  "TCDB" = tmp.tcdb,
  "TCDom" = tmp.tcdom
) %>% 
  map(function(i) {
    colnames(i) <- c("query", "transporter_family", 
                     "transporter_query_start", "transporter_query_end")
    i
  }) %>% 
  bind_rows() %>% 
  left_join(tcdb_meta, by = "transporter_family") %>% 
  select(query, transporter_family, transporter_name, everything()) %>% 
  group_by(query) %>% 
  summarise(
    transporter_family = map(list(transporter_family), wrap),
    transporter_name = map(list(transporter_name), wrap),
    transporter_query_from = map(list(transporter_query_start), wrap),
    transporter_query_to = map(list(transporter_query_end), wrap)
  )

dedup_kegg <- kegg %>% 
  mutate(kegg_ko = strsplit(kegg_ko, "/")) %>% 
  unnest(kegg_ko) %>% 
  left_join(kegg_ko_def, by = join_by("kegg_ko" == "kegg_entry")) %>% 
  group_by(qseqid) %>% 
  summarise(
    kegg_subject = map(list(sseqid), wrap),
    kegg_ko = map(list(kegg_ko), wrap),
    kegg_symbol = map(list(kegg_symbol), wrap),
    kegg_definition = map(list(kegg_name), wrap),
    kegg_query_from = map(list(qstart), wrap),
    kegg_query_to = map(list(qend), wrap)
  ) %>% 
  rename(query = qseqid)
  
dedup_signal_transduction <- stp %>% 
  group_by(query) %>% 
  summarise(
    signal_transduction_subject = map(list(target), wrap),
    signal_transduction_query_from = map(list(query_from), wrap),
    signal_transduction_query_to = map(list(query_to), wrap)
  )

dedup_transcription_factor <- tf %>% 
  group_by(query) %>% 
  summarise(
    transcription_factor_subject = map(list(target), wrap),
    transcription_factor_query_from = map(list(query_from), wrap),
    transcription_factor_query_to = map(list(query_to), wrap)
  )

rm(list = ls(pattern = "tmp"))

## Compile deduplicated tables ----
curated_ann <- mget(ls(pattern = "dedup_.*")) %>% 
  reduce(function(i, j) full_join(i, j, by = "query")) %>% 
  left_join(ann, ., by = "query") %>% 
  # Kofam EC was made before wrap() was written. Clean up that column.
  mutate(kofam_ec = map(kofam_ec, wrap)) %>% 
  arrange(
    as.numeric(str_replace(query, "Ww(\\d+)_.*", "\\1")),
    as.numeric(str_replace(query, ".*NODE_(\\d+)_.*", "\\1")),
    as.numeric(str_extract(query, "\\d+$"))
  )

rm(list = ls(pattern = "dedup_.*"))

# Save image ----
save.image("finding_caz.RData")
save(curated_ann, file = "data/curated_annotations.RData")

# Scratch space ----

# Attempt to update EC numbers in subfam_map
# Genbank ID - EC in characterized_cazy is ground truth

tmp.chrcaz <- characterized_cazy %>% 
  filter(!is.na(genbank)) %>% 
  mutate(genbank = map(strsplit(genbank, " "), \(i) i[!is.na(i)]),
         ec_number = strsplit(ec_number, " ")) %>% 
  select(family, genbank, ec_number)


update_subfam_ec <- function(genbank_id, ec_number) {
  # inputs refer to columns in subfam_map
  # subfam_family <- gsub("([A-Z]{2,3}\\d+)_e\\d+", "\\1", subfamily)
  genbank_id <- unlist(strsplit(genbank_id, ","))
  
  matched_ec <- tmp.chrcaz %>% 
    filter(map_lgl(genbank, function(i) any(i %in% genbank_id))) %>% 
    pull(ec_number) %>% 
    unlist() %>% 
    unique()
  
  return(paste(matched_ec, collapse = " "))
  
  # if (length(matched_ec) < 1) {
  #   return("-")
  # } else {
  #   if (all(matched_ec == ec_number)) {
  #     return(ec_number)
  #   } else {
  #     paste(matched_ec, collapse = " ")
  #   }
  # }
  
}

plan(multisession, workers = availableCores())
tmp.update <- subfam_map %>% 
  mutate(
    updated_ec = future_pmap(
      list(genbank, subfamily),
      function(gb, sf, ec) update_subfam_ec(gb, sf)
    )
  )
plan(sequential)

