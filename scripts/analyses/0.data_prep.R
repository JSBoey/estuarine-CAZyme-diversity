# Prepare data

# This is run once to generate the correct .RData file for this project.
# If there are changes to the raw data files for import, then run this again.

# Data source
srctype <- c("WGS", "WTS") %>% 
  setNames(nm = .)

# Numeric data
COUNTS <- lapply(srctype, \(s) fread(glue("data/{s}_clean_count.tsv.gz")))
CONTIGCOV <- lapply(srctype, \(s) fread(glue("data/{s}.contig_coverage.tsv.gz")))
RLEN <- fread("data/read_lengths.csv.gz") %>%
  split(., .$srctype)

# Annotations
ANNOTATION <- fread("data/curated_annotation_table.tsv.gz")
CHECKM <- fread("data/hq_bins.checkm_data.txt.gz")
TAXONOMY <- fread("data/bin_taxonomy_final.txt.gz")
CAZYME <- as.data.table(readRDS("data/curated_cazymes_substrates.rds"))
METADATA <- fread("data/sample_metadata.txt")
SUMMARY <- lapply(srctype, \(s) fread(glue("data/{s}_count.tsv.summary")))
SUBSTRATE <- as.data.table(
  readxl::read_excel("data/curated_dbcansub_substrate_map_2.xlsx")
)

# Clean data ----

# ANNOTATION: Remove '_pred' suffix from bins
ANNOTATION$bin <- gsub("_pred$", "", ANNOTATION$bin)

# CONTIGCOV: Split coverage and depth tables, fix column names, and convert to
#            matrix
CONTIGCOV <- lapply(CONTIGCOV, \(dt) {
  names(dt)[1] <- "contig"
  l <- list(
    "covbases" = dt[, .SD, .SDcols = patterns("contig|covbases")],
    "numreads" = dt[, .SD, .SDcols = patterns("contig|numreads")]
  )
  l <- lapply(l, \(DT) {
    names(DT) <- gsub("W[GT]S\\.|\\.(numreads|covbases)", "", names(DT))
    return(DT)
  })
  return(l)
}) %>% 
  unlist(recursive = FALSE) %>% 
  lapply(as.matrix, rownames = "contig")

# CHECKM: Replace column name delimiters
setnames(CHECKM, \(nm) gsub(" ", "_", tolower(nm)))
setnames(CHECKM, old = "rename", new = "bin")

# COUNTS: Convert to matrix and order columns using METADATA$sample
COUNTS <- lapply(COUNTS, \(dt) {
  nm <- names(dt)[names(dt) %in% METADATA$sample]
  as.matrix(dt[, .SD, .SDcols = nm], rownames = dt$Geneid)
})

# SUMMARY: Convert to matrix and remove paths from column names
SUMMARY <- lapply(SUMMARY, \(dt) {
  m <- as.matrix(dt, rownames = "Status")
  colnames(m) <- gsub(".*\\/|\\.bam|W[GT]S\\.", "", colnames(m))
  return(m)
})

# TAXONOMY: Expand taxonomic levels and remove irrelevant columns
tax_level <- c(
  "domain", 
  "phylum", 
  "class", 
  "order", 
  "family", 
  "genus", 
  "species"
)

TAXONOMY <- TAXONOMY[
  , (tax_level) := tstrsplit(gtdb_taxonomy, ";", fill = NA)
][
  , (tax_level) := lapply(.SD, \(s) gsub("[a-z]__", "", s)), 
  .SDcols = tax_level
][
  , c("gtdb_taxonomy", "taxa") := NULL
]

# Standardised objects ----

# Habitats

habitat_colour <- c("sediment" = "#F89812", "water" = "#25A7F8")
habitat_shape = c("sediment" = 17, "water" = 19)

# Salinity: Viridis default, freshwater is purple and marine is yellow

# CAZyme colour
cazyme_colour <- c("GH" = "#CE1235", "PL" = "#59C9A5")

# Derived data ----

bin2contig <- get_bin2contig()

# Genome size estimates
CHECKM[
  , assembly_size := sapply(bin, \(x) sum(bin2contig[bin == x]$contig_size))
][
  , genome_size := (assembly_size * 100/completeness) - (assembly_size * contamination/100)
]

# Bin coverage 
# Based on Lander & Waterman (1988) doi:10.1016/0888-7543(88)90007-9
# (r * l) / g 
# r: Total reads mapped to a genome per sample
# l: Average read length
# g: Genome size
BINCOV <- map2(
  keep_at(CONTIGCOV, \(x) grepl("numreads", x)),
  get_average_read_length(),
  \(m, rlen) {
    
    DT <- as.data.table(m, keep.rownames = "contig") %>% 
      merge(bin2contig, ., all.x = TRUE, by = "contig")
    
    DT[
      , lapply(.SD, sum), # r
      .SDcols = patterns("contig_size|Sed|Filt"), by = "bin"
    ][
      , (names(rlen)) := Map(`*`, .SD, rlen), # rl
      .SDcols = names(rlen)
    ][
      , lapply(.SD, \(x) x / contig_size), # rl/g
      .SDcols = patterns("Sed|Filt"), by = "bin"
    ] %>% 
      as.matrix(rownames = "bin")
  }
)

# For relative abundance of MAGs, refer to calc_bin_relative_abundance() in 0.functions.R

# Clean up ----
rm(srctype, tax_level, bin2contig)
gc(verbose = TRUE, full = TRUE)

# Write out ----
save.image()
