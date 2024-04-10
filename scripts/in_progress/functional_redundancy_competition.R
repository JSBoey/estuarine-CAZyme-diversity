# Functional redundancy in endolytic CAZymes

# Environment ----
library(data.table)
library(matrixStats)
library(tidyverse)
library(furrr)

# Import data ----
ANNOTATION <- fread("results/curated_annotation_table.tsv")
CAZYME     <- readRDS("data/curated_cazymes_substrates.rds")
METADATA   <- fread("data/sample_metadata.txt")
BINCOV     <- fread("data/bin_coverage.txt")
TPM        <- readRDS("results/transformed_counts.rds")

# Mapping files 
brenda <- brendaDb::ReadBrenda(brendaDb::DownloadBrenda())
substrate <- readxl::read_xlsx("data/curated_dbcansub_substrate_map_2.xlsx")

# Functions ----
unwrap <- function(list) {
  # Creates a unique vector from a list or list-like object
  unique(unlist(list))
}

splitHabitat <- function(m, habitat) {
  # Splits count/TPM matrix by pattern
  m[, grepl(habitat, colnames(m))]
}

getCAZyBySubstrate <- function(x, action, CAZyme_class) {
  # Output should be a list of:
  # 1. EC Numbers
  # 2. CAZyme families
  # 3. Specified action mode of the enzyme
  
  # Substrate
  cols <- grep("Substrate", names(substrate), value = TRUE)
  f1 <- substrate[
    , Reduce(`|`, lapply(.SD, \(a) grepl(x, a))),
    .SDcols = cols
  ]
  sb <- substrate[f1]
  
  # Endolytic or Exolytic
  ACTION <- c('endolytic', 'exolytic', 'all')
  action <- ACTION[pmatch(action, ACTION)]
  
  if (action == "endolytic") sb <- sb[Endolytic == "Yes"]
  if (action == "exolytic") sb <- sb[is.na(Endolytic)]
  
  # CAZyme class
  p1 <- paste0(CAZyme_class, collapse = "|")
  sb <- sb[sapply(Family, \(x) grepl(p1, x))]
  
  # Results
  list(
    "EC"     = unwrap(sb[, "EC_Number"]),
    "Family" = unwrap(sb[, "Family"]),
    "Action" = action
  )
}

findCAZymes <- function(EC) {
  # Output should be a list of:
  # 1. substrate from substrate table
  # 2. hits from CAZYME
  # 3. vector of latent bins split by habitat
  # 4. vector of active bins split by habitat
  # 5. matrix of genes with capacity split by habitat
  # 6. matrix of transcripts split by habitat
  
  # EC can be a result from getCAZyBySubstrate or a vector of EC Numbers
  if (is.list(EC)) {
    ec <- EC$EC
    fm <- EC$Family
  } else {
    ec <- EC
  }
  
  pool <- unwrap(c(substrate$EC_Number, CAZYME$EC))
  if (!all(ec %in% pool)) 
    stop("EC numbers do not exist in mapping files. Please check.")
  
  # Subset substrates
  dat1 <- substrate[EC_Number %chin% ec]
  
  # Subset hit list
  cols <- c('bin', 'node', 'dbcan_label', 'tcdb_label', 'sulfatlas_label',
            'signalp_prediction', 'transcription_regulation')
  dat2 <- CAZYME[sapply(EC, \(x) any(x %in% ec))]
  mags <- ANNOTATION[node %chin% dat2$NODE, ..cols]
  dat2 <- dat2[mags, on = c("NODE" = "node")]
  
  # Additional Family-based subset as required
  if (exists("fm")) {
    dat1 <- dat1[Family %chin% fm]
    dat2 <- dat2[FAMILY %chin% fm]
  }
  
  # Gene and transcripts by habitat
  tpm <- lapply(TPM[grepl("TPM", names(TPM))], \(m) {
    m <- m[unwrap(dat2$NODE), ]
    # Define habitats
    h <- c("water" = "Filt", "sediment" = "Sed")
    lapply(h, \(i) splitHabitat(m, i))
  })
  
  # Latent and active bins
  lab <- rapply(tpm, \(m) {
    rn <- rownames(m[rowSums(m) > 0, ])
    unwrap(ANNOTATION[node %in% rn, "bin"])
  }, how = 'replace')
  
  # Output
  list(
    "substrates"  = dat1,
    "hits"        = dat2,
    "bins_latent" = lab$wgs_TPM,
    "bins_active" = lab$wts_TPM,
    "genes"       = tpm$wgs_TPM,
    "transcript"  = tpm$wts_TPM
  )
}

aggregateStats <- function(m) {
  flist <- list(
    "total" = sum,
    "average" = \(x) mean(x[x > 0]),
    "stdev" = \(x) sd(x[x > 0]),
    "min" = \(x) min(x[x > 0]),
    "max" = max,
    "entropy" = vegan::diversity,
    "hits" = \(x) sum(x > 0)
  )
  X <- as.data.frame(m)
  as.data.frame(lapply(flist, \(f) sapply(X, f)))
}
# Clean data ----
# Convert substrate and CAZYME to data.table
lapply(c("CAZYME", "substrate"), \(x) setDT(get(x)))

# Define endolytic CAZymes ----
endolytic <- unique(
  substrate[Endolytic == "Yes"], 
  by = c("Substrate_high_level", "Family", "EC_Number")
)

# Substrates of interest ----
SOI <- list(
  'cellulose'       = findCAZymes('3.2.1.4'),
  'starch/glycogen' = findCAZymes('3.2.1.1'),
  'xylan'           = findCAZymes('3.2.1.8'),
  'beta-mannan'     = findCAZymes('3.2.1.78'),
  'laminarin'       = findCAZymes(getCAZyBySubstrate('laminarin', 'endo', c('GH', 'PL'))),
  'pectin'          = findCAZymes(getCAZyBySubstrate('pectin', 'endo', c('GH', 'PL'))),
  'chitin'          = findCAZymes(getCAZyBySubstrate('chitin', 'endo', c('GH', 'PL'))),
  'chitosan'        = findCAZymes(getCAZyBySubstrate('chitosan', 'endo', c('GH', 'PL'))),
  'alginate'        = findCAZymes(getCAZyBySubstrate('alginate', 'endo', c('GH', 'PL'))),
  'ulvan'           = findCAZymes(getCAZyBySubstrate('ulvan', 'endo', c('GH', 'PL')))
)

HELPER <- c('cellulose', 'chitin', 'pectin', 'xylan')
names(HELPER) <- HELPER
HELPER <- lapply(HELPER, \(x) {
  findCAZymes(getCAZyBySubstrate(x, 'all', c('CE', 'AA')))
})

# Active v. latent bins
actVLat <- lapply(
  set_names(c("bins_latent", "bins_active"), \(x) str_remove(x, "bins_")),
  \(header) sapply(SOI, \(l) length(unwrap(l[[header]])))
) %>% 
  as.data.frame() %>% 
  mutate(ratio = round(active/latent, 3))

lapply(SOI, \(l) aggregateStats(l$transcript$water))

# Terminal enzymes ----
terminal_enzymes <- c(
  "β-glucosidase" = "3.2.1.21", 
  "α-glucosidase" = "3.2.1.20",
  "β-galactosidase" = "3.2.1.23",
  "α-mannosidase" = "3.2.1.24",
  "β-mannosidase" = "3.2.1.25",
  "β-xylosidase" = "3.2.1.37",
  "β-D-fucosidase" = "3.2.1.38",
  "α-L-rhamnosidase" = "3.2.1.40",
  "α-L-fucosidase" = "3.2.1.51",
  "β-N-acetylhexosaminidase" = "3.2.1.52",
  "α-L-arabinofuranosidase" = "3.2.1.55"
)
TERMINAL <- lapply(terminal_enzymes, findCAZymes)

lapply(
  set_names(c("bins_latent", "bins_active"), \(x) str_remove(x, "bins_")),
  \(header) sapply(TERMINAL, \(l) length(unwrap(l[[header]])))
) %>% 
  as.data.frame() %>% 
  mutate(ratio = round(active/latent, 3))

lapply(TERMINAL, \(l) aggregateStats(l$transcript$water))

# Latency and expression per substrate per MAG ----
latentExpressedBins <- lapply(SOI, \(l) {
  # Extract only node and bin from hits table
  h <- unique(l$hits[, c("NODE", "bin")])
  
  # Sum of latent and expressed per bin
  lapply(l[c("genes", "transcript")], \(ll) {
    lapply(ll, \(m) {
      # Data table with distinct nodes
      M <- unique(as.data.table(m, keep.rownames = "NODE"))
      # Join bin information for summarised data
      S <- M[h, on = "NODE"][, NODE := NULL][, lapply(.SD, sum), by = "bin"]
      Snum <- names(S)[sapply(S, is.numeric)]
      S[rowSums(S[, ..Snum]) > 0]
    })
  }) %>% 
    transpose() %>% 
    lapply(., \(l3) bind_rows(l3, .id = "hit_type")[order(bin)])
}) %>% 
  transpose() %>% 
  lapply(., \(l4) {
    x <- bind_rows(l4, .id = "substrate")
    setorderv(x, c("bin", "substrate", "hit_type"))
  })

bin_expression <- lapply(SOI, \(l) {
  # Extract only node and bin from hits table
  h <- unique(l$hits[, c("NODE", "bin")])
  
  # Create sum per bin
  lapply(l$transcript, \(m) {
    # Remember to create data.table with only distinct rows
    M <- unique(as.data.table(m, keep.rownames = "NODE"))
    # Join bin information for summarised data
    S <- M[h, on = "NODE"][, NODE := NULL][, lapply(.SD, sum), by = "bin"]
    Snum <- names(S)[sapply(S, is.numeric)]
    S[rowSums(S[, ..Snum]) > 0]
  })
}) %>% 
  transpose() %>% 
  lapply(., \(l) bind_rows(l, .id = "substrate")[order(bin)])

## Correlation with multiple substrates
cross_SOI <- lapply(bin_expression, \(dat) {
  # Split by substrate
  by_substrate <- split(dat, dat$bin) %>% 
    keep(~ nrow(.x) > 1)
  
  return(by_substrate)
})

# Distribution ----


# Save objects ----
save.image("rdata/functional_redundancy_competition.RData")

# Scratch space ----
test <- lapply(latentExpressedBins, \(l) {
  a <- l[hit_type == "genes", .N, by = "bin"]
  b <- l[hit_type == "transcript", .N, by = "bin"]
  b[a, on = "bin"] %>% 
    setnames(., c("N", "i.N"), c("active", "latent"))
})

lapply(SOI, \(l) {
  l[grepl("bins", names(l))] %>% 
    unlist(recursive = F) %>% 
    lapply(., length) %>% 
    as.data.frame()
}) %>% 
  rbindlist(use.names = TRUE, idcol = "substrate")

test_2 <- lapply(latentExpressedBins, \(l) {
  l2 <- lapply(c("genes", "transcript") %>% set_names(., .), \(ht) {
    d1 <- l[
      hit_type == ht, c("substrate", "bin")
    ][
      , value := 1
    ]
  })
  
  d2 <- Reduce(rbind, l2)[
    , lapply(.SD, sum), .SDcols = "value", by = c("substrate", "bin")
  ]
  d3 <- dcast(d2, bin ~ substrate, value.var = "value", fill = 0)
  d4 <- d3[
    , expressed := rowSums(.SD == 2), .SDcols = alginate:xylan
  ][
    , latent := rowSums(.SD == 1), .SDcols = alginate:xylan
  ]
  
  d4
})


# Redundant (review to remove) ----