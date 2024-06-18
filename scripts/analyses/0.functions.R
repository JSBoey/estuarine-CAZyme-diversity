# Function library

# This script holds all the convenience and custom functions for this project.

# Convenience functions ----
unwrap <- function(x) unique(unlist(x))

threshold_matrix <- function(x, cutoff) {
  x[x < cutoff] <- 0
  x
}

make_pattern <- function(x, bounded = FALSE) {
  if (isTRUE(bounded))
    x <- gsub("^|$", "\\\\b", x, perl = TRUE)
  
  paste(x, collapse = "|")
}

# If matching tables downstream, make sure to check row name order if 
# remove.zero.rowsum = TRUE
split_matrix <- function(x,
                         columns = NULL,
                         rows = NULL,
                         remove_zero_rowsum = FALSE) {
  if (!is.matrix(x)) {
    stop("x must be a matrix")
  }
  
  if (!is.numeric(x) & isTRUE(remove_zero_rowsum)) {
    stop("x must be a numeric matrix if remove_zero_rowsum = TRUE")
  }
  
  if (is.null(columns)) columns <- colnames(x)
  if (is.null(rows)) rows <- rownames(x)
  
  if (!is.null(columns) & length(columns) == 1) {
    columns <- grep(columns, colnames(x), value = TRUE)
  }
  
  if (!is.null(rows) & length(rows) == 1) {
    rows <- grep(rows, rownames(x), value = TRUE)
  }
  
  x <- x[rows, columns]
  
  if (remove_zero_rowsum) {
    x <- x[rowSums(x) > 0, ]
  }
  
  x
}

lapply_at_depth <- function(LIST, DEPTH, FUN) {
  if (DEPTH == 1) {
    lapply(LIST, FUN)
  } else {
    lapply(LIST, function(l) lapply_at_depth(l, DEPTH - 1, FUN))
  }
}

# Convenience extractors ----
get_GHPL <- function() unwrap(CAZYME$NODE[grepl("(GT|PL)", CAZYME$FAMILY)])

get_contig_length <- function(x) {
  as.numeric(gsub(".*length_(\\d+)_.*", "\\1", x))
}

get_bin2contig <- function() {
  ANNOTATION[
    , .(bin, node)
  ][
    , `:=`(
      contig = gsub("_\\d+$", "", node),
      contig_size = get_contig_length(node)
    )
  ][
    , node := NULL
  ] %>% 
    unique() %>% 
    na.omit()
}

get_average_read_length <- function() {
  lapply(RLEN, \(dt) {
    dt <- dt[
      , .(x = sum(as.numeric(read_length) * count) / sum(count)), 
      by = sample
    ]
    x <- setNames(dt$x, dt$sample)
    x
  })
}

get_procrustes_axes <- function(x) {
  rot <- x$rotation
  rot_m <- numeric(2)
  for (i in seq_len(nrow(rot))) {
    rot_m[i] <- rot[i, 2]/rot[i, 1]
  }
  data.frame(
    intercept = numeric(2),
    slope = rot_m
  )
}

get_gene_length <- function() {
  g <- ANNOTATION %$% (end - start + 1)
  names(g) <- ANNOTATION$node
  
  g
}

get_breadth_coverage <- function(x) {
  DT <- merge(
    get_bin2contig(),
    as.data.table(x, keep.rownames = "contig"),
    by = "contig",
    all.x = TRUE
  )
  cols <- colnames(x)
  DT[
    , lapply(.SD, sum),
    by = "bin", .SDcols = c("contig_size", cols)
  ][
    , (cols) := Map(\(i) i / contig_size, .SD),
    .SDcols = cols
  ][
    , .SD, .SDcols = c("bin", cols)
  ] %>% 
    as.matrix(rownames = "bin")
}

get_EC_by_substrate <- function(substrate, 
                                enzyme_class = NULL, 
                                action = c("all", "endolytic", "exolytic")) {
  # Output:
  # 1. EC numbers
  # 2. CAZyme families
  # 3. Mode of action
  
  cols <- grep("substrate",
               names(SUBSTRATE),
               value = TRUE,
               ignore.case = TRUE)
  pattern <- paste0(substrate, collapse = "|")
  substrate_match <- SUBSTRATE[ # Output: Boolean OR across multiple columns
    , Reduce(`|`, lapply(.SD, \(j) grepl(pattern, j))),
    .SDcols = cols
  ]
  
  pattern_2 <- paste0(enzyme_class, collapse = "|")
  family_match <- grepl(pattern_2, SUBSTRATE$Family)
  
  DT <- SUBSTRATE[substrate_match & family_match]
  
  action <- match.arg(action)
  
  DT <- if (action == "endolytic") {
    DT[Endolytic == "Yes"]
  } else if (action == "exolytic") {
    DT[is.na(Endolytic)]
  } else {
    DT
  }
  
  as.list(DT[, c("EC_Number", "Family")])
  
}

get_cazyme_nodes <- function(ec,
                             additional_column = NULL) {
  # Parse output from get_EC_from_substrate
  if (is.list(ec)) {
    ec_number <- make_pattern(ec$EC_Number, bounded = TRUE)
    cazy_family <- ec$Family
  } else {
    ec_number <- ec
  }
  
  DT <- CAZYME[grepl(ec_number, EC)]
  
  if (exists("cazy_family")) 
    DT[FAMILY %in% cazy_family]
  
  if (!all(additional_column %in% names(ANNOTATION)))
    stop("Please pick a column from ANNOTATION")
  
  cols <- c("bin", "node", additional_column)
  
  merge(
    ANNOTATION[, ..cols], 
    DT[, c("NODE", "SUBFAMILY", "EC")],
    all.y = TRUE, 
    by.x = "node", 
    by.y = "NODE"
  )
  
}

# Normalisation ----

# Transcripts per million
# Formula as in Zhao, Ye, and Stanton (2020) doi:10.1261/rna.074922.120
# This is a suitable standardisation for ordination, but not for differential
# expression analysis.
calc_TPM <- function(x, g, metric = c("tpm", "rpk")) {
  if (is.null(names(g))) {
    stop("g must be a named numeric vector")
  }
  
  g <- g[order(match(names(g), rownames(x)))] 
  
  rpk <- sweep(x, 1, g/1e3, "/")
  
  if (metric == "rpk") return(rpk)
  
  sweep(rpk, 2, colSums(rpk)/1e6, "/")
}

# MAG relative abundance (read coverage) estimated based on 
# Probst et al. (2018) doi:10.1038/s41564-017-0098-y
# Ar = (Nm / Ns) * (r * l) / g
# Ar: Relative abundance of a genome
# Nm: Maximum number of reads of all metagenome samples
# Ns: Total reads for a sample
# r: Total reads mapped to a genome per sample
# l: Average read length
# g: Genome size
# Here, coverage should be a matrix of (r * l) / g

calc_bin_relative_abundance <- function(coverage,
                                        community = c("bin", "all")) {
  community <- match.arg(community)
  
  if (community == "bin") {
    Ns <- colSums(CONTIGCOV$WGS.numreads)
  }
  
  if (community == "all") {
    Ns_dt <- RLEN$WGS[
      , .(depth = sum(as.numeric(read_length) * as.numeric(count))),
      by = "sample"
    ]
    Ns <- Ns_dt %$% setNames(depth, sample)
  }
  
  sweep(coverage, 2, max(Ns)/Ns, "*")
}

# MAG relative activity
# Mickol et al. (2021) doi:10.1128/AEM.01676-21
# Ac = (r / Nr) / (g / Ng)
# r: Transcript reads mapped to each sample per MAG
# Nr: Total transcript reads mapped to each sample
# g: Genome size
# Ng: Total size of all genomes

calc_bin_relative_activity <- function(WTS_reads) {
  
  DT <- as.data.table(WTS_reads, keep.rownames = "contig") %>% 
    merge(get_bin2contig(), ., all.x = TRUE, by = "contig")
  
  samples <- colnames(WTS_reads)
  
  DT[
    , lapply(.SD, sum), # per bin mapped reads and contig size
    .SDcols = c(samples, "contig_size"), by = bin
  ][
    , (c(samples, "contig_size")) := Map(proportions, .SD), # r/Nr; g/Ng
    .SDcols = c(samples, "contig_size")
  ][
    , lapply(.SD, \(j) j / contig_size),
    .SDcols = samples, by = bin
  ] %>% 
    as.matrix(., rownames = "bin")

}

# General per group proportions/normalisation
# For WTS per MAG, this represents the MAG's "investment" in the function
# For WTS normalised by MAG abundance, this is the MAG's contribution to ecosystem function relative to the community.
# Thanks Kim Handley and Mike Hoggard for the concepts

calc_per_bin_proportions <- function(matrix,
                                     table,
                                     join_column,
                                     sum_column,
                                     normalise = NULL) {
  
  DT <- as.data.table(matrix, keep.rownames = join_column) %>% 
    merge(table, ., all.x = TRUE, by = join_column)
  
  cols <- colnames(matrix)
  
  if (is.null(normalise)) {
    
    f <- function(x) {
      pp <- proportions(x)
      pp[is.na(pp)] <- 0
      pp
    }
    
    proportion_matrix <- DT[
      , (cols) := Map(f, .SD),
      .SDcols = cols, by = by_key,
      env = list(
        by_key = sum_column
      )
    ] 
    
    # Remove summation grouping
    # proportion_matrix <- proportion_matrix[
    #   , .SD, .SDcols = c(join_column, colnames(matrix))
    # ] %>% 
    #   as.matrix(rownames = join_column)
    
    return(proportion_matrix)
  }
  
  normalise <- normalise[, cols]
  
  matrix_list <- split(DT, by = sum_column, keep.by = FALSE) %>% 
    lapply(as.matrix, rownames = join_column)
  
  for (i in seq_along(matrix_list)) {
    A <- sweep(matrix_list[[i]], 
               MARGIN = 2, 
               STATS = normalise[names(matrix_list)[i], ],
               FUN = "/")
    A[is.na(A) | is.infinite(A)] <- 0
    matrix_list[[i]] <- A
  }
  
  Reduce("rbind", matrix_list)
  
}

## Unused code
## For reference

# getCAZyBySubstrate <- function(x, action, CAZyme_class) {
#   # Output should be a list of:
#   # 1. EC Numbers
#   # 2. CAZyme families
#   # 3. Specified action mode of the enzyme
#   
#   # Substrate
#   cols <- grep("Substrate", names(SUBSTRATE), value = TRUE)
#   f1 <- SUBSTRATE[
#     , Reduce(`|`, lapply(.SD, \(a) grepl(x, a))),
#     .SDcols = cols
#   ]
#   sb <- SUBSTRATE[f1]
#   
#   # Endolytic or Exolytic
#   ACTION <- c('endolytic', 'exolytic', 'all')
#   action <- ACTION[pmatch(action, ACTION)]
#   
#   if (action == "endolytic") sb <- sb[Endolytic == "Yes"]
#   if (action == "exolytic") sb <- sb[is.na(Endolytic)]
#   
#   # CAZyme class
#   p1 <- paste0(CAZyme_class, collapse = "|")
#   sb <- sb[sapply(Family, \(x) grepl(p1, x))]
#   
#   # Results
#   list(
#     "EC"     = unwrap(sb[, "EC_Number"]),
#     "Family" = unwrap(sb[, "Family"]),
#     "Action" = action
#   )
# }
# 
# findCAZymes <- function(EC) {
#   # Output should be a list of:
#   # 1. substrate from substrate table
#   # 2. hits from CAZYME
#   # 3. vector of latent bins split by habitat
#   # 4. vector of active bins split by habitat
#   # 5. matrix of genes with capacity split by habitat
#   # 6. matrix of transcripts split by habitat
#   
#   # EC can be a result from getCAZyBySubstrate or a vector of EC Numbers
#   if (is.list(EC)) {
#     ec <- EC$EC
#     fm <- EC$Family
#   } else {
#     ec <- EC
#   }
#   
#   pool <- unwrap(c(SUBSTRATE$EC_Number, CAZYME$EC))
#   if (!all(ec %in% pool)) 
#     stop("EC numbers do not exist in mapping files. Please check.")
#   
#   # Subset substrates
#   dat1 <- SUBSTRATE[EC_Number %chin% ec]
#   
#   # Subset hit list
#   cols <- c('bin', 'node', 'dbcan_label', 'tcdb_label', 'sulfatlas_label',
#             'signalp_prediction', 'transcription_regulation')
#   dat2 <- CAZYME[sapply(EC, \(x) any(x %in% ec))]
#   mags <- ANNOTATION[node %chin% dat2$NODE, ..cols]
#   dat2 <- dat2[mags, on = c("NODE" = "node")]
#   
#   # Additional Family-based subset as required
#   if (exists("fm")) {
#     dat1 <- dat1[Family %chin% fm]
#     dat2 <- dat2[FAMILY %chin% fm]
#   }
#   
#   # Gene and transcripts by habitat
#   dat3 <- lapply(PARTITION, \(l) l[["tpm"]]) %>%
#     lapplyAtDepth(2, \(m) {
#       DT <- as.data.table(m, keep.rownames = "NODE")[dat2, on = "NODE"]
#       f <- DT[, rowSums(!is.na(.SD)) > 0, .SDcols = patterns("^Sed|^Filt")]
#       DT[f]
#     })
#   
#   # Output
#   return(
#     list(
#       "substrates" = dat1,
#       "WGS"        = lapply(dat3, \(l) l[["WGS"]]),
#       "WTS"        = lapply(dat3, \(l) l[["WTS"]])
#     )
#   )
#   gc()
# }


