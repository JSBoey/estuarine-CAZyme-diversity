# Function library

# This script holds all the convenience and custom functions for this project.

# Convenience functions ----
unwrap <- function(x) unique(unlist(x))

threshold_matrix <- function(x, cutoff) {
  x[x < cutoff] <- 0
  x
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

# Normalisation ----

# Formula as in Zhao, Ye, and Stanton (2020) doi:10.1261/rna.074922.120
# This is a suitable standardisation for ordination, but not for differential
# expression analysis.
calc_TPM <- function(x, g) {
  if (is.null(names(g))) {
    stop("g must be a named numeric vector")
  }
  
  g <- g[order(match(names(g), rownames(x)))] 
  
  rpk <- sweep(x, 1, g/1e3, "/")
  sweep(rpk, 2, colSums(rpk)/1e6, "/")
}




