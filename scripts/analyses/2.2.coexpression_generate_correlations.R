#!/usr/bin/env Rscript

# Massively parallel correlations with empirical significance 
# determined using permutations and p-value 95% confidence intervals obtained
# using Wilson score method. 

# NOTE: Some functions are duplicated from 0.functions.R as this needs to be 
#       sourced at HPC side.

# Libraries ----
dynamicRequire <- function(libs, quietly = TRUE) {
  
  libs_to_install <- character(0)
  
  for (lib in libs) {
    if(!require(lib, character.only = TRUE, quietly = quietly)) {
      libs_to_install <- c(libs_to_install, lib)
    }
  }
  
  if (length(libs_to_install) > 0) {
    install.packages(libs_to_install)
    sapply(libs_to_install, \(lib) {
      library(lib, character.only = TRUE, quietly = quietly)
    })
  }
  
  sapply(libs_to_install, \(lib) {
    library(lib, character.only = TRUE, quietly = quietly)
  })
  
  invisible(libs)
}

required_packages <- c("Rfast", 
                       "binom",
                       "future.apply", 
                       "data.table", 
                       "purrr", 
                       "optparse")
dynamicRequire(required_packages)

# Functions ----
# Create permutations from a vector
permat <- function(x, k) replicate(k, x[sample.int(length(x), replace = FALSE)])

# Spearman correlation with significance by permutation
perm_sp_cor <- function(x, y, k) {
  # Remove shared zeroes
  i <- rowSums(cbind(x, y)) > 0
  n <- sum(i)
  x <- Rank(x[i])
  y <- Rank(y[i])
  
  # Observed correlation
  t0 <- cor(x, y)
  
  # Generate k permutations
  pY <- permat(y, k)
  
  # Generate permuted correlations
  t1 <- as.vector(cor(x, pY))
  
  # Significance
  ns <- sum(abs(t1) >= abs(t0)) 
  p_ci <- binom.confint(ns, k, conf.level = .95, methods = "wilson")
  
  # Result
  list("rho" = t0, 
       "n" = n, 
       "p" = p_ci$mean, 
       "p_lower" = p_ci$lower, 
       "p_upper" = p_ci$upper)
}

# Split outputs into chunks per thread
parcor <- function(matrix, 
                   k = 1e3,
                   workers = 1,
                   temp.path = "./tmp/") {
  # Init
  if (workers > 1) {
    plan(multisession, workers = workers)
    cat("Running with multisession futures\n")
  } else {
    plan(sequential)
    cat("Running with sequential futures\n")
  }
  on.exit(plan(sequential))
  nr <- nrow(matrix)
  
  future_lapply(
    seq_len(nr - 1),
    \(i) {
      # Open temporary outputs
      pid <- Sys.getpid()
      if (!dir.exists(temp.path)) dir.create(temp.path)
      temp.file <- paste0(temp.path, "/", pid)
      
      # Correlations
      results <- lapply((i+1):nr, \(j) {
        append(
          list("node1" = rownames(matrix)[i],
               "node2" = rownames(matrix)[j]),
          perm_sp_cor(matrix[i, ], matrix[j, ], k = k)
        )
      })
      results <- as.data.table(list_transpose(results))
      
      # Filter
      
      fwrite(results,
             file = temp.file,
             append = TRUE)
    },
    future.seed = TRUE,
    future.scheduling = Inf
  )
}


# Parse options ----
option_list <- list(
  make_option("--in-file", type = "character", 
              help = "A TSV or CSV file where the first column is a unique identifier"),
  make_option("--out-file", type = "character", default = "./parallel_permcor.csv",
              help = "Output file [defaults to %default]"),
  make_option("--threads", type = "integer", default = future::availableCores()[[1]],
              help = "Number of threads for parallel processing [default = all available cores]"),
  make_option("--permutations", type = "integer", default = 999,
              help = "Number of permutations to estimate significance [default = %default]"),
  make_option("--tmp-dir", type = "character", default = paste0(getwd(), "/tmp/parallel_permcor"),
              help = "Temporary directory [default = ./tmp/parallel_permcor]"),
  make_option("--keep-tmp-dir", type = "logical", 
              action = "store_true", default = FALSE,
              help = "Keep temporary outputs")
)

prog_description <- "Calculates row-wise pairwise correlations and evaluates significance using permutations."

opt <- parse_args(object = OptionParser(option_list = option_list, 
                                        description = prog_description),
                  convert_hyphens_to_underscores = TRUE)

# Main ----
## Validate data
if (is.null(opt$in_file))
  stop("Please supply an input file delimited by [,\t |;:]")

m <- data.table:::as.matrix.data.table(fread(opt$in_file), rownames = 1)

if (!is.numeric(m))
  stop("Input file must be coercible to a numeric matrix")

parcor(matrix = m, 
       k = opt$permutations, 
       workers = opt$threads, 
       temp.path = opt$tmp_dir)

## Collate
tmp_out <- list.files(path = opt$tmp_dir,
                      full.names = TRUE)

for (tmpfile in tmp_out) {
  fwrite(x = fread(tmpfile), file = opt$out_file, append = TRUE)
}

## Remove temporary files
if (isFALSE(opt$keep_tmp_dir)) 
  system(paste("rm -rf", opt$tmp_dir))
