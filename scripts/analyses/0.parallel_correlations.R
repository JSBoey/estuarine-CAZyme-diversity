#!/usr/bin/env Rscript

# Functions for massively parallel correlations with empirical significance 
# determined using permutations

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
                       "future.apply", 
                       "data.table", 
                       "purrr", 
                       "optparse")
suppressPackageStartupMessages(dynamicRequire(required_packages))

# Functions ----
rank_matrix <- function(mat,
                        parallel = FALSE,
                        cores = 0) {
  matrix(
    rowRanks(mat, parallel = parallel, cores = cores),
    nrow = nrow(mat),
    ncol = ncol(mat),
    dimnames = dimnames(mat)
  )
}

parallel_permcor <- function(mat, 
                             path,
                             workers = 1,
                             k = 999, 
                             p.threshold = 2,
                             rho.threshold = 2) {
  # Initialise
  if (workers > 1) {
    plan(multisession, workers = workers)
    cat("Running with multisession futures\n")
  } else {
    plan(sequential)
    cat("Running with sequential futures\n")
  }
  
  on.exit(plan(sequential))
  
  if (!dir.exists(path)) dir.create(path, recursive = TRUE)
  cat("Temporary outputs are here:\n")
  cat("  ", path)
  
  nr <- nrow(mat)
  
  # Parallel correlations
  future_lapply(seq_len(nr - 1), \(i) {
    
    # Output
    pid <- Sys.getpid()
    filename <- paste0(path, "/parallel_permcor_tmp.", pid, ".csv")
    
    # Correlations (Pearson)
    res <- lapply((i + 1):nr, \(j) {
      
      S <- permcor(mat[i, ], mat[j, ], k)
      
      if (abs(S[[1]]) < rho.threshold & S[[2]] >= p.threshold) {
        
        S <- NULL
        
      } else {
        
        c(
          "node1"  = rownames(mat)[i],
          "node2"  = rownames(mat)[j],
          "rho"    = S[[1]],
          "perm.p" = S[[2]]
        )
        
      }
      
    })
    
    # Convert list to data.table
    res <- res[sapply(res, length) > 0]
    res_dt <- as.data.table(list_transpose(res))
    
    # Progress message
    if (i %% 50 == 0) message(paste(i, "rows processed"))
    if (i == nr - 1) message("Completed")
    
    # Write output
    fwrite(res_dt, file = filename, append = TRUE)
    
  }, 
  future.seed = TRUE, 
  future.scheduling = Inf)
  
}

# Parse options ----

option_list <- list(
  make_option("--in-file", type = "character", 
              help = "A TSV or CSV file where the first column is a unique identifier"),
  make_option("--out-file", type = "character", default = "./parallel_permcor.csv",
              help = "Output file [defaults to %default]"),
  make_option("--method", type = "character", default = "pearson",
              help = "Coefficient type (pearson or spearman) [default = %default]"),
  make_option("--rho-threshold", type = "double", default = 2,
              help = "Correlations coefficients lower than this value will be filtered out. [default = %default]"),
  make_option("--p-threshold", type = "double", default = 2,
              help = "Correlations with significance greater than or equal to this value will be filtered out [default = %default]"),
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

m <- as.matrix.data.table(fread(opt$in_file), rownames = 1)

if (!is.numeric(m))
  stop("Input file must be coercible to a numeric matrix")

if (!(opt$method %in% c("pearson", "spearman")))
  stop("--method accepts only 'pearson' or 'spearman'")

## Variables
parallel <- ifelse(opt$threads > 1, TRUE, FALSE)
cores <- ifelse(opt$threads > 1, opt$threads, 0)

if (opt$method == "spearman") {
  m <- rank_matrix(m, parallel = parallel, cores = cores)
}

parallel_permcor(mat = m,
                 path = opt$tmp_dir, 
                 workers = cores,
                 k = opt$permutations, 
                 p.threshold = opt$p_threshold,
                 rho.threshold = opt$rho_threshold)

## Collate
tmp_out <- list.files(path = opt$tmp_dir, 
                      pattern = "parallel_permcor_tmp.*.csv$",
                      full.names = TRUE)

for (tmpfile in tmp_out) {
  fwrite(x = fread(tmpfile), file = opt$out_file, append = TRUE)
}

## Remove temporary files
if (isFALSE(opt$keep_tmp_dir)) 
  system(paste("rm -rf", opt$tmp_dir))
