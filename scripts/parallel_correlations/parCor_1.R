# Parallelised correlations 1: Calculate pairwise correlations

# Get functions
source("parCor_functions.R")

# Read inputs
# args[1]: numeric matrix
# args[2]: prefix for output (add "/" to point to a directory)
args <- commandArgs(trailingOnly=TRUE)
dat <- as.matrix(
  fread(file = args[1], data.table = TRUE),
  rownames = TRUE
)

# Variables
k <- nrow(dat)/5e3
rho <- 0.6
method <- "spearman"
nc <- availableCores()

ofn <- sub("\\.tsv$", "", basename(args[1]))

# Calculate correlations ----
dat_cor <- parCor(x = dat, k = k, rho = rho, method = method, nc = nc)
fwrite(dat_cor, file = paste0(args[2], ofn, ".cor.tsv"), sep = "\t")

## Create temporary directory ----
tmp_dir <- paste0("tmp/parCor/", ofn)
dir.create(tmp_dir, recursive = TRUE)

## Chunk correlations into temporary files ----
cat("Chunking outputs into 1,000,000 rows/file for parCor_2.R\n")
lapply(
  chunkData(dat_cor, k = nrow(dat_cor)/1e6), function(l) {
    chkname <- paste0(tmp_dir, "/", digest::digest(l))
    fwrite(l, chkname)
  }
)
cat("Completed part 1\n")

### END ###
