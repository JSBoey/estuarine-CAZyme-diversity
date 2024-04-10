# # Parallelised correlation with the option of permutation test

# Description ----
# Performs pairwise correlations with "chunking" and parallel processing
# Can filter correlations by coefficient magnitude and/or significance
# Significance can be estimated by permutation with stopping criteria using "simctest" (Gandy, 2009)
# Confidence intervals for significance and permutation stops can be estimated

# Style and package choices ----
# Base R for all code within functions (avoids dependencies)
# future parallel backends via doFuture invoked by foreach API

# Packages ----
pkgs <- c("simctest", "glue", "data.table", "doFuture", "digest")
dynamicRequire <- function(pkg) {
  if (!require(pkg, character.only = TRUE)) install.packages(pkg)
  library(pkg, character.only = TRUE)
}
sapply(pkgs, dynamicRequire)

# Enumerate combinations ----
# x : numeric matrix
ncomb <- function(x) {
  n <- nrow(x)
  prod(n, n-1, 0.5) # Number of rows
}

# Chunk data ----
# k : number of chunks
chunkData <- function(x, k) {
  rowID   <- seq_len(nrow(x)) # Row indices
  breaks  <- pretty(rowID, k) # Breaks points 
  f       <- cut(rowID, breaks) # Break by
  lID     <- split(rowID, f) # Split rows by f
  
  cat(glue("Split data into {length(lID)} chunks"), "\n")
  
  lapply(lID, function(i) x[i, ]) # List of chunks
}

# Within-matrix correlation with filtering ----
# rho : correlation coefficient filtering threshold
selfCor <- function(x, rho, method) {
  y <- cor(t(x), method = method)
  
  y_ut <- y[upper.tri(y, diag = FALSE)]
  ind <- which(upper.tri(y, diag = FALSE), arr.ind = TRUE)
  
  res <- data.table(
    "gene1" = rownames(y)[ind[, 1]],
    "gene2" = colnames(y)[ind[, 2]],
    "estimate" = y_ut
  )
  
  res[abs(estimate) >= rho]
}

# Between-matrix correlations wtih filtering ----
# x1, x2 : numeric matrices 1 and 2
xCor <- function(x1, x2, rho, method) {
  y <- cor(t(x1), t(x2), method = method)
  
  res <- expand.grid("gene1" = rownames(y), "gene2" = colnames(y), 
                     stringsAsFactors = FALSE)
  
  as.data.table(res)[, estimate := as.vector(y)][abs(estimate) >= rho]
}

# Combinatorial chunked correlation with filtering ----
## Depends on:
## - chunkData()
## - selfCor()
## - xCor()
# nc : number of parallel processes
parCor <- function(x, k, rho, method, nc = NULL) {
  cat(glue("Matrix has {nrow(x)} genes and {ncomb(x)} pairwise correlations"), "\n")
  
  # Chunk data
  l <- chunkData(x, k)
  
  # If num.cores = NULL, use all available cores
  nc <- if (is.null(nc)) availableCores() else nc
  if (Sys.info()[1] == "Linux") {
    plan(multicore, workers = nc)
  } else {
    plan(multisession, workers = nc)
  }
  
  
  cat(glue("Registered {nc} parallel processes"), "\n")
  
  K <- length(l)
  
  # Within-chunk correlation
  cat(glue("Calculating correlations for {K} matrices"), "\n")
  uu <- foreach(u = l, .combine = 'rbind') %dofuture% {
    selfCor(u, rho = rho, method = method)
  }
  
  # Between-chunk correlation
  nCk <- prod(K, K-1, 0.5) # Number of pairwise matrices
  cat(glue("Calculating correlations for {nCk} matrix pairs"), "\n")
  uv <- foreach(u = 1:(K-1), .combine = 'rbind') %:%
    foreach(v = (u+1):K, .combine = 'rbind') %dofuture% {
      xCor(l[[u]], l[[v]], rho = rho, method = method)
    }
  
  cat("Closing parallel processes\n")
  plan(sequential)
  
  cat("Binding rows\n")
  rbind(uu, uv)
}

# Correlation p-values by asymptotic approximation ----
# r : Result from parCor()
# ksize : Number of rows per chunk to process
# a : Significance threshold where P(rho != 0) < a
parAsymCorP <- function(x, r, ksize, method, a = 0.05, nc = NULL) {
  # If num.cores = NULL, use all available cores
  nc <- if (is.null(nc)) availableCores() else nc
  if (Sys.info()[1] == "Linux") {
    plan(multicore, workers = nc)
  } else {
    plan(multisession, workers = nc)
  }
  
  cat(glue("Registered {nc} parallel processes"), "\n")
  
  # Approximate asymptotic significance
  cat(glue("Calculating {nrow(r)} P-values"), "\n")
  opts.future <- list(chunk.size = ksize)
  asymP <- foreach(
    u = r$gene1, v = r$gene2, 
    .options.future = opts.future) %dofuture% {
      cor.test(
        x[rownames(x) == u, ], x[rownames(x) == v, ], method = method
      )$p.value
    }
  
  cat("Closing parallel processes\n")
  plan(sequential)
  
  cat("Combining and filtering results\n")
  as.data.table(r)[, p := unlist(asymP)][p < a]
  
}

# Correlation p-values by permutation ----
# level : Significance threshold for simctest
# epsilon : Resampling risk
# CI_level : Confidence interval level
parPermCorP <- function(x, r, ksize, method, a = 0.05, nc = NULL, np = 1e5,
                        epsilon = 1e5, CI_level = 0.95) {
  
  if (a >= 1 || a <= 0) stop("a must be between 0 and 1")
  
  # If num.cores = NULL, use all available cores
  nc <- if (is.null(nc)) availableCores() else nc
  if (Sys.info()[1] == "Linux") {
    plan(multicore, workers = nc)
  } else {
    plan(multisession, workers = nc)
  }
  
  cat(glue("Registered {nc} parallel processes"), "\n")
  
  # Permute P-values
  cat(glue("Calculating {nrow(r)} P-values"), "\n")
  opts.future <- list(chunk.size = ksize, seed = TRUE)
  simP <- foreach(
    u = r$gene1, v = r$gene2, rho = r$estimate,
    .combine = 'rbind', .options.future = opts.future
  ) %dofuture% {
    # Extract variables
    X <- x[rownames(x) == u, ]
    Y <- x[rownames(x) == v, ]
    
    # Generator function
    gen <- function() abs(cor(X, sample(Y), method = method)) >= abs(rho)
    
    # Permutations
    pXY <- simctest(gen, a, epsilon, np)
    p <- if (is.na(pXY@p.value)) pXY@pos/pXY@steps else pXY@p.value
    bounds <- if (is.na(pXY@p.value)) getbounds(pXY) else c(NA, NA)
    ci <- confint(pXY, level = CI_level)
    
    # Results
    data.table(
      "p" = p,
      "permutations" = pXY@steps,
      "lower_bound" = bounds[1],
      "upper_bound" = bounds[2],
      "lower_CI" = ci[1],
      "upper_CI" = ci[2]
    )
  }
  
  cat("Closing parallel processes\n")
  plan(sequential)
  
  cat("Combining and filtering results\n")
  cbind(r, simP)[p < a]
  
}

# Deploy ----
## Read arguments ----
args = commandArgs(trailingOnly=TRUE)

## Import data ----
dat <- as.matrix(
  fread(file = args[1], data.table = TRUE),
  rownames = TRUE
)

## Set variables ----
k <- nrow(dat)/500     # Should have about 500 rows per chunk
rho <- 0.3             # Lenient coefficient threshold for permutation
a <- 0.05              # Significance threshold
method <- "spearman"   # Spearman correlation
nc <- availableCores() # Use all available cores
ksize <- 500           # Number of pairwise correlations to process per chunk
np <- 1e5              # Number of permutations
epsilon <- 1e-5        # Resampling risk
CI_level <- 0.95       # 95% confidence interval

## Calculate correlations ----
dat_cor <- parCor(x = dat, k = k, rho = rho, method = method, nc = nc)
fwrite(dat_cor, file = paste0(args[2], ".cor.tsv"), sep = "\t")

## Create temporary directory ----
tmp_dir <- paste0("tmp_", Sys.getenv("SLURM_JOB_ID"))
dir.create(tmp_dir)

## Chunk correlations into temporary files ----
lapply(
  chunkData(dat_cor, k = nrow(dat_cor)/1e5), function(l) {
    chkname <- paste0(tmp_dir, "/", digest::digest(l))
    fwrite(l, chkname)
  }
)
rm(dat_cor)
gc()

## Calculate significance per chunk ----
# Read in chunked files
chkcor <- list.files(path = tmp_dir, full.names = TRUE)

lapply(
  chkcor, function(fn) {
    icor <- fread(file = fn)
    # Asymptotic p-values
    asymp <- parAsymCorP(
      x = dat, r = icor, ksize = ksize, method = method, a = a, nc = nc
    )
    fwrite(asymp, file = paste0(fn, ".asymCorP.csv"))
    # Permutation p-values
    permp <- parPermCorP(
      x = dat, r = icor, ksize = ksize, method = method, a = a, nc = nc, 
      np = np, epsilon = epsilon, CI_level = CI_level
    )
    # Write out chunked results
    fwrite(permp, file = paste0(fn, ".permCorP.csv"))
  }
)

## Combine chunks and write out results ----
lapply(
  c(".asymCorP.csv", ".permCorP.csv"), function(ext) {
    files <- list.files(path = tmp_dir, pattern = ext, full.names = TRUE)
    ofn <- paste0(args[2], ext)
    fwrite(
      rbindlist(
        lapply(files, fread, data.table = TRUE)
      ), 
      file = ofn, sep = "\t"
    )
  }
)

# 
# 
# dat_corAsymP <- parAsymCorP(
#   x = dat, r = dat_cor, ksize = ksize, method = method, a = a, nc = nc
# )
# fwrite(dat_corAsymP, file = paste0(args[2], ".cor_AsymP.tsv"), sep = "\t")
# 
# dat_corPermP <- parPermCorP(
#   x = dat, r = dat_cor, ksize = ksize, method = method, a = a, nc = nc,
#   np = np, epsilon = epsilon, CI_level = CI_level
# )
# fwrite(dat_corAsymP, file = paste0(args[2], ".cor_PermP.tsv"), sep = "\t")

### END ###
