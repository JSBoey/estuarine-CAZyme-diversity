# Parallelised correlation with the option of permutation test

# Description ----
# Performs pairwise correlations with "chunking" and parallel processing
# Can filter correlations by coefficient magnitude and/or significance
# Significance can be estimated by permutation with stopping criteria using "simctest" (Gandy, 2009)
# Confidence intervals for significance and permutation stops can be estimated

# Environment ----
pkgs <- c(
  "stats", "simctest", "snow", "doSNOW", "foreach", "parallel", "glue", "furrr",
  "progressr", "tidyverse"
)
lapply(
  pkgs, \(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x)
    }
  }
)

# Enumerate combinations ----
ncomb <- function(x) {
  if (!(is.matrix(x) || is.list(x))) stop("x must be a numeric matrix or list")
  n <- if (is.matrix(x)) nrow(x) else length(x)
  nCk <- prod(n, n-1, 0.5)
  return(nCk)
}

# Chunk data ----
chunkData <- function(x, k) {
  
  n <- nrow(x)
  row_id <- seq_len(n)
  breaks <- pretty(row_id, k)
  f <- cut(row_id, breaks = breaks)
  l_id <- split(row_id, f)
  l <- lapply(l_id, \(i) x[i, ])
  
  cat(glue("Data split into {length(l)} chunks"), "\n")
  
  return(l)
  
}

# Within matrix correlation with filtering ----
# 1. Correlation matrix within chunks
# 2. Only retain pairs that have coef >= rho
selfCor <- function(x, rho, method) {
  
  y <- cor(t(x), method = method)
  
  y_ut <- y[upper.tri(y, diag = FALSE)]
  ind <- which(upper.tri(y, diag = FALSE), arr.ind = TRUE)
  
  res <- data.frame(
    "gene1" = rownames(y)[ind[, 1]],
    "gene2" = colnames(y)[ind[, 2]],
    "estimate" = y_ut
  )
  
  res <- res[abs(res$estimate) >= rho, ]
  
  return(res)
  
}

# Between-matrix correlations with filtering ----
crossCor <- function(x1, x2, rho, method) {
  
  y <- cor(t(x1), t(x2), method = method)
  
  res <- expand.grid(
    "gene1" = rownames(y), "gene2" = colnames(y), 
    stringsAsFactors = FALSE
  )
  res$estimate <- as.vector(y)
  
  res <- res[abs(res$estimate) >= rho, ]
  
  return(res)
  
}

# Chunked combinatorial correlation with filtering ----
combCor <- function(l, rho = 0.3, method = "spearman",
                    num.cores = NULL) {
  
  # Register multi-processor back-end
  num.cores <- if (is.null(num.cores)) detectCores() - 1 else num.cores
  cl <- makeCluster(num.cores)
  registerDoSNOW(cl)
  
  cat(glue("Registered {num.cores} parallel processes"), "\n")
  
  # Variables
  k <- length(l) # Number of chunks
  nCk <- ncomb(l) # Number of pairwise chunks to compare
  
  # Within chunks
  cat(glue("Calculating {k} within-chunk correlations"), "\n")
  pb <- txtProgressBar(max = k, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  
  self <- foreach(u = 1:k, .combine = 'bind_rows', .export = "selfCor", 
                  .options.snow = opts) %dopar% {
    selfCor(l[[u]], rho = rho, method = method)
  }
  
  # Between chunks
  cat("\n")
  cat(glue("Calculating {nCk} between-chunk correlations"), "\n")
  pb <- txtProgressBar(max = nCk, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  
  btwn <- foreach(u = 1:(k-1), .combine = 'bind_rows', .export = "crossCor",
                  .options.snow = opts) %:% 
    foreach(v = (u+1):k, .combine = 'bind_rows') %dopar% {
      crossCor(l[[u]], l[[v]], rho = rho, method = method)
    }
  
  cat("\n")
  cat("Closing parallel processors", "\n")
  stopCluster(cl)
  
  cat("Row binding data frame", "\n")
  res <- dplyr::bind_rows(self, btwn)
  
  cat(glue("Complete, {nrow(res)} valid correlations found"), "\n")
  
  gc()
  
  return(res)
  
}

# Estimate significance by approximation with filtering ----
testCorP <- function(x, r, k, method, a, num.cores = NULL) {
  
  # Chunk data
  l <- chunkData(r, k)
  
  # Plan futures
  if (is.null(num.cores)) num.cores <- detectCores() - 1
  plan(multisession, workers = num.cores)
  
  # Progress updates
  pb <- progressor(steps = length(l))
  
  # Parallel P-value estimation by chunk
  p <- future_map(l, \(u) {
    pb()
    X <- map(u$gene1, \(ni) x[rownames(x) == ni, ])
    Y <- map(u$gene2, \(nj) x[rownames(x) == nj, ])
    map2_dbl(X, Y, \(i, j) cor.test(i, j, method = method)$p.value)
  })
  
  res <- mutate(r, "p" = unlist(p)) %>% 
    filter(p < a)
  
  gc()
  plan(sequential)
  
  return(res)
  
}

# Empirical P-values using permutations
simCorP <- function(x, r, k, method, a, num.perm, num.cores = NULL, 
                    level, epsilon = 1e-3, CI_level = 0.95) {
  
  # Chunk data
  l <- chunkData(r, k)
  
  # Plan futures
  if (is.null(num.cores)) num.cores <- detectCores() - 1
  plan(multisession, workers = num.cores)
  
  cat(glue("Registered {num.cores} processes", "\n"))
  
  # Progress updates
  pb <- progressor(steps = length(l))
  
  # Parallel P-value estimation by permutation with simctest
  cat("Calculating empirical P-values")
  fopt <- furrr_options(seed = TRUE)
  pl <- future_map(l, .options = fopt, \(u) {
    
    pb()
    
    X <- map(u$gene1, \(ni) x[rownames(x) == ni, ])
    Y <- map(u$gene2, \(nj) x[rownames(x) == nj, ])
    
    pXY <- pmap(list(X, Y, u$estimate), \(i, j, rho) {
      
      # Generator function
      gen <- function() abs(cor(i, sample(j), method = method)) >= rho
      
      # Permute p-values
      simP <- simctest(gen, level = level, epsilon = epsilon, maxsteps = num.perm)
      bounds <- if (is.na(simP@p.value)) getbounds(simP) else c(NA, NA)
      ci <- confint(simP, level = CI_level)
      
      # Results
      sim_res <- data.frame(
        "p" = simP@pos/simP@steps,
        "permutations" = simP@steps,
        "lower_bound" = bounds[1],
        "upper_bound" = bounds[2],
        "lower_CI" = ci[1],
        "upper_CI" = ci[2]
      )
      
      return(sim_res)
      
    })
    
    sig <- bind_cols(u, bind_rows(pXY))
    return(sig)
    
  })
  
  res <- bind_rows(pl) %>% 
    filter(p < a)
  
  gc()
  plan(sequential)
  
  return(res)
}

# Test code ----
library(tidyverse)
load("parallel_cor.test_data.RData")

medium_size <- sample(nrow(one_third_data$water), size = 500)
medium_data <- one_third_data$water[medium_size, ]

test_medium <- chunkData(medium_data, 6) %>% 
  combCor(rho = 0.6)

with_progress({
  test_medium_approx <- testCorP(
    x = medium_data, r = test_medium, k = nrow(test_medium)/2e3, 
    method = "spearman", a = 1
    )
})

system.time(
with_progress({
  test_medium_sim <- simCorP(
    x = medium_data, r = test_medium, k = nrow(test_medium)/2e3, method = "spearman", a = 2, 
    num.perm = 1e5, level = 0.05, epsilon = 1e-5, CI_level = 0.95
  )
})
)

test_medium_sig[sample(nrow(test_medium_sig), 1), ]
test_X <- medium_data[rownames(medium_data) == "NODE_326_length_75820_cov_47.484749_56", ]
test_Y <- medium_data[rownames(medium_data) == "NODE_2_length_742964_cov_26.862788_58", ]
cor.test(test_X, test_Y, method = "spearman")$p.value


plan(sequential)
gc()
