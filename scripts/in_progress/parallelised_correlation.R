# Parallelised correlation with the option of permutation test

# Description ----
# This script performs pairwise correlations for a large combinations by parallel processing. It can be tuned to retain correlations that meet a coefficient or significance threshold. It can also perform permutation test of significance using the "simctest" package to provide confident intervals for empirical p-values. This is useful for cases of many rank ties in Spearman correlation.

# Performs pairwise correlations with parallel processing.
# Can filter correlations by coefficient magnitude and/or significance
# Significance can be estimated by permutation with stopping criteria using "simctest" (Gandy, 2009)
# Confidence intervals for significance and permutation stops can be estimated

# Number of valid combinations ----
nComb <- function(n) (n * (n - 1)) / 2

# Main function ----
parGeneCor <- function(x, method = "spearman",
                       sig.threshold = 0.05, cor.threshold = 0.3,
                       num.cores = NULL,
                       permutation.method = "simctest", 
                       num.permutation = 1e5) {
  # Libraries
  req.packages <- c(
    "simctest", "glue", "stats", "foreach", "parallel", "doParallel", "doSNOW"
  )
  lapply(req.packages, function(x) {
    if (!require(x, character.only = TRUE)) install.packages(x)
    require(x, character.only = TRUE)
  })
  
  # Input error
  if (!(is.numeric(x) && is.matrix(x))) stop("x must be a numeric matrix")
  
  # Warning
  
  # Variables
  n <- nrow(x) # Number of genes
  n.comp <- nComb(n) # Maximum pairwise comparisons
  PERM.METHOD <- c("none", "traditional", "simctest")
  perm.method <- match.arg(permutation.method, PERM.METHOD)
  
  # Create cluster
  num.cores <- if (is.null(num.cores)) detectCores() - 1 else num.cores
  cl <- makeCluster(num.cores)
  registerDoSNOW(cl)
  
  # Message
  cat(
    glue(
      "
      Matrix has {n} genes with {n.comp} pairwise comparisons
      Running {num.cores} processes in parallel
      Calculating {method} correlations with thresholds:
      \t|Correlation coefficient| >= {cor.threshold}
      \tP-value < {sig.threshold}
      "
    ), "\n"
  )
  
  if (perm.method == "none") {
    cat(
      glue(
        "No permutations requested"
      ), "\n"
    )
  }
  
  if (perm.method %in% c("traditional", "simctest")) {
    cat(
      glue(
        "
        Correlation significance estimated using {perm.method} method with \\
        {num.permutation} permutations
        "
      ), "\n"
    )
  }
  
  # Progress report
  pb <- txtProgressBar(max = n.comp, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  
  # Loop
  output <- 
    foreach(i = 1:(n-1), .combine = 'rbind', .packages = req.packages[1:3], 
            .options.snow = opts) %:%
      foreach(j = (i+1):n, .combine = 'rbind') %dopar% {
        
        # Extract variables
        X <- x[i, ]
        Y <- x[j, ]
        nX <- rownames(x)[i]
        nY <- rownames(x)[j]
        
        # Compute correlation
        t.stat <- cor.test(X, Y, method = method)
        
        if (perm.method == "none") p <- t.stat$p.value
        
        # Traditional permutation
        if (perm.method == "traditional") {
          T.stat <- replicate(
            num.permutation, cor(X, sample(Y), method = method)
          )
          p <- mean(abs(T.stat) >= abs(t.stat$estimate))
          # Standard error and (Wald) 95% confidence intervals
          z <- qnorm(1 - (0.05 / 2))
          p.se <- sqrt(p * (1 - p) / num.permutation)
          p.lCI <- p - (z * p.se)
          p.uCI <- p + (z * p.se)
        }
        
        # simctest permutation
        if (perm.method == "simctest") {
          gen <- function() {
            perm.rho <- cor(X, sample(Y), method = method)
            return(abs(perm.rho) >= abs(t.stat$estimate))
          }
          T.stat <- simctest(
            gen, level = sig.threshold, maxsteps = num.permutation
          )
          
          # Getting p-values and related statistics
          if (is.na(T.stat@p.value)) {
            p <- T.stat@pos / T.stat@steps
            p.lB <- getbounds(T.stat)[1]
            p.uB <- getbounds(T.stat)[2]
          } else {
            p <- T.stat@p.value
            p.lB <- NA
            p.uB <- NA
          }
          
          p.lCI <- confint(T.stat)[1]
          p.uCI <- confint(T.stat)[2]
        }
        
        # Result data frame
        if (p < sig.threshold & abs(t.stat$estimate) >= cor.threshold) {
          res <- data.frame("gene1" = nX,
                            "gene2" = nY,
                            "coef" = t.stat$estimate,
                            "sig" = p)
          
          if (perm.method == "traditional") {
            res$p.SE <- p.se
            res$p.lowerCI <- p.lCI
            res$p.upperCI <- p.uCI
          }
          
          if (perm.method == "simctest") {
            res$p.lowerbound <- p.lB
            res$p.upperbound <- p.uB
            res$p.lowerCI <- p.lCI
            res$p.upperCI <- p.uCI
          }
        } else {
          res <- NULL
        }
        
        # Final output
        return(res)
        
      }
  
  close(pb)
  stopCluster(cl)
  
  return(output)
  
}


# Testing ----
test_cor <- parGeneCor(test_data, permutation.method = "none", num.cores = 4)
test_cor <- parGeneCor(test_data[1:5,], permutation.method = "trad", num.cores = 4)
system.time(
test_cor <- parGeneCor(test_data, permutation.method = "sim", num.cores = 4)
)

test_wts_sediment_sum_none <- 
  parGeneCor(one_third_data$sediment_sum[1:1000, ], permutation.method = "none")


