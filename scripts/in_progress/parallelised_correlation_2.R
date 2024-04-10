# Parallelised correlation with the option of permutation test

# Description ----
# This script performs pairwise correlations for a large combinations by parallel processing. It can be tuned to retain correlations that meet a coefficient or significance threshold. It can also perform permutation test of significance using the "simctest" package to provide confident intervals for empirical p-values. This is useful for cases of many rank ties in Spearman correlation.

# Performs pairwise correlations with parallel processing.
# Can filter correlations by coefficient magnitude and/or significance
# Significance can be estimated by permutation with stopping criteria using "simctest" (Gandy, 2009)
# Confidence intervals for significance and permutation stops can be estimated

# Traditional permutation ----
resampleCorP <- function(x, y, method, num.permutation, t.stat) {
  T.stat <- replicate(num.permutation, cor(x, sample(y), method = method))
  p.hat <- mean(abs(T.stat) >= abs(t.stat))
  return(p.hat)
}

# Wilson score CI ----
calcWilsonCI <- function(p, n, level) {
  a <- 1 - level
  z <- -qnorm(a/2)
  z2_n <- (z^2)/n
  pq_n <- (p * (1 - p))/n
  
  denom <- 1 + z2_n
  top.left <- p + 0.5 * z2_n
  top.right <- z * sqrt(pq_n + (0.25 * z2_n))
  
  ci <- (top.left + c(-1, 1) * top.right) / denom
  ci[1] <- if (ci[1] < 0) 0 else ci[1]
  label <- paste(c(a/2, 1 - a/2) * 100, "%")
  ci <- matrix(ci, nrow = 1, ncol = 2, dimnames = list("p.value", label))
  
  return(ci)
}

# Number of valid combinations ----
nComb <- function(n) (n * (n - 1)) / 2

# Main function ----
parGeneCor <- function(x, method = "spearman",
                       sig.threshold = 0.05, cor.threshold = 0.3,
                       num.cores = NULL, confint.level = 0.95,
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
            .options.snow = opts, .export = c("resampleCorP", "calcWilsonCI")) %:%
    foreach(j = (i+1):n, .combine = 'rbind') %dopar% {
      
      # Extract variables
      X <- x[i, ]
      Y <- x[j, ]
      nX <- rownames(x)[i]
      nY <- rownames(x)[j]
      
      # Compute correlation
      t.stat <- cor.test(X, Y, method = method)
      
      # Stop condition based on coefficient
      if (abs(t.stat$estimate) < cor.threshold) {
        res <- NULL
      } else {
        
        res <- data.frame(
          "gene1" = nX, "gene2" = nY, "estimate" = t.stat$estimate
        )
        
        # Approximate p-value without permutations
        if (perm.method == "none") {
          p <- t.stat$p.value
          
          outdata <- data.frame(
            "p" = p, 
            "permutations" = "none"
          )
        }
        
        # Traditional permutation
        if (perm.method == "traditional") {
          p <- resampleCorP(X, Y, method, num.permutation, t.stat$estimate)
          p.ci <- calcWilsonCI(p, num.permutation, confint.level)
          
          outdata <- data.frame(
            "p" = p,
            "p.lowerCI" = p.ci[1],
            "p.upperCI" = p.ci[2],
            "permutations" = num.permutation
          )
        }
        
        # SIMCTEST permutation
        if (perm.method == "simctest") {
          
          gen <- function() {
            # Inherit everything from function environment
            rho <- cor(X, sample(Y), method = method)
            return(abs(rho) >= abs(t.stat$estimate))
          }
          
          T.stat <- 
            simctest(gen, level = sig.threshold, maxsteps = num.permutation)
          
          p <- T.stat@pos/T.stat@steps
          p.ci <- confint(T.stat, level = confint.level)
          p.bound <- if (is.na(T.stat@p.value)) getbounds(T.stat) else c(NA, NA)
          
          outdata <- data.frame(
            "p" = p,
            "p.lowerCI" = p.ci[1],
            "p.upperCI" = p.ci[2],
            "p.lowerBound" = p.bound[1],
            "p.upperBound" = p.bound[2],
            "permutations" = num.permutation
          )
        }
        
        # Construct output
        if (p >= sig.threshold) {
          res <- NULL
        } else {
          res <- cbind(res, outdata)
        }
        
      }
      
      # Final output
      return(res)
      
    }
  
  close(pb)
  stopCluster(cl)
  
  return(output)
  
}

# Testing ----
load("testing_parallalised_cor.RData")

no_perm <- parGeneCor(one_third_data$water[1:500, ], permutation.method = "sim")
