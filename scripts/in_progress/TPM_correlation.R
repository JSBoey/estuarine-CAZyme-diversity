# Initial correlation analyses with TPM normalised counts

# Environment ----
library(edgeR)
library(DESeq2)
library(tidyverse)

source("scripts/in_progress/normalise_counts.R")

# Import data ----
wgs <- read_tsv("results/WGS_clean_count.tsv.gz")
wts <- read_tsv("results/WTS_clean_count.tsv.gz")

# Clean data ----
gene_length <- wts$End - wts$Start + 1
names(gene_length) <- wts$Geneid

# Split tables ----
count_data <- list(
  "water" = select(wts, Geneid, starts_with("Filt")),
  "sediment" = select(wts, Geneid, starts_with("Sed"))
) %>% 
  map(\(x) as.matrix(column_to_rownames(x, "Geneid")))

# Create sediment sum table ----
sediment_patterns <- colnames(count_data$sediment) %>% 
  str_remove("_\\d+$") %>% 
  unique()

count_data$sediment_sum <- matrix(
  nrow = nrow(count_data$sediment), 
  ncol = length(sediment_patterns), 
  dimnames = list(rownames(count_data$sediment), sediment_patterns)
)

for (i in seq_along(sediment_patterns)) {
  
  pattern <- sediment_patterns[i]
  count_data$sediment_sum[, pattern] <- rowSums(
    count_data$sediment[, grep(pattern, colnames(count_data$sediment))]
  )
  
  rm(i, pattern)
  
}

# Normalise data ----
tpm_data <- map(count_data, \(x) norm_TPM(x, gene_length = gene_length))

# Data 1: Remove non-expressed ----
tpm_data <- map(tpm_data, \(x) x[rowSums(x) > 0, ])

# Data 2: Expressed in >= a 1/3 of samples ----
one_third_data <- map(tpm_data, \(x) x[rowSums(x > 0) >= ncol(x)/3, ])

# Data 3: Singletons ----
singletons <- map(tpm_data, \(x) x[rowSums(x > 0) == 1, ])

# Gene-wise correlation with filtering criteria ----
calc_comparisons <- function(x) {
  n <- nrow(x)
  (n * (n - 1))/2
}

test_data <- one_third_data$water[c(1:24, 8975:9e3), ]

test_truth <- broom::tidy(Hmisc::rcorr(t(test_data), type = "spearman")) %>% 
  filter(p.value < 0.05) %>% 
  select(-n) %>% 
  dplyr::rename(
    "gene1" = "column1",
    "gene2" = "column2",
    "p_value" = "p.value"
  )

gene_cor <- function(x, method = "spearman", sig_cutoff = 0.05, 
                     cor_cutoff = NULL) {
  
  if (!(is.numeric(x) & is.matrix(x))) {
    stop("Argument x must be numeric matrix")
  }
  
  # Maximum number of comparisons
  M <- calc_comparisons(x)
  
  # List index counter for correlation that meet threshold
  a <- 1
  
  # Counter for each iteration
  b <- 1
  
  # Number of genes
  n <- nrow(x)
  
  # Print messages
  cat(str_glue("Matrix has {n} genes and {M} unique comparisons"), "\n")
  
  # Initiate output
  output <- list("X" = NA_character_, 
                 "Y" = NA_character_, 
                 "estimate" = NA_real_, 
                 "p_value" = NA_real_)
  
  # Outer loop: 1 to nrow(x)
  for (i in seq_len(n)) {
    # Condition to stop if i reaches end of matrix or counter reaches all valid 
    # comparisons
    if (i == n || b > M) {
      break
    }
    
    # Inner loop: i+1 to nrow(x)
    for (j in seq(i+1, n)) {
      
      xi <- x[i, ]
      xj <- x[j, ]
      nm_xi <- rownames(x)[i]
      nm_xj <- rownames(x)[j]
      
      # Correlation test
      test <- cor.test(xi, xj, method = method)
      
      # Add to list if valid
      if (test$p.value < 0.05) {
        output$X[a] <- nm_xi
        output$Y[a] <- nm_xj
        output$estimate[a] <- test$estimate
        output$p_value[a] <- test$p.value
        a <- a + 1
      }
      
      b <- b + 1
      
      if (b %% 100 == 0) {
        cat(str_glue("Processeed {b} valid comparisons"), "\n")
      }
      
    }
    
  }
  
  cat(str_glue("Found {(a-1)} associations that pass threshold"), "\n")
  
  as.data.frame(output)
  
}

gene_cor_parallel <- function(x, method = "spearman", sig_cutoff = 0.05, 
                              cor_cutoff = NULL, ncores = NULL, np = NULL) {
  
  if (!(is.numeric(x) & is.matrix(x))) {
    stop("Argument x must be numeric matrix")
  }
  
  # if (!(is.numeric(permute) | is.null(permute))) {
  #   stop("Argument 'permute' should be a numeric value of the number of permutations required or NULL")
  # }
  
  n <- nrow(x) # Number of genes
  M <- (n * (n - 1)) / 2 # Maximum number of unique combinations
  
  cat(str_glue("Matrix has {n} genes and {M} unique comparisons"), "\n")
  
  # Register parallel cores
  require(foreach)
  require(parallel)
  require(doParallel)
  
  if (is.null(ncores)) ncores <- detectCores() - 2
  cl <- makeCluster(ncores, 'PSOCK')
  registerDoParallel(cl)
  
  output <- 
    foreach(i = 1:(n-1), .packages = c('stats'), .combine = 'rbind') %:%
    foreach(j = (i+1):n, .combine = 'rbind') %dopar% {
      
      cr <- cor.test(x[i, ], x[j, ], method = method)
      
      if (is.numeric(np)) {
        # Permute coefficient
        cr_perm <- sapply(1:np, \(k) cor(x[i, ], sample(x[j, ])))
        # Significance
        cr$p.value <- (sum(abs(cr_perm) >= abs(cr$estimate)) + 1) / (np + 1)
      }
      
      if (cr$p.value < sig_cutoff) {
        data.frame(
          "gene1" = rownames(x)[i],
          "gene2" = rownames(x)[j],
          "estimate" = cr$estimate,
          "p_value" = cr$p.value
        )
      }
      
    }
  
  stopCluster(cl)
  
  cat(str_glue("Found {(nrow(output))} comparisons that pass thresholds"), "\n")
  
  output
  
}

test_cor <- gene_cor(test_data)
test_cor_perm <- gene_cor_parallel(test_data, np = 599)

test_cor_perm_2 <- map(1:5, \(x) gene_cor_parallel(test_data, np = 1999))

