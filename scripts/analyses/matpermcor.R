# Packages
library(Rfast)
library(future.apply)
library(data.table)
library(purrr)

# Set variables
workers <- availableCores()
args <- commandArgs(trailingOnly = TRUE)

# Import data
matrix_list <- readRDS(args[1])
cat(paste("Imported", args[1]), "\n")

# Function
corPermMat2 <- function(mat, path, workers = 1,
                        k = 999, p.threshold = 2, rho.threshold = 2) {
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
  cat("Outputs are here:\n")
  cat("  ", path)
  
  nr <- nrow(mat)
  
  # Parallel correlations
  future_lapply(seq_len(nr - 1), \(i) {
    
    # Output
    pid <- Sys.getpid()
    filename <- paste0(path, "/corPermMat2_output.", pid, ".csv")
    
    # Correlations
    res <- lapply((i + 1):nr, \(j) {
      
      S <- permcor(mat[i, ], mat[j, ], k)
      
      if (abs(S[[1]]) < rho.threshold & S[[2]] > p.threshold) {
        
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
    if (i %% 50 == 0) message(paste(i, "transcripts processed"))
    if (i == nr - 1) message("Completed")
    
    # Write output
    fwrite(res_dt, file = filename, append = TRUE)
    
  }, future.seed = TRUE, future.scheduling = Inf)
  
}

# Run
for (n in seq_len(length(matrix_list))) {
  path <- paste0(getwd(), "/", names(matrix_list)[n])
  corPermMat2(matrix_list[[n]], path = path,
              workers = workers, p.threshold = 0.05)
}

