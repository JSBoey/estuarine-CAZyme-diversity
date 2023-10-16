# Parallelised correlations 3: Calculate significance by permutation

# Get functions
source("parCor_functions.R")

# Read inputs
# args[1]: numeric matrix (must be the same one in parCor_1.R)
# args[2]: temporary output from parCor_1.R
args <- commandArgs(trailingOnly=TRUE)
dat <- as.matrix(
  fread(file = args[1], data.table = TRUE),
  rownames = TRUE
)
r <- fread(file = args[2], data.table = TRUE)

# Set variables
ksize <- 1000
method <- "spearman"
a <- 0.05
nc <- availableCores()
np <- 1e4
epsilon <- 1e-5
CI_level <- 0.95

odn <- paste0(dirname(args[2]), "/permCorP")
if (!dir.exists(odn)) dir.create(odn, recursive = TRUE)
ofn <- paste0(odn, "/", basename(args[2]), ".permCorP")

# Approximate significance 
cat(glue("Permuting significance for chunk: {args[2]}"), "\n")
sig <- parPermCorP(
  x = dat, r = r, ksize = ksize, method = method, a = a, nc = nc, 
  np = np, epsilon = epsilon, CI_level = CI_level
)
cat(glue("Exporting results to {ofn}"), "\n")
fwrite(sig, file = ofn)

cat("Completed part 3\n")

### END ###
