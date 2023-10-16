# Parallelised correlations 2: Calculate significance by approximation

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
method <- "spearman"
a <- 0.05
nc <- availableCores()

odn <- paste0(dirname(args[2]), "/asymCorP")
if (!dir.exists(odn)) dir.create(odn, recursive = TRUE)
ofn <- paste0(odn, "/", basename(args[2]), ".asymCorP")

# Approximate significance 
cat(glue("Approximating significance for chunk: {args[2]}"), "\n")
sig <- parAsymCorP(
  x = dat, r = r, method = method, a = a, nc = nc
)
cat(glue("Exporting results to {ofn}"), "\n")
fwrite(sig, file = ofn)

cat("Completed part 2\n")

### END ###
