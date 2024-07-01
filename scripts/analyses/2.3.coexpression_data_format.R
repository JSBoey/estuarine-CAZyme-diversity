#!/usr/bin/env Rscript
# 2.2. Co-expression analysis - Format transformation
# The data stored as compressed CSV are still difficult to use. Convert to 
# Parquet for accessibility.

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

libs <- c("arrow", "dplyr", "optparse")
dynamicRequire(libs)

# Parse options ----
option_list <- list(
  make_option("--in_file", type = "character",
              help = "CSV file, first row is header"),
  make_option("--parquet_path", type = "character",
              help = "Output directory for Parquet dataset"),
  make_option("--cpus", type = "integer", default = 8,
              help = "Number of CPUs [default: %default]")
)

prog_description <- "Converts CSV file into Parquet."

opt <- parse_args(object = OptionParser(option_list = option_list, 
                                        description = prog_description),
                  convert_hyphens_to_underscores = TRUE)

# Main ----
dts <- open_csv_dataset(sources = opt$in_file)
nr <- nrow(dts)

max_rows_per_file <- ifelse(nr/1000 <= 5e6, 5e6, pretty(nr/1000)[2])
set_cpu_count(opt$cpus)

write_dataset(dts,
              opt$parquet_path, 
              max_rows_per_file = max_rows_per_file)
