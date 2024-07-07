#!/usr/bin/env Rscript
# Coexpression analyses - Build and decompose graph

# Description ----
# Builds a global network based on hard threshold similarities (as edge list)
# then decomposes them into connected components. This allows for processing of
# each component individually.

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

libs <- c("arrow", "dplyr", "optparse", "igraph", "future.apply")
dynamicRequire(libs)

rm(libs)

# Programme options ----
option_list <- list(
  make_option("--dataset", type = "character",
              help = "Parquet/Arrow dataset directory path"),
  make_option("--out_path", type = "character",
              help = "Output directory"),
  make_option("--weighted", type = "logical",
              action = "store_true", default = FALSE,
              help = "Create weighted graph"),
  make_option("--correlation_cutoff", type = "double",
              help = "Cutoff for correlation coefficient (hard threshold)"),
  make_option("--significance_cutoff", type = "double",
              help = "Cutoff for P-value (hard threshold)")
)

prog_description <- "Builds neighbourhood graphs from Parquet/Arrow dataset"

# Parse options ----
opt <- parse_args(object = OptionParser(option_list = option_list, 
                                        description = prog_description),
                  convert_hyphens_to_underscores = TRUE)

# Main ----
if (!dir.exists(opt$out_path))
  dir.create(opt$out_path, recursive=TRUE)

## Build graph
if (cpu_count() != availableCores()) 
  set_cpu_count(availableCores())

DT <- open_dataset(opt$dataset)

cat(paste0("Imported ", opt$dataset, " containing ", nrow(DT), " rows\n"))
cat(paste0("Filtering data for correlations >= ", opt$correlation_cutoff,
           " and p-values < ", opt$significance_cutoff, "\n"))

E <- DT %>%
  filter(rho >= opt$correlation_cutoff & p < opt$significance_cutoff) %>%
  collect()

cat(paste0("Retained ", nrow(E), " after filtering\n"))

G <- graph_from_data_frame(E[, c("node1", "node2")], 
                           directed = FALSE)

cat(paste0("Built global graph with ",
           length(V(G)), " vertices and ",
           length(E(G)), " edges\n"))

if (isTRUE(opt$weighted)) {
  z <- min(abs(E$rho)[abs(E$rho) > 0])
  w <- z + sqrt(1 - E$rho)
  G <- set_edge_attr(G, "weight", value = w)
  cat(paste0("Edge weights set with constant ", w, "\n"))
}

## Decompose
cm <- components(G)[c("csize", "no")]
cmf <- as.data.frame(table(cm$csize))
colnames(cmf) <- c("Size", "Frequency")
cmf$Size <- as.numeric(as.character(cmf$Size))
cat(paste0("Found ", cm$no, " components with distributions:\n"))
print(cmf)

dG <- decompose(G, min.vertices = 50)
dG <- dG[order(sapply(dG, length), decreasing = TRUE)]
cat(paste0("Extracted ", length(dG), " components with at least 50 vertices\n"))

sG <- G - igraph:::union.igraph(dG)
sG <- delete.vertices(sG, degree(sG) == 0)
cat(paste0("Extracted a subgraph with ",
           length(V(sG)), " vertices and ",
           length(E(sG)), " edges",
           " consisting of ", components(sG)$no, 
           " components with less than 50 vertice containing\n"))

## Write graphs
write_graph(G, paste0(opt$out_path, "/global.graphml"), format = "graphml")
cat(paste0("Global graph written to ", opt$out_path, "\n"))

for (i in seq_along(dG)) {
  g <- dG[[i]]
  filename <- paste0(opt$out_path, "/component", i, 
                     "_v", length(V(g)),
                     "e", length(E(g)))
  write_graph(g, filename, format = "graphml")
  cat(paste0("Component ", i, " written\n"))
}

write_graph(sG,
            paste0(opt$out_path,
                   "/subgraphs_v50_g",
                   components(sG)$no,
                   ".graphml"),
            format = "graphml")
cat(paste0("Subgraph written\nDone!\n"))

