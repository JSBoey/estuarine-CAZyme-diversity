#!/usr/bin/env Rscript
# Co-expression analyses - Graph metrics

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

libs <- c("arrow", "dplyr", "optparse", "igraph", "future.apply", "data.table")
dynamicRequire(libs)

rm(libs)

# Programme options ----
option_list <- list(
  make_option("--graph", type = "character",
              help = "Path to graphs in GraphML format."),
  make_option("--out_path", type = "character",
              help = "Output directory")
)

prog_description <- "Calculates centralities and graph topologies and find communities from input graph"

# Parse options ----
opt <- parse_args(object = OptionParser(option_list = option_list, 
                                        description = prog_description),
                  convert_hyphens_to_underscores = TRUE)

# Main ----
## Future topology
workers <- availableCores()

if (availableCores("multicore") > 1L) {
  plan(multicore, workers = workers)
} else if (workers > 1){
  plan(multisession, workers = workers)
} else {
  plan(sequential)
}

on.exit(plan(sequential))

## Read graph
G <- read_graph(opt$graph, format = "graphml")

cat(paste0("Imported graph with ", 
           vcount(G), " vertices and ", 
           ecount(G), " edges\n"))

## Shared metrics
comp <- components(G)
eigcen <- eigen_centrality(G)

cat(paste0("Graph has ", comp$no, 
           " components and eigenvalue of ", eigcen$value, "\n"))

## Global metrics
flist_global_metrics <- list(
  "degree_assortivity" = function(g) assortativity_degree(g),
  "global_transitivity" = function(g) transitivity(g, type = "global"),
  "mean_distance" = function(g) mean_distance(g),
  "diameter" = function(g) diameter(g),
  "edge_density" = function(g) edge_density(g)
)

## Centralities
flist_centralities <- list(
  "degree" = function(g) degree(g),
  "strength" = function(g) strength(g),
  "betweenness" = function(g) betweenness(g)
)

## Clustering
flist_clustering <- list(
  "infomap" = function(g) cluster_infomap(g),
  "leiden" = function(g) {
    r <- quantile(strength(g))[2] / (gorder(g) - 1)
    cluster_leiden(g, resolution_parameter = r)
  },
  "clauset_newman_moore" = function(g) cluster_fast_greedy(g)
)

## Compile
global_metrics <- future_lapply(
  flist_global_metrics, \(f) f(G),
  future.packages = "igraph",
  future.scheduling = Inf
)

global_metrics <- append(global_metrics, 
                         list("number_of_components" = comp$no,
                              "size_of_largest_component" = max(comp$csize))) 

paste0(names(global_metrics), "\t", global_metrics) %>% 
  writeLines(con = paste0(out_path, "/global_metrics.txt"))

cat(paste0("Computed global metrics\n"))

centrality <- future_lapply(flist_centralities, \(f) f(G), 
                            future.packages = "igraph",
                            future.scheduling = Inf)
centrality <- append(centrality, list("eigenvector" = eigcen$vector))
as.data.table(as.data.frame(centrality), keep.rownames = "node") %>% 
  fwrite(file = paste0(opt$out_path, "/centrality.csv"))

cat(paste0("Computed centralities\n"))

clustering <- future_lapply(flist_clustering, \(f) f(G),
                            future.packages = "igraph",
                            future.scheduling = Inf,
                            future.seed = TRUE)

sapply(clustering, membership) %>% 
  as.data.table(keep.rownames = "node") %>% 
  fwrite(file = paste0(opt$out_path, "/cluster_membership.csv"))

cat(paste0("Computed clusters\n"))
cat(paste0("Tasks complete. Outputs are in ", opt$out_path, "\n"))
