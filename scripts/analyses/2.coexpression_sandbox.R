dt <- open_csv_dataset(sources = "sandbox/parcor_water")
# glimpse(dt)

dt <- open_csv_dataset("sandbox/parcor_water/")
E <- dt %>% 
  filter((node1 %in% get_GHPL() | node2 %in% get_GHPL()) & rho > 0.8 & p < .1) %>% 
  collect() %>% 
  as.data.table()

E <- E[p.adjust(p, "fdr") < 0.05]

G <- graph_from_data_frame(E[, c("node1", "node2")], directed = FALSE)
E(G)$weight <- 1/E$rho

components(G)$no

dG <- decompose(G)
dG <- dG[order(sapply(dG, \(i) length(V(i))), decreasing = T)]
# Ignore stars
is_star <- function(g) sum(degree(g) == 1) == (length(V(g)) - 1)
dG <- dG[!sapply(dG, is_star) & sapply(dG, \(i) length(V(i))) > 10]

for (i in seq_along(dG)) {
  V(dG[[i]])$articulation <- ifelse(V(dG[[i]]) %in% articulation_points(dG[[i]]), T, F)
  V(dG[[i]])$degree <- degree(dG[[i]])
  V(dG[[i]])$eigen <- eigen_centrality(dG[[i]])$vector
  V(dG[[i]])$betweenness <- betweenness(dG[[i]], normalized = TRUE)
  V(dG[[i]])$CAZyme <- names(V(dG[[i]])) %in% get_GHPL()
}

par(mfrow = c(2, 2))
for (i in 1:4) {
  color <- viridis::viridis(length(V(dG[[i]])))
  color <- color[order(V(dG[[i]])$betweenness)]
  plot(dG[[i]], layout = layout_with_fr(dG[[i]]),
       vertext.size = 3,
       vertex.color = color)
  rm(color)
}

cl.G <- cluster_fast_greedy(G)
cm.G <- membership(cl.G)


list(
  "degree_assortativity" = assortativity_degree(G, directed = FALSE),
  "transitivity" = transitivity(G),
  "mean_distance" = mean_distance(G, directed = FALSE),
  "diameter" = diameter(G, directed = FALSE),
  "edge_density" = edge_density(G),
  "number_of_components" = comp$no,
  "eigenvalue" = eigencen$value
) |>
  paste0(names(.), "\t", .) |> 
  writeLines(con = "test.txt")

ivertex <- intersect(names(V(G)), get_GHPL()) %>% 
  setNames(nm = .)
nb <- names(neighbors(G, v = ivertex))
nbG <- induced_subgraph(G, unique(c(ivertex, nb)))
plot(nbG, vertex.size = 1, vertex.label = NA, layout = layout_with_(nbG))

nb_graph <- make_neighborhood_graph(G, order = 1, nodes = ivertex)
nb_graph <- nb_graph[order(sapply(nb_graph, length), decreasing = T)]
plot(nb_graph[[27]], vertex.size = 3, vertex.label = NA)

data.frame(
  "node" = names(V(G)),
  "degree" = degree(G, normalized = FALSE),
  "strength" = strength(G),
  "betweenness" = betweenness(G, normalized = FALSE, directed = FALSE),
  "eigenvector" = eigencen$vector,
  "component_membership" = comp$membership,
  "community_membership" = membership(cluster_fast_greedy(G)),
  "articulation_point" = names(V(G)) %in% names(articulation_points(G))
) %>%
  write.csv(file = "vertex_attributes.csv", row.names = FALSE)

flist <- list(
  "node" = function(x) names(V(x)),
  "degree" = function(x) degree(x, normalized = FALSE),
  "strength" = function(x) strength(x),
  "betweenness" = function(x) betweenness(x, normalized = FALSE, directed = FALSE),
  "community_membership" = function(x) membership(cluster_fast_greedy(x)),
  "articulation_point" = function(x) names(V(x)) %in% names(articulation_points(x))
)

plan(multisession, workers = 4)
future_lapply(
  flist, \(f) f(G),
  future.packages = c("igraph")
) %>% 
  as.data.frame()

future_lapply(
  flist, \(f) system.time(f(G)),
  future.packages = c("igraph")
)


vertex_attributes_solo <- lapply(
  flist, \(f) f(G)
)

microbenchmark::microbenchmark(
  "flist_multi" = future_lapply(flist, \(f) f(G), future.packages = "igraph"),
  "flist_solo" = lapply(flist, \(f) f(G)),
  times = 10
)


# About correlations
# Negative correaltion without shared zeroes but positive with shared zeroes
ax <- c(0, 0, 6, 7, 5, 8, 9)
ay <- c(0, 0, 7, 9, 8, 6, 5)

# Weak positive w/ shared, Strong negative w/o shared
bx <- c(0, 0, 0, 1, 5, 8, 9)
by <- c(0, 0, 7, 9, 8, 6, 0)


cx <- c(0, 0, 0, 1, 5, 8, 9)
cy <- c(0, 0, 9, 7, 3, 1, 0)