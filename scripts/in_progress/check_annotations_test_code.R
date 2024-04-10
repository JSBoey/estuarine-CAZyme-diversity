str_subset(unique(cazy$target), "\\|[1-7]\\.") %>% 
  length()

str_detect(unique(cazy$target), "[1-7]\\.") %>% 
  length()

overlaps <- group_by(dbcansub, query) %>% 
  filter(n() > 1) %>% 
  arrange(query, query_start, query_end) %>% 
  select(query, target, query_length, query_start, query_end) %>% 
  mutate(
    overlap = round((query_start - lag(query_end)) / query_length, 2)
  ) %>% 
  filter(if_any(overlap, ~ . < 0))

dbcansub_substrate %>% 
  filter(Family == "GH23" & EC_Number == "3.2.1.17")

setdiff(
  str_subset(unique(str_remove(DBCAN$family, "_.*")), "GT", negate = T),
  unique(dbcansub_substrate$Family)
)

str_subset(cazy_map$ec_number, " ")

test_data <- tibble(
  a = 1:10,
  b = LETTERS[1:10]
)
test_data$c <- map(1:nrow(test_data), ~ paste(test_data$a[.x], test_data$b[.x]))

dash_cazy_map <- filter(cazy_map, map_lgl(ec_number, ~ any(str_detect(.x, "-"))))
dash_dbcansub_substrate <- filter(dbcansub_substrate, str_detect(EC_Number, "-"))

# Impossible case
a <- filter(DBCANSUB, query == "hellp")$ec %>% 
  unwrap()
getHitStatus(intersect(a, poolEC("GH24")))

# No EC case
b <- filter(DBCANSUB, query == "NODE_10072_length_12146_cov_6.879168_6")$ec %>% 
  unwrap()
getHitStatus(intersect(b, a))

# EC hit case
filter(DBCANSUB, query == "NODE_100680_length_4117_cov_6.015516_4")$ec %>% 
  unwrap() %>% 
  getHitStatus()

# Multi EC hit case
yay <- filter(DBCANSUB, query == "NODE_100837_length_3505_cov_4.034279_2")$ec %>% 
  unwrap()



# This works
test_intersect <- full_join(DBCAN, DBCANSUB, by = c('query', 'family')) %>% 
  drop_na(subfamily, cluster)

plan(multisession, workers = 6)
test_intersect$ec_overlap <- future_map(1:nrow(test_intersect), \(i) {
  ec_pool <- poolEC(test_intersect$subfamily[[i]])
  intersect(test_intersect$ec[[i]], ec_pool)
})
future:::ClusterRegistry('stop')

test_overlap <- slice_sample(DBCAN, prop = 0.01)

plan(multisession, workers = availableCores()/2)

microbenchmark::microbenchmark(
  {
    test_overlap$ec <- future_map(
      1:nrow(test_overlap), \(i) {
        a <- overlapEC(
          node = test_overlap$query[[i]], 
          family = test_overlap$family[[i]], 
          subfamily = test_overlap$subfamily[[i]]
        )
        if (all(is.na(a))) a <- getKOEC(node = test_overlap$query[[i]])
        if (is.null(a)) a <- NA_character_
        return(a)
      }
    )
  }, times = 20
)

plan(sequential)

future:::ClusterRegistry('stop')


# These nodes are CAZymes without CLUSTERS
no_cluster_cazyme <- setdiff(unique(DBCAN$query), unique(DBCANSUB$query))

# These nodes are CAZymes without best EC CAZy BLAST hits
no_best_ec_cazyme <- setdiff(unique(DBCAN$query), unique(CAZY$best_ec$query))

# These nodes have no CLUSTERS and best EC CAZy BLAST hits
no_ec <- intersect(no_best_ec_cazyme, no_cluster_cazyme)

filter(DBCAN, query %in% no_ec) %>% View()

# Try strategies for parallel processing
test <- slice_sample(DBCAN, n = 100)

plan(multisession, workers = 6)

microbenchmark::microbenchmark(
  {
    test$ec1 <- future_pmap(
      list(test$query, test$family, test$subfamily), 
      ~ assignEC(..1, ..2, ..3)
    )
  }, times = 10
)

# Unit: seconds
# expr      min       lq
# {     test$ec1 <- future_pmap(list(test$query, test$family, test$subfamily),          ~assignEC(..1, ..2, ..3)) } 129.0458 131.2588
# mean   median       uq      max neval
# 132.8799 132.8418 133.9644 138.1316    10

microbenchmark::microbenchmark(
  {
    test <- test %>% 
      mutate(
        ec2 = future_pmap(list(query, family, subfamily), ~ assignEC(..1, ..2, ..3))
      )
  }, times = 10
)

plan(sequential)

filter(
  dbcansub_cazy_substrate,
  ecs %in% "4.2.2.-" & 
    map_lgl(cazy_protein_id, ~ any(.x %in% "hello"))
)$substrates %>%
  unwrap()


# Benchmarking full function list
TMP.NODE <- 'NODE_1003_length_50341_cov_2.472282_29'
TMP.FAMILY <- 'GT5'
TMP.SUBFAMILY <- 'GT5'

TMP.F <- function(NODE, FAMILY, SUBFAMILY) {
  EC <- assignEC(NODE, FAMILY, SUBFAMILY)
  G1 <- getHitsCAZY(NODE, SUBFAMILY, 'best')
  G2 <- getHitsCAZY(NODE, SUBFAMILY, 'best_ec')
  CL <- getHitsDBCANSUB(NODE, FAMILY)
  SB <- getSubstrate(FAMILY, SUBFAMILY, CL, G2, EC)
  return(
    list(
      'NODE' = NODE,
      'FAMILY' = FAMILY,
      'SUBFAMILY' = SUBFAMILY,
      'DBCANSUB_CLUSTER' = CL,
      'CAZY' = G1,
      'CAZY_EC' = G2,
      'EC' = EC,
      'SUBSTRATE' = SB
    )
  )
}

profvis::profvis(assignEC(TMP.NODE, TMP.FAMILY, TMP.SUBFAMILY))
profvis::profvis(
  {EC <- assignEC(TMP.NODE, TMP.FAMILY, TMP.SUBFAMILY)
   G1 <- getHitsCAZY(TMP.NODE, TMP.SUBFAMILY, 'best')
   G2 <- getHitsCAZY(TMP.NODE, TMP.SUBFAMILY, 'best_ec')
   CL <- getHitsDBCANSUB(TMP.NODE, TMP.FAMILY)
   SB <- getSubstrate(TMP.FAMILY, TMP.SUBFAMILY, CL, G2, EC)}
)
microbenchmark:::autoplot.microbenchmark(mc)

