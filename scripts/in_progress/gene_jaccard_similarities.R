# Gene/transcript-wise Jaccard similarity

# Which genes and transcripts overlap significantly in terms of their 
# distributions?

# Environment ----
library(jaccard)
library(vegan)
library(tidyverse)

# Import data ----
ann <- read_tsv("results/allbins_pred.annotation_table.tsv.gz")
data <- map(list.files("results/", ".*_clean_count.tsv.gz", full.names = T), 
            \(x) {
              read_tsv(x)
            }) %>% 
  set_names(c("wgs", "wts"))

# Preprocessing ----
# 1. Split data into water and sediment
data <- map(data, \(x) {
  water <- select(x, Geneid, contains(c("Filt")))
  sediment <- select(x, Geneid, contains(c("Sed")))
  list("water" = water, "sediment" = sediment)
})

data <- list_flatten(data)

# 2. Convert matrices into binary matrix
data <- map(data, \(x) {
  X <- as.matrix(column_to_rownames(x, "Geneid"))
  X[X > 0] <- 1
  
  return(X)
})

# 3. Retain genes/transcript that are present in more than 20% of samples
data <- map(data, \(x) {
  n <- ncol(x) / 5
  x[rowSums(x) > n, ]
})

# (Optional) 4. Retain only CAZymes (GH and PL only)
cazy <- map(data, \(x) {
  cazymes <- ann$node[!is.na(ann$dbcan_label)]
  x[rownames(x) %in% cazymes, ]
})

# Calculate Jaccard similarity ----
calc_Jaccard <- function(x) {
  1 - vegdist(x, method = "jaccard")
}

# Calculating for CAZymes in the transcriptome
jac <- map(cazy, calc_Jaccard)

test_data <- cazy$wts_sediment[rowSums(cazy$wts_sediment) < ncol(cazy$wts_sediment), ]

jac_test <- jaccard.test.pairwise(test_data, method = "asymptotic")
# Check intersections 




