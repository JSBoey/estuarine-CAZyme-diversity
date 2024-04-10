# Correlation analyses for sediment transcriptomes

# Environment ----
library(edgeR)
library(Hmisc)
library(tidyverse)

# Import data ----
ann <- read_tsv("results/allbins_pred.annotation_table.tsv.gz")

wts <- list(
  "all" = read_tsv("results/WTS_clean_count.tsv.gz"),
  "mean" = read_tsv("results/WTS_sed_mean.tsv.gz"),
  "sum" = read_tsv("results/WTS_sed_sum.tsv.gz")
)

# Preprocessing ----
# 1. Filter to only GH and PL
cazy <- ann$node[str_detect(ann$dbcan_label, "GH|PL")]

wts_cazy <- map(wts, \(x) {
  select(x, Geneid, starts_with("Sed")) %>% 
    filter(Geneid %in% cazy)
})

# 2. Convert to matrix ----
wts_cazy <- map(wts_cazy, \(x) {
  column_to_rownames(x, "Geneid") %>% 
    as.matrix()
})

# 3. Filter data ----
# Thresholds:
# - Present in more than 1/3 of samples
wts_cazy <- map(wts_cazy, \(x) {
  nSample <- ncol(x) / 3
  x[rowSums(x > 0) >= nSample, ]
})

# Check matrix density and number of transcripts retained
walk(wts_cazy, \(x) {
  mDensity <- round(sum(x > 0) * 100 / length(x), 2)
  nTranscript <- nrow(x)
  
  cat(str_glue("Matrix density: {mDensity}"), "\n", 
      str_glue("Transcripts: {nTranscript}"), "\n")
})

# 4. Normalise data ----
# For this, I will use geTMM with TMMwsp and output as CTF
norm_wts_cazy <- map(wts_cazy, \(x) {
  # Rejoin with metadata
  metadata <- select(wts$all, !starts_with(c("Filt", "Sed"))) %>% 
    filter(Geneid %in% rownames(x))
  
  joined_df <- as.data.frame(x) %>% 
    rownames_to_column("Geneid") %>% 
    left_join(metadata, ., by = "Geneid")
  
  # Normalise
  norm_tmm(joined_df, method = "TMMwsp", output = "ctf", gene_length_correction = T)
})

# 5. Correlation matrices ----
cor_wts_cazy <- map(norm_wts_cazy, \(x) {
  count_matrix <- split_table(x)$count
  rcorr(count_matrix, type = "pearson")
})


