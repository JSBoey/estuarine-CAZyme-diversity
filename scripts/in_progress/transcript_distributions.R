# Determine threshold for filtering out genes

# Not all transcripts are useful for correlation type analyses. In this dataset,
# sparsity is a real issue in generating useful gene-wise correlations. This 
# explores gene distributions and aims to define a threshold for filtering
# transcripts for downstream analyses.

# Load libraries ----
library(vegan)
library(tidyverse)

# Import data ----
ann <- read_tsv("results/allbins_pred.annotation_table.tsv.gz")
wts <- read_tsv("results/WTS_clean_count.tsv.gz")
wts_sum <-read_tsv("results/WGS_sed_sum.tsv.gz")

# Split and filter data ----
# 1. Split data into water and sediment samples
# 2. Coerce to count matrix
# 3. Remove transcripts with no reads

data_matrix <- map(c("water" = "Filt", "sediment" = "Sed"), \(x) {
  
  # Split
  a <- select(wts, Geneid, contains(x))
  
  # Matrix
  b <- as.matrix(column_to_rownames(a, "Geneid"))
  
  # Remove zero sums
  b[rowSums(b) > 0, ]
  
})

# Append summed counts to data_matrix ----
data_matrix$sum_sediment <- wts_sum %>% 
  select(Geneid, contains("Sed")) %>% 
  column_to_rownames("Geneid") %>% 
  as.matrix() %>% 
  .[rowSums(.) > 0, ]

# Transcripts retained based on sample presence ----
# How many transcripts are retained if they are filtered based on how many 
# samples they must be present in?

sample_fractions <- seq(0.1, 0.9, by = 0.1)

threshold_data <- tibble(
  "sample_type" = rep(names(data_matrix), each = length(sample_fractions)),
  "sample_fraction_removed" = rep(sample_fractions, times = length(data_matrix)),
  "transcripts_retained" = NA,
  "transcripts_retained_percentage" = NA,
  "average_transcript_diversity" = NA,
  "sd_transcript_diversity" = NA,
  "log10_correlation_size" = NA
)

for (i in seq_len(nrow(threshold_data))) {
  
  # Extract relevant matrix
  sample_type <- threshold_data$sample_type[i]
  sample_fraction_removed <- threshold_data$sample_fraction_removed[i]
  M <- pluck(data_matrix, sample_type)
  
  # Filter based on presence threshold
  threshold <- ncol(M) * sample_fraction_removed
  fM <- M[rowSums(M > 0) > threshold, ]
  
  # Results
  threshold_data$transcripts_retained[i] <- nrow(fM)
  threshold_data$transcripts_retained_percentage[i] <- nrow(fM) * 100 / nrow(M)
  H <- diversity(t(fM))
  threshold_data$average_transcript_diversity[i] <- mean(H)
  threshold_data$sd_transcript_diversity[i] <- sd(H)
  threshold_data$log10_correlation_size[i] <- log10((nrow(fM)^2 - nrow(fM)) / 2)
  
  # Remove things
  rm(i, sample_type, sample_fraction_removed, M, threshold, fM, H)
  
}

# Plotting
ggplot(data = threshold_data, aes(x = sample_fraction_removed, colour = sample_type)) +
  geom_line(aes(y = transcripts_retained))

ggplot(data = threshold_data, aes(x = sample_fraction_removed, colour = sample_type)) +
  geom_line(aes(y = transcripts_retained_percentage))

ggplot(data = threshold_data, aes(x = sample_fraction_removed, y = average_transcript_diversity, colour = sample_type)) +
  geom_point() +
  geom_errorbar(aes(ymin = average_transcript_diversity - sd_transcript_diversity,
                    ymax = average_transcript_diversity + sd_transcript_diversity), width = 0.005) +
  geom_line()

ggplot(data = threshold_data, aes(x = sample_fraction_removed, y = log10_correlation_size, colour = sample_type)) +
  geom_line()

# "Lonely" transcripts ----
# How many transcripts are only present in one sample?
lonely_transcripts <- map_dbl(data_matrix, \(x) {
  nrow(x[rowSums(x > 0) == 1, ])
})

# How large must the correlation matrix be if lonely transcripts were filtered
# out?
correlation_matrix_size <- map2_dbl(data_matrix, lonely_transcripts, \(x, y) {
  # Transcripts retained 
  a <- nrow(x) - y
  
  # Correlation matrix size in log10
  log10(((a^2) - a) / 2)
})


