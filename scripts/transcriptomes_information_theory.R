# Transcriptomes through information theory

# Preamble ----
# This is a collection of functions derived from the work of Martinez & Reyes-Valdéz (2008). They describe the specificity of a gene to a specific transcriptomic profile and its potential divergence from the global transcriptome using information theory by introducing variations of the Shannon entropy.

# Notation ----
# There are a set of samples j = (1, 2, ..., N)
# There are a set of genes i = (1, 2, ..., g)

# Environment ----
library(vegan)

# Average frequency of a gene across a matrix of samples ----
aveGeneFreq <- function(count_matrix) {
  # Per sample frequency table 
  freq_table <- proportions(count_matrix, 2)
  # Number of samples
  N <- ncol(freq_table)
  # Average gene frequencies across samples
  p_i <- rowSums(freq_table) / N
  
  return(p_i)
}

# Gene specificity ----
# A measure of how specific the transcription profile of a gene is to a particular sample or set of samples. It takes a value of 0 if the gene is evenly transcribed across all samples; it has a maxima of log2(N) if the gene is only transcribed in one sample. 

geneSpecificity <- function(count_matrix, log.base = 2) {
  # Per sample frequency table
  freq_table <- proportions(count_matrix, 2)
  # Number of samples
  N <- ncol(freq_table)
  # Per gene average frequency across samples
  p_i <- aveGeneFreq(count_matrix = count_matrix)
  # Calculate frequency ratios
  r <- sweep(freq_table, 1, p_i, "/")
  # Gene specificity
  S_i <- rowSums(r * log(r, base = log.base), na.rm = TRUE) / N
  
  return(S_i)
}

# Sample specificity ----
# A measure of how specific the transcriptomic profile of the sample is compared to other samples given the specificity of genes transcribed across the samples. It takes a value of 0 if all genes expressed have the same transcription profiles everywhere else; it has maxima of log2(t) when the transcription profile of the sample is not found in any other samples.

sampleSpecificity <- function(count_matrix, log.base = 2) {
  # Per sample frequency table
  freq_table <- proportions(count_matrix, 2)
  # Gene specificity
  S_i <- geneSpecificity(count_matrix = count_matrix, log.base = log.base)
  # Average gene specificity
  aS_i <- sweep(freq_table, 1, S_i, "*")
  # Sample specificity
  d_j <- colSums(aS_i, na.rm = FALSE)
  
  return(d_j)
}

# Average transcript log frequencies of a sample ----
aveSampleFreq <- function(count_matrix, log.base = 2) {
  # Per sample frequency table
  freq_table <- proportions(count_matrix, 2)
  # Log-scaled per gene average frequency across samples
  lp_i <- log(aveGeneFreq(count_matrix = count_matrix), base = log.base)
  # Per sample average transcript log frequencies
  H_Rj <- - colSums(sweep(freq_table, 1, lp_i, "*"), na.rm = T)
  
  return(H_Rj)
}

# Kullback-Leiber divergence of a sample ---- 
KLDivSample <- function(count_matrix, log.base = 2) {
  H_Rj <- aveSampleFreq(count_matrix = count_matrix, log.base = log.base)
  H_j <- diversity(count_matrix, MARGIN = 2, index = "shannon")
  
  D_j <- H_Rj - H_j
  
  return(D_j)
}


