# Differential expression

# Preamble ----
# The analysis here will follow suggestions and methods derived from Zhang et al. (2021) 10.1093/bioinformatics/btab327 
# Firstly, it is important to note that we are lucky to have paired metagenome/transcriptome data. This makes the following model accessible to us:
# RNA ~ DNA + var1 + var2 + ...
# In microbial ecology, log-transformation is common and easy to interpret. 
# The other concern is the non-linearity of transcripts that might be highly expressed in intermediate salinities as opposed to fresh/saline extremes.

# Environment ----
library(matrixStats)
library(gamlss)
library(brms)
library(tidyverse)
library(broom)

# Import data ----
load("results/transformed_counts.RData")
metadata <- read_tsv("data/sample_metadata.txt")
annotation <- read_tsv("results/allbins_pred.annotation_table.tsv.gz")

# Equalise samples ----
transform_count <- transform_count[grepl("TPM", names(transform_count))]
transform_count$wgs_TPM <- transform_count$wgs_TPM[,
  colnames(transform_count$wgs_TPM) %in% colnames(transform_count$wts_TPM)
]
names(transform_count) <- str_remove(names(transform_count), "_TPM")

# Split data ----
split_data <- map(transform_count, function(MAT) {
  sample_type <- split(metadata$sample, f = metadata$type)
  map(sample_type, function(nm) MAT[, colnames(MAT) %in% nm])
}) %>% 
  transpose()

# Filter data ----
cazy_nodes <- annotation %>% 
  filter(str_detect(dbcan_label, "GH|PL")) %>% 
  pull(node)

clean_matrices <- map_depth(split_data, 2, function(MAT) {
  # Retain rows that are non-zero sum and are CAZymes
  MAT[rowSums(MAT) > 0 & rownames(MAT) %in% cazy_nodes, ]
}) %>% 
  map(function(LIST) {
    # Retain rows in WGS that are present in WTS
    LIST$wgs <- LIST$wgs[rownames(LIST$wgs) %in% rownames(LIST$wts),]
    LIST$wts <- LIST$wts[rownames(LIST$wts) %in% rownames(LIST$wgs),]
    
    LIST
  })

# Transform, clean, and nest data ----
clean_data <- map_depth(clean_matrices, 2, function(MAT) {
  as.data.frame(MAT) %>% 
    as_tibble(rownames = "node") %>% 
    pivot_longer(-node, names_to = "sample", values_to = "tpm")
}) %>% 
  map(reduce, left_join, by = c("node", "sample"), suffix = c(".wgs", ".wts")) %>% 
  map(function(tb) {
    left_join(tb, select(metadata, sample, salinity)) %>% 
      nest_by(node)
  }) %>% 
  bind_rows(.id = "type")

# Model: Logistic regression ----
# Here's a simpler question: Is expression related to salinity? 
# This can be modeled using a logistic regression.
binary_data <- clean_data %>% 
  mutate(
    data = list(
      mutate(data, tpm.wts = ifelse(tpm.wts > 0, 1, 0)) %>% 
        rename("wts" = "tpm.wts")
    )
  )

logis_model <- binary_data %>%
  filter(sum(data$wts) > 1) %>% # Remove single response
  transmute(
    type, node,
    model = list(
      glm(wts ~ tpm.wgs + salinity, data, family = binomial(link = 'logit'))
    ),
    model_statistics = list(glance(model)),
    model_results = list(tidy(model)),
    model_terms = list(tidy(anova(model, test = "Chisq")))
  )

logis_model_statistics <- logis_model %>% 
  select(type, node, model_statistics) %>% 
  unnest(cols = model_statistics)

logis_model_results <- logis_model %>% 
  select(type, node, model_results) %>% 
  unnest(cols = model_results)

logis_model_terms <- logis_model %>% 
  select(type, node, model_terms) %>% 
  unnest(cols = model_terms)

logis_model_terms_sig_salinity <- logis_model_terms %>% 
  filter(term == 'salinity' & p.value < 0.05)

logis_model_results_sig_salinity <- logis_model_results %>% 
  filter(term == "salinity" & node %in% logis_model_terms_sig_salinity$node)

binary_data %>% 
  filter(node == "NODE_953_length_40353_cov_8.418771_8" & type == "water") %>% 
  unnest(data) %>% 
  pivot_longer(c(tpm.wgs, salinity), names_to = "f", values_to = "v") %>% 
  ggplot(aes(x = v, y = wts)) +
  geom_point() +
  stat_smooth(method = "glm", method.args = list(family = binomial("logit"))) +
  labs(title = "NODE_953_length_40353_cov_8.418771_8") +
  facet_wrap(~ f, scales = "free")

rank_clean <- clean_data %>% 
  mutate(
    data = list(mutate(data, r = rank(tpm.wts)))
  ) %>% 
  unnest(data) %>% 
  group_by(node, type) %>% 
  filter(r == max(r) & between(salinity, 1, 30) & type == "water")

clean_data %>% 
  filter(node == "NODE_100_length_138186_cov_4.921704_45" & type == "sediment") %>% 
  unnest(data) %>% 
  pivot_longer(c(tpm.wgs, salinity), names_to = "f", values_to = "v") %>% 
  ggplot(aes(x = v, y = rank(tpm.wts))) +
  labs(title = "NODE_100_length_138186_cov_4.921704_45") +
  facet_wrap(~ f, scales = "free")

brackish_example <- clean_data %>% 
  filter(node == "NODE_105_length_130382_cov_18.754032_71" & type == "water") %>% 
  unnest(data)

brackish_model <- glm(
  rank(tpm.wts) ~ tpm.wgs + poly(salinity, 3), 
  family = Gamma(link = "log"), data = brackish_example
)

anova(brackish_model, test = "F")

tidy(brackish_model)

# Model 2: Beta regression ----
# Omic data is compositional, thus it is technically bounded. The fact that we use TPM here only changes the closure from unity to 1 million. If we convert TPM values back to a closure of unity, that means we have sum-scaled RPK. I am less inclined to use raw count values given the biases of read recruitment due to gene length. Therefore, the nucleic acid measures will be in sum-scaled RPK. 
# The bounds for a beta-regression is for values between, but not including, 0 and 1. However, due to the nature of sequencing and, well, biology/ecology, there will be 0s. Hence, it might be good to opt for a beta-distributed regression.

## Formula
f <- tpm.wts ~ tpm.wgs + salinity

## Frequentist Beta regression
freq_ziβ <- clean_data[1:9, ] %>% 
  mutate(
    model = list(glm(f, family = BEZI(), data = data)),
    model_statistics = list(broom::glance(model)),
    model_results = list(broom::tidy(model))
  )

## Bayesian Beta regression
bayes_ziβ <- clean_data[1:10, ] %>% 
  mutate(
    model = list(
      brm(
        bf(tpm.wts ~ salinity, phi ~ salinity, zi ~ tpm.wgs),
        data = data, 
        family = zero_inflated_beta(),
        chains = 4, iter = 2000, warmup = 1000,
        cores = 4, seed = 1234,
        file = "model_beta_zi"
      )
    )
  )


