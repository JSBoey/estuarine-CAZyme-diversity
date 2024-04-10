# Make sure to load data from differential_expression.R first!

library(energy)

partial_correlation <- clean_data %>% 
  mutate(
    pcor_kd = list(
      pcor.test(data$tpm.wts, data$salinity, data$tpm.wgs, method = "kendall")
    )
  )

partial_correlation_results <- partial_correlation %>% 
  dplyr::select(type, node, pcor_kd) %>% 
  unnest(pcor_kd) %>% 
  filter(p.value < 0.05)

distance_correlation <- clean_data %>% 
  mutate(
    dcor = list(
      pdcor.test(data$tpm.wts, data$salinity, data$tpm.wgs, R = 199)
    )
  )

distance_correlation_results <- distance_correlation %>% 
  dplyr::select(type, node, dcor) %>% 
  mutate(
    dcor = list(tidy(dcor))
  ) %>% 
  unnest(dcor) %>% 
  filter(p.value < 0.05)

diff_test <- distance_correlation_results %>% 
  filter(!node %in% partial_correlation_results$node) %>% 
  arrange(estimate)

plots <- vector("list", length = nrow(diff_test))
for (i in 1:nrow(diff_test)) {
  noi <- diff_test$node[i]
  typ <- diff_test$type[i]
  est <- diff_test$estimate[i]
  pd <- clean_data %>% 
    filter(node == noi & type == typ) %>% 
    unnest(data) %>% 
    pivot_longer(c(tpm.wgs, salinity), names_to = "f", values_to = "v")
  plots[[i]] <- ggplot(pd, aes(x = v, y = tpm.wts)) +
    geom_point() +
    labs(title = paste(noi, typ, est)) +
    facet_wrap(~ f, scales = "free")
}
pdf("Distance_correlation.pdf", paper = "a4r")
plots
dev.off()



noi <- "NODE_10_length_501835_cov_18.196130_36"
clean_data %>% 
  filter(node == noi & type == "water") %>% 
  unnest(data) %>% 
  pivot_longer(c(tpm.wgs, salinity), names_to = "f", values_to = "v") %>% 
  ggplot(aes(x = v, y = tpm.wts)) +
  geom_point() +
  labs(title = noi) +
  facet_wrap(~ f, scales = "free")
