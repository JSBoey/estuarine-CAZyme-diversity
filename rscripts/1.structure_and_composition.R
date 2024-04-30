# Ordinations

# 0. Set up

# Normalisation is important to ensure distances are properly calculated
g <- ANNOTATION %$%
  (end - start + 1)
names(g) <- ANNOTATION$node

# 1. Structure and composition of GH and PL encoding genes and transcripts
#    in the environment.

distance_method <- c("structure" = "bray", 
                     "composition" = "jaccard")

present_threshold <- 2 # 0 if less than present_threshold

distance_matrix <- lapply(
  distance_method, \(method) {
    binary <- ifelse(method == "jaccard", TRUE, FALSE)
    lapply(COUNTS, \(x) {
      if (binary)
        x[x < present_threshold] <- 0
      else
        x <- calc_TPM(x, g)
      
      x <- split_matrix(x, rows = get_GHPL(), remove_zero_rowsum = TRUE)
      vegdist(t(x), method = method, binary = binary)
    })
  }
)

nmds <- lapply_at_depth(
  distance_matrix, 2, 
  \(x) metaMDS(x, try = 999, trymax = 1999, autotransform = FALSE)
)

nmds_scores <- lapply_at_depth(nmds, 2, \(x) scores(x, tidy = TRUE)) %>% 
  unlist(recursive = FALSE) %>% 
  rbindlist(idcol = "data_type") %>% 
  merge(., METADATA[, .(sample, type, salinity)], 
        by.x = "label", by.y = "sample")

nmds_scores[
  , c("dist_type", "source_type") := tstrsplit(data_type, ".", fixed = TRUE)
]

nmds_plot <- ggplot(
  nmds_scores, 
  aes(x = NMDS1, y = NMDS2, colour = salinity, shape = type)
) +
  geom_point() +
  scale_colour_viridis_c() +
  scale_shape_manual(
    values = habitat_shape,
    labels = stringr::str_to_title
  ) +
  labs(
    colour = "Salinity (ppt)",
    shape = "Habitat substrate"
  ) +
  facet_wrap(
    dist_type ~ source_type,
    scales = "free"
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    legend.position = "right",
    aspect.ratio = 1,
    axis.title = element_blank()
  )

nmds_plot

walk(
  c("svg", "png", "tiff"),
  \(x) {
    filename <- glue("figures/ordination.{x}")
    ggsave(filename, nmds_plot, width = 2.5, height = 2, dpi = 600, scale = 2.5)
  }
)

sink("statistics/ordinations.txt")
cat("#### Ordination stress values ####\n")
lapply_at_depth(nmds, 2, \(i) i$stress)
cat("\n#### PERMANOVA ####\n")
lapply_at_depth(distance_matrix, 2, \(i) {
  sample_data <- METADATA[sample %in% labels(i)]
  adonis2(i ~ type + salinity, sample_data, permutations = 999, by = "margin")
})
cat("\n#### Beta-dispersion ####\n")
lapply_at_depth(distance_matrix, 2, \(i) {
  sample_data <- METADATA[sample %in% labels(i)]
  permutest(betadisper(i, sample_data$type, bias.adjust = TRUE))
})
cat("Results were computed from 1.structure_and_composition.R")
cat("\n#### EOF ####\n")
sink()

# 2. Association between genes and transcripts across the environment

## Qualitative

procrustes_test <- lapply(
  nmds, \(ord) {
    overlap_samples <- Reduce(intersect, lapply(ord, \(i) rownames(i$points)))
    scr <- lapply(ord, \(i) split_matrix(scores(i), rows = overlap_samples))
    protest(scr$WGS, scr$WTS)
  }
)

procrustes_scores <- lapply(
  procrustes_test, \(x) {
    axes <- get_procrustes_axes(x)
    data <- as.data.table(cbind(x$X, x$Yrot), keep.rownames = "sample")
    names(data)[-1] <- c("x_gene", "y_gene", "x_transcript", "y_transcript")
    data <- merge(data, METADATA[, .(sample, type, salinity)], 
                  by = "sample", all.x = TRUE)
    list("axes" = axes, "coordinates" = data)
  }
)

procrustes_plot <- map2(
  procrustes_scores, names(procrustes_scores),
  \(x, nm) {
    
    line_colour <- "grey70"
    segment_colour <- "navy"
    
    ggplot() +
      geom_hline(yintercept = 0, linetype = 2, colour = line_colour) +
      geom_vline(xintercept = 0, linetype = 2, colour = line_colour) +
      geom_abline(
        data = x$axes,
        mapping = aes(intercept = intercept, slope = slope),
        colour = line_colour
      ) +
      geom_segment(
        data = x$coordinates,
        mapping = aes(x = x_gene,
                      y = y_gene,
                      xend = x_transcript,
                      yend = y_transcript),
        colour = segment_colour,
        arrow = arrow(
          angle = 30,
          length = unit(3, "points")
        )
      ) +
      geom_point(
        data = x$coordinates,
        mapping = aes(x = x_gene,
                      y = y_gene,
                      colour = salinity,
                      shape = type)
      ) +
      labs(
        title = stringr::str_to_title(nm),
        colour = "Salinity (ppt)",
        shape = "Habitat substrate"
      ) +
      scale_colour_viridis_c() +
      scale_shape_manual(
        values = habitat_shape,
        labels = stringr::str_to_title
      ) +
      theme_bw() +
      theme(
        aspect.ratio = 1,
        panel.grid = element_blank(),
        axis.title = element_blank()
      )
    
  }
)

sink("statistics/procrustes_test.txt")
cat("#### Symmetric Procrustes test ####\n")
procrustes_test
cat("Results were computed from 1.structure_and_composition.R")
cat("\n#### EOF ####\n")
sink()

# Quantitative

procrustes_resids <- lapply(procrustes_test, \(x) {
  residual <- residuals(x)
  data.table(
    "sample" = names(residual),
    "procrustes_residual" = residual
  ) %>% 
    merge(., METADATA[, .(sample, type, salinity)],
          by = "sample", all.x = TRUE)
})

residual_model <- lapply(procrustes_resids, \(dt) {
  dt <- dt[
    , type := factor(type, 
                     levels = c("water", "sediment"), 
                     ordered = FALSE)
  ]
  
  LM <- lm(procrustes_residual ~ type * salinity, dt)
  ANOVA <- Anova(LM, type = "II")
  HOMOGENEITY <- leveneTest(procrustes_residual ~ type, dt)
  
  list(
    "regression" = LM,
    "anova" = ANOVA,
    "homogeneity" = HOMOGENEITY
  )
})

calc_outlier_threshold <- function(x) quantile(x, .75) + 1.5 * IQR(x)

residual_model_sans_outlier <- lapply(procrustes_resids, \(dt) {
  outlier <- calc_outlier_threshold(dt$procrustes_residual)
  dt <- dt[
    procrustes_residual <= outlier, 
    type := factor(type, 
                   levels = c("water", "sediment"), 
                   ordered = FALSE)
  ]
  
  LM <- lm(procrustes_residual ~ type + salinity, dt)
  ANOVA <- Anova(LM, type = "II")
  HOMOGENEITY <- leveneTest(procrustes_residual ~ type, dt)
  
  list(
    "regression" = LM,
    "anova" = ANOVA,
    "homogeneity" = HOMOGENEITY
  )
})

# Models are not different

plot_procrustes_model <- function(data, model, title) {
  coefs <- round(coef(model$regression), 3)
  caption <- glue(
    "Slopes: salinity = {coefs[3]}; habitat_water = {coefs[2]}"
  )
  ggplot(data, aes(x = salinity,
                   y = procrustes_residual,
                   colour = type)) +
    geom_point() +
    geom_smooth(method = "lm", se = FALSE) +
    scale_colour_manual(values = habitat_colour) +
    labs(
      x = "Salinity (ppt)",
      y = "Procrustes residuals",
      colour = "Habitat substrate",
      title = title,
      caption = caption
    ) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      aspect.ratio = 1
    )
}

procrustes_model_plot <- mapply(plot_procrustes_model, 
                                procrustes_resids, 
                                residual_model, 
                                names(residual_model),
                                SIMPLIFY = FALSE)

# 3. Compositional differences between genes and transcripts.

binary_ghpl <- lapply(COUNTS, \(m) {
  m <- split_matrix(m, rows = get_GHPL())
  ifelse(m > present_threshold, 1, 0)
})

binary_ghpl_order <- function() {
  h <- vegdist(binary_ghpl$WGS, method = "jaccard") %>%
    hclust(method = "ward.D2")
  h$labels[h$order]
}

ghpl_composition_heatmap <- lapply(binary_ghpl, \(m) {
  as.data.table(m, keep.rownames = "node") %>% 
    melt(id.vars = "node", 
         measure.vars = patterns("Sed|Filt"),
         variable.name = "sample", 
         value.name = "value")
}) %>% 
  rbindlist(idcol = "src") %>% 
  .[sample %in% colnames(COUNTS$WTS)] %>% 
  ggplot() +
    geom_tile(
      aes(x = sample,
          y = factor(node, levels = binary_ghpl_order()),
          fill = as.factor(value))
    ) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_fill_viridis_d(option = "rocket", 
                         labels = c("0" = "Absent", "1" = "Present")) +
    facet_grid(cols = vars(src), 
               space = "fixed", 
               scales = "fixed") +
    labs(fill = "Gene/Transcripts detected") +
    theme_bw() +
    theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      legend.position = "bottom"
    )

replicate_metadata <- function() {
  dt <- copy(METADATA[sample %in% colnames(COUNTS$WTS), 
                      .(sample, type, salinity)])
  rbindlist(list("WGS" = dt, "WTS" = dt), idcol = "src")
}

type_meta_heatmap <- ggplot(
  replicate_metadata(),
  aes(x = sample, y = "type", fill = type)
) +
  geom_tile() +
  scale_x_discrete(expand = c(0, 0)) +
  scale_fill_manual(values = habitat_colour) +
  facet_grid(cols = vars(src)) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "bottom",
    strip.background = element_blank(),
    strip.text = element_blank()
  )

salinity_meta_heatmap <- ggplot(
  replicate_metadata(),
  aes(x = sample, y = "salinity", fill = salinity)
) +
  geom_tile() +
  scale_x_discrete(expand = c(0, 0)) +
  scale_fill_viridis_c() +
  facet_grid(cols = vars(src)) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "bottom",
    strip.background = element_blank(),
    strip.text = element_blank()
  )

ghpl_composition_heatmap / type_meta_heatmap / salinity_meta_heatmap +
  plot_layout(nrow = 3, heights = c(20, 1, 1), guides = 'collect') &
  theme(legend.position = "bottom")

# Addressing comments for the section

# Hence, although freshwater taxa appear able to use a diverse range of
# carbohydrate substrates (on average ???? +/- SD genes per MAG, Fig ????), 
# results indicate a preference for particular substrates over others, and 
# therefore uneven transcription of associated genes.
richness <- sapply(binary_ghpl, \(i) colSums(i[, 1:30])) %>% 
  as.data.table(keep.rownames = "sample") %>% 
  merge(., METADATA[, .(sample, type, salinity)], all.x = TRUE, by = "sample") %>% 
  .[, ratio := WTS/WGS]
lm_richness <- glm(ratio ~ salinity + type, family = "Gamma", richness)
summary(lm_richness)
car::Anova(lm_richness, type = 2)

melt(richness, 
     id.vars = "sample", 
     measure.vars = c("WGS", "WTS"), 
     variable.name = "")

caz_per_mag <- lapply(binary_ghpl, \(i) {
  as.data.table(i, keep.rownames = "node") %>% 
    merge(., ANNOTATION[, .(node, bin)], all.x = TRUE, by = "node") %>% 
    .[, lapply(.SD, sum), by = "bin", .SDcols = patterns("Sed|Filt")] %>% 
    as.matrix(rownames = "bin")
})

lapply(caz_per_mag, \(i) {
  i <- ifelse(i > 0, i, NA)
  data.table(
    "average" = colMeans(i, na.rm = T)[1:30],
    "stdev" = colVars(i, na.rm = T, std = T)[1:30],
    "n" = colSums(i > 0, na.rm = T)[1:30]
  ) %>% 
    .[, covar := stdev/average]
})
