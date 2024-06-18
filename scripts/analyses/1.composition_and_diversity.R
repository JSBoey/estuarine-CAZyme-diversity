# Ordinations and diversity

# 0. Set up

# Normalisation is important to ensure distances are properly calculated
g <- ANNOTATION %$%
  (end - start + 1)
names(g) <- ANNOTATION$node

# 1. Structure and composition of GH and PL encoding genes and transcripts
#    in the environment.

distance_method <- c("structure" = "bray", 
                     "composition" = "jaccard")

present_threshold <- 1 # 0 if less than present_threshold

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

# walk(
#   c("svg", "png", "tiff"),
#   \(x) {
#     filename <- glue("figures/ordination.{x}")
#     ggsave(filename, nmds_plot, width = 2.5, height = 2, dpi = 600, scale = 2.5)
#   }
# )

# sink("statistics/ordinations.txt")
# cat("#### Ordination stress values ####\n")
# lapply_at_depth(nmds, 2, \(i) i$stress)
# cat("\n#### PERMANOVA ####\n")
# lapply_at_depth(distance_matrix, 2, \(i) {
#   sample_data <- METADATA[sample %in% labels(i)]
#   adonis2(i ~ type + salinity, sample_data, permutations = 999, by = "margin")
# })
# cat("\n#### Beta-dispersion ####\n")
# lapply_at_depth(distance_matrix, 2, \(i) {
#   sample_data <- METADATA[sample %in% labels(i)]
#   permutest(betadisper(i, sample_data$type, bias.adjust = TRUE))
# })
# cat("Results were computed from 1.structure_and_composition.R")
# cat("\n#### EOF ####\n")
# sink()

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
    
    nm <- ifelse(nm == "composition", "Presence/absence", "Abundance-weighted")
    
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
                      shape = type),
        size = 2
      ) +
      labs(
        title = stringr::str_to_sentence(nm),
        colour = "Salinity (ppt)",
        shape = "Habitat substrate"
      ) +
      scale_colour_viridis_c() +
      scale_shape_manual(
        values = habitat_shape,
        labels = stringr::str_to_sentence
      ) +
      theme_bw() +
      theme(
        aspect.ratio = 1,
        panel.grid = element_blank(),
        axis.title = element_blank()
      )
    
  }
)

procrustes_plot_together <- wrap_plots(procrustes_plot, ncol = 2) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

# sink("statistics/procrustes_test.txt")
# cat("#### Symmetric Procrustes test ####\n\n")
# procrustes_test
# cat("Results were computed from 1.structure_and_composition.R\n")
# cat("\n#### EOF ####\n")
# sink()
# 
# walk(
#   c("svg", "png", "tiff"),
#   \(x) {
#     filename <- glue("figures/procrustes.{x}")
#     ggsave(filename, procrustes_plot_together, 
#            width = 2.5, height = 2, dpi = 600, scale = 2.5)
#   }
# )


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
  DT <- dt[
    , type := factor(type, 
                     levels = c("water", "sediment"), 
                     ordered = FALSE)
  ]
  
  LM <- lm(procrustes_residual ~ type * salinity, DT)
  ANOVA <- Anova(LM, type = "II")
  HOMOGENEITY <- leveneTest(procrustes_residual ~ type, DT)
  
  list(
    "regression" = LM,
    "regression_summary" = summary(LM),
    "anova" = ANOVA,
    "homogeneity" = HOMOGENEITY
  )
})

calc_outlier_threshold <- function(x) quantile(x, .75) + (1.5 * IQR(x))

residual_model_sans_outlier <- lapply(procrustes_resids, \(dt) {
  outlier <- calc_outlier_threshold(dt$procrustes_residual)
  DT <- dt[
    procrustes_residual <= outlier
  ][
    , type := factor(type, 
                     levels = c("water", "sediment"), 
                     ordered = FALSE)
  ]
  
  LM <- lm(procrustes_residual ~ type * salinity, DT)
  ANOVA <- Anova(LM, type = "II")
  HOMOGENEITY <- leveneTest(procrustes_residual ~ type, DT)
  
  list(
    "regression" = LM,
    "regression_summary" = summary(LM),
    "anova" = ANOVA,
    "homogeneity" = HOMOGENEITY
  )
})

# Models are not that different, in that only salinity was significant
# Outlier Sed.S2_2 was found in the composition model

residual_model_update <- lapply(procrustes_resids, \(dt) {
  lm(procrustes_residual ~ salinity, dt)
})

residual_model_sans_outlier_update <- lapply(procrustes_resids, \(dt) {
  outlier <- calc_outlier_threshold(dt$procrustes_residual)
  dt <- dt[procrustes_residual <= outlier]
  lm(procrustes_residual ~ salinity, dt)
})

plot_residual_model_structure <- ggplot(
  procrustes_resids$structure,
  aes(x = salinity, y = procrustes_residual, colour = type)
) +
  geom_point() +
  geom_abline(
    intercept = coef(residual_model_update$structure)[[1]],
    slope = coef(residual_model_update$structure)[[2]]
  ) +
  scale_colour_manual(values = habitat_colour, labels = stringr::str_to_title) +
  labs(x = "Salinity (ppt)", 
       y = "Procrustes residuals", 
       colour = "Habitat substrate",
       title = "Abundance-weighted") +
  theme_bw() + 
  theme(
    aspect.ratio = 1,
    panel.grid = element_blank()
  )

outlier_coord <- lapply(procrustes_resids, \(dt) {
  outlier <- calc_outlier_threshold(dt$procrustes_residual)
  dt[procrustes_residual >= outlier]
})

plot_residual_model_composition <- ggplot(
  procrustes_resids$composition,
  aes(x = salinity, y = procrustes_residual, colour = type)
) +
  geom_point() +
  geom_abline(
    intercept = coef(residual_model_update$composition)[[1]],
    slope = coef(residual_model_update$composition)[[2]]
  ) +
  geom_abline(
    intercept = coef(residual_model_sans_outlier_update$composition)[[1]],
    slope = coef(residual_model_sans_outlier_update$composition)[[2]],
    linetype = 2
  ) +
  scale_colour_manual(values = habitat_colour, labels = stringr::str_to_title) +
  labs(x = "Salinity (ppt)", 
       y = "Procrustes residuals", 
       colour = "Habitat substrate",
       title = "Presence/absence") +
  annotate(
    geom = "text",
    x = outlier_coord$composition$salinity * 15,
    y = outlier_coord$composition$procrustes_residual,
    label = outlier_coord$composition$sample, 
    size = 3
  ) +
  theme_bw() + 
  theme(
    aspect.ratio = 1,
    panel.grid = element_blank()
  )

plot_procrustes_residuals <- 
  (plot_residual_model_structure | plot_residual_model_composition) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

# walk(
#   c("svg", "png", "tiff"),
#   \(x) {
#     filename <- glue("figures/procrustes_residual_model.{x}")
#     ggsave(filename, plot_procrustes_residuals, 
#            width = 2.5, height = 2, dpi = 600, scale = 2.5)
#   }
# )

# 3. Uneven transcription relative to genomic capabilities(?)

alpha_diversity <- lapply(COUNTS, \(m) {
  m <- split_matrix(m, rows = get_GHPL())
  dt <- data.table(
    "sample" = colnames(m),
    "q" = specnumber(m, MARGIN = 2),
    "N2" = diversity(m, index = "invsimpson", MARGIN = 2)
  )
}) %>% 
  rbindlist(idcol = "srctype") %>% 
  dcast(sample ~ srctype, value.var = c("q", "N2")) %>% 
  merge(., METADATA[, .(sample, type, salinity)], all.x = T, by = "sample") %>% 
  .[, type := factor(type, levels = c("water", "sediment"))]

plot_alpha_diversity <- lapply(
  c("q", "N2"),
  \(i) {
    wgs <- paste0(i, "_WGS")
    wts <- paste0(i, "_WTS")
    
    trans <- ifelse(i %in% c("q", "N2"), "log10", "identity")
    label_y <- dplyr::case_when(
      i == "q" ~ "Richness",
      i == "N2" ~ "Inverse Simpson's",
      i == "E" ~ "Evenness"
    )
    
    gg <- ggplot(alpha_diversity, aes(x = glue("{sample} ({salinity})"))) +
      geom_segment(aes(y = get(wgs), 
                       yend = get(wts))) +
      geom_point(aes(y = get(wgs)),
                 colour = "cornflowerblue") +
      geom_point(aes(y = get(wts)),
                 colour = "firebrick") +
      scale_y_continuous(transform = trans) +
      scale_colour_viridis_c() +
      labs(x = "Salinity (ppt)",
           y = label_y) +
      facet_grid(cols = vars(type),
                 scales = "free_x",
                 space = "free_x",
                 labeller = labeller(type = stringr::str_to_title)) +
      theme_bw() +
      theme(panel.grid = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))
    
    if (i == "N2") {
      gg
    } else {
      gg +
        theme(axis.text.x = element_blank(),
              axis.title.x = element_blank(),
              axis.ticks.x = element_blank())
    }
  }
)

wrap_plots_alpha_diversity <- 
  wrap_plots(plot_alpha_diversity, ncol = 1, nrow = 2, guides = 'collect')

# walk(
#   c("svg", "png", "tiff"),
#   \(x) {
#     filename <- glue("figures/alpha_diversity.{x}")
#     ggsave(filename, wrap_plots_alpha_diversity,
#            width = 2.5, height = 2.5, dpi = 600, scale = 2.5)
#   }
# )

# Pooled alpha diversity

pooled_alpha_diversity <- lapply(COUNTS, \(m) {
  m <- split_matrix(m, rows = get_GHPL(), remove_zero_rowsum = TRUE)
  g <- METADATA$type[1:ncol(m)]
  
  data.frame(
    "q" = specnumber(t(m), groups = g, MARGIN = 1),
    "N2" = diversity(m, index = "invsimpson", groups = g, MARGIN = 2)
  ) %>% 
    as.data.table(keep.rownames = "habitat")
  
}) %>% 
  rbindlist(idcol = 'srctype')

# Rank abundance

rank_abundance <- lapply(COUNTS, \(m) {
  m <- split_matrix(m, rows = get_GHPL(), remove_zero_rowsum = TRUE)
  m <- proportions(m, 2)
  
  dt <- as.data.table(m, keep.rownames = "node") %>% 
    melt(id.vars = "node",
         value.name = "p",
         variable.name = "sample")
  
  dt <- dt[
    p > 0, .SD[order(p, decreasing = TRUE)], 
    by = "sample"
  ][
    , i := 1:.N,
    by = "sample"
  ]
  
  merge(dt, METADATA[, .(sample, type, salinity)],
        all.x = TRUE, by = "sample")
}) %>% 
  rbindlist(idcol = "srctype") %>% 
  .[, srctype := as.factor(srctype)]

split_rank_abundance <- split(rank_abundance, rank_abundance$sample)

pmax <- max(rank_abundance$p)
pmin <- min(rank_abundance$p)
imax <- max(rank_abundance$i)

rac_list <- lapply(split_rank_abundance, \(dt) {
  nm <- glue("{unique(dt$sample)} ({unique(dt$salinity)} ppt)")
  # ylim <- round(log10(c(pmin, pmax + 0.1)), 1)
  # xlim <- round(c(0, imax+500)/1000, 0)
  
  breaks_x <- function(x) c(1, median(x), max(x))
  breaks_y <- function(y) c(min(y), max(y))
  
  ggplot(dt, 
         aes(x = i,
             y = p,
             colour = srctype)) +
    geom_point() +
    facet_wrap(~ srctype,
               scales = "free",
               ncol = 2,
               nrow = 1) +
    scale_y_continuous(breaks = breaks_y, 
                       transform = "log10",
                       labels = label_number(accuracy = 0.01)) +
    scale_x_continuous(breaks = breaks_x,
                       labels = label_number(accuracy = 1,
                                             big.mark = "")) +
    scale_colour_manual(values = c("WGS" = "cornflowerblue",
                                   "WTS" = "firebrick")) +
    labs(title = nm,
         x = "Rank",
         y = "Relative abundance") +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.title = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5),
          # axis.ticks = element_blank(),
          strip.background = element_blank(),
          strip.text = element_blank()
    )
})

# rank_abundance_curve <- 
#   wrap_plots(rac_list[colnames(COUNTS$WTS)], guides = "collect")

rank_abundance_curve_water <- 
  wrap_plots(rac_list[grep("Filt", colnames(COUNTS$WTS), value = TRUE)], 
             guides = "collect")

rank_abundance_curve_sediment <- 
  wrap_plots(rac_list[grep("Sed", colnames(COUNTS$WTS), value = TRUE)], 
             guides = "collect")

# walk(
#   c("svg", "png", "tiff"),
#   \(x) {
#     filename <- glue("figures/rank_abundance_curve_water.{x}")
#     ggsave(filename, rank_abundance_curve_water,
#            width = 5.5, height = 3, dpi = 600, scale = 1.5)
#   }
# )
# 
# walk(
#   c("svg", "png", "tiff"),
#   \(x) {
#     filename <- glue("figures/rank_abundance_curve_sediment.{x}")
#     ggsave(filename, rank_abundance_curve_sediment,
#            width = 5, height = 3, dpi = 600, scale = 2.5)
#   }
# )




# The end

# Testing ground

pa <- lapply(procrustes_resids, \(dt) {
  dt <- merge(dt, alpha_diversity, all.x = T, by = "sample")
  dt[
    , `:=`(q_ratio = q_WTS/q_WGS,
           N2_ratio = N2_WTS/N2_WGS,
           E_ratio = E_WTS/E_WGS)
  ]
  # numeric_cols <- colnames(dt)[sapply(dt, is.numeric)]
  # dt[, .SD, .SDcols = c("sample", numeric_cols)] %>% 
  #   as.matrix(rownames = "sample")
})

pa_cor <- lapply(pa, \(m) cor(m, method = "spearman"))
















# map2(rank_abundance, names(rank_abundance), \(x, nm) {
#   nm <- ifelse(nm == "WTS", "transcripts", "genes")
#   plot(x, log = "y", main = glue("Rank abundance curves for GH and PL encoding {nm}"))
# })
# 
# test <- lapply(COUNTS, \(m) {
#   m <- split_matrix(m, rows = get_GHPL(), remove_zero_rowsum = TRUE)
#   p <- proportions(m, 2)
#   r <- apply(p, 2, as.rad)
#   
#   l <- list(
#     "proportions" = p,
#     "ranks" = r
#   ) %>% 
#     lapply(as.data.table, keep.rownames = "node")
#   
#   DT <- rbindlist(l, idcol = "measure") %>% 
#     melt(id.vars = c("node", "measure"),
#          measure.vars = patterns("Sed|Filt"),
#          variable.name = "sample",
#          value.name = "abundance") %>%
#     dcast(sample + node ~ measure, value.var = "abundance")
# })
# 
# ggplot(test$WGS, aes(x = ranks, y = proportions)) +
#   geom_point() +
#   facet_wrap(~ sample)

# Addressing comments for the section

# Hence, although freshwater taxa appear able to use a diverse range of
# carbohydrate substrates (on average ???? +/- SD genes per MAG, Fig ????), 
# results indicate a preference for particular substrates over others, and 
# therefore uneven transcription of associated genes.
# richness <- sapply(binary_ghpl, \(i) colSums(i[, 1:30])) %>% 
#   as.data.table(keep.rownames = "sample") %>% 
#   merge(., METADATA[, .(sample, type, salinity)], all.x = TRUE, by = "sample") %>% 
#   .[, ratio := WTS/WGS]
# lm_richness <- glm(ratio ~ salinity + type, family = "Gamma", richness)
# summary(lm_richness)
# car::Anova(lm_richness, type = 2)
# 
# melt(richness, 
#      id.vars = "sample", 
#      measure.vars = c("WGS", "WTS"), 
#      variable.name = "")
# 
# caz_per_mag <- lapply(binary_ghpl, \(i) {
#   as.data.table(i, keep.rownames = "node") %>% 
#     merge(., ANNOTATION[, .(node, bin)], all.x = TRUE, by = "node") %>% 
#     .[, lapply(.SD, sum), by = "bin", .SDcols = patterns("Sed|Filt")] %>% 
#     as.matrix(rownames = "bin")
# })
# 
# lapply(caz_per_mag, \(i) {
#   i <- ifelse(i > 0, i, NA)
#   data.table(
#     "average" = colMeans(i, na.rm = T)[1:30],
#     "stdev" = colVars(i, na.rm = T, std = T)[1:30],
#     "n" = colSums(i > 0, na.rm = T)[1:30]
#   ) %>% 
#     .[, covar := stdev/average]
# })



