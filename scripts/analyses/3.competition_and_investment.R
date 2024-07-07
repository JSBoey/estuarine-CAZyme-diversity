# Resource access: Functional redundancy and competition
library(ggdendro)

cols <- colnames(COUNTS$WTS)

senv <- c("Filt", "Sed") %>% 
  setNames(nm = .)

fig.fmt <- c('png', 'tiff', 'svg')

present_bins <- ifelse(
  get_breadth_coverage(CONTIGCOV$WGS.covbases) < 0.1,
  0, 1
)

if (all(rownames(COUNTS$WTS) == names(get_gene_length()))) {
  tpm <- calc_TPM(threshold_matrix(COUNTS$WTS, 5), get_gene_length())
}

## Zero genes with < 5 transcript reads
# tpm[COUNTS$WTS < 5] <- 0

# Unnest CAZYME
cdt <- CAZYME[
  grepl("(GH|PL)", FAMILY), c("NODE", "FAMILY", "EC")
][
  , .(NODE = rep(NODE, lengths(EC)),
      FAMILY = rep(FAMILY, lengths(EC)),
      EC = unlist(EC))
][
  , EC := ifelse(is.na(EC), "", EC)
][
  , cazyme := paste0(FAMILY, " | ", EC)
]
names(cdt) <- tolower(names(cdt))
cdt <- merge(cdt, ANNOTATION[, c("node", "bin")], by = "node")

# Per MAG investment (normalised by genome coverage)
inv <- calc_per_bin_proportions(tpm,
                                ANNOTATION[, c("node", "bin")],
                                "node",
                                "bin",
                                BINCOV$WGS.numreads[, cols])

# Functional profile
# Tian et al. (2020) doi:10.1038/s41467-020-19940-1
# Aggregated by CAZyme family and EC across samples (whole)
# Hclust ordering by Euclidean distance
G <- cdt %$% table(bin, cazyme)
pv <- proportions(BINCOV$WGS.numreads * present_bins, 2)
pv <- pv[match(rownames(G), rownames(pv)), ]

all(rownames(pv) == rownames(G))

k <- 1 / colSums(t(G) %*% pv)
fv <- t(k * (t(pv) %*% G))

fv.env <- c("Filt", "Sed") %>% 
  setNames(., .) %>% 
  lapply(\(ev) {
  m <- split_matrix(fv, columns = ev, remove_zero_rowsum = TRUE)
  d <- vegdist(m, "euclidean")
  h <- hclust(d, method = "complete")
  lvl <- h$labels[h$order]
  dt <- as.data.table(m, keep.rownames = "label") %>% 
    melt(measure.vars = colnames(m),
         variable.name = "sample",
         value.name = "rel_f")
  dt_top <- unique(dt[order(-rel_f), head(.SD, 10), by = "sample"]$label)
  # dt[, label := factor(label, levels = lvl)]
  dt[label %in% dt_top]
})

fv.col <- RColorBrewer::brewer.pal(n = 8, name = "Set1")
fv.plt <- lapply(fv.env, \(dt) {
  cazymes <- unique(dt$label)
  fill <- rep(fv.col, length.out = length(cazymes)) %>% 
    setNames(cazymes)
  ggplot(dt) +
    geom_col(aes(x = sample, y = rel_f, fill = label)) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0), 
                       limits = c(0, 0.65)) +
    scale_fill_manual(values = fill,
                      name = "CAZyme") +
    theme_bw() +
    theme(legend.position = "right",
          axis.title = element_blank(),
          panel.grid = element_blank(),
          axis.text = element_text(size = 7),
          axis.text.x = element_text(angle = 90,
                                     hjust = 1, 
                                     vjust = 0.5))
})

fv.pltwp <- wrap_plots(fv.plt, widths = c(7, 20)) +
  plot_annotation(caption = "CAZyme gene profile per sample")

walk(fig.fmt, \(fmt) {
  ggsave(filename = glue("figures/cazyme_gene_profile.{fmt}"),
         plot = fv.pltwp, 
         width = 297, 
         height = 210, 
         units = 'mm', 
         scale = 1)
})


## Bin to function matrix
fv.bmt <- lapply(fv.env, \(dt) {
  caz <- unique(dt$label)
  scol <- levels(dt$sample)
  
  d1 <- as.data.table(present_bins[, scol], keep.rownames = "bin") %>% 
    melt(id.vars = 'bin',
         variable.name = 'sample',
         value.name = 'presence') %>% 
    .[presence > 0]
  
  d2 <- as.data.table(G, keep.rownames = 'bin')[N > 0 & cazyme %in% caz]
  
  DT <- merge(d1, d2, by = "bin", all = TRUE, allow.cartesian = TRUE)
  DT <- DT[
    !is.na(presence)
  ][
    , .(nbin = .N), by = .(cazyme, sample)
  ]
})

fv.bmt.plt <- lapply(fv.bmt, \(dt) {
  ggplot(dt) +
    geom_tile(aes(x = sample, 
                  y = cazyme,
                  fill = nbin),
              alpha = 0.4) +
    geom_text(aes(x = sample,
                  y = cazyme,
                  label = nbin,
                  colour = nbin),
              size = 2) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    scale_colour_gradient(low = "navy", high = "red") +
    scale_fill_viridis_c() +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_line(colour = "black"),
          aspect.ratio = 1,
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          text = element_text(size = 7),
          axis.text.x = element_text(angle = 90, 
                                     hjust = 1, 
                                     vjust = 0.5))
})
ggsave(
  "sandbox/bin2cazyme2sample.matrix.tiff",
  wrap_plots(fv.bmt.plt), 
  width = 297, height = 210, units = "mm"
)

# Spatial distribution of most transcriptionally active CAZyme genes
# Heuristic: Substrate favoured at the community level
# Aggregate sum at family-EC level per MAG
czdst <- merge(cdt, 
               as.data.table(tpm, keep.rownames = "node"),
               by = "node") %>% 
  .[, lapply(.SD, sum), by = "cazyme", .SDcols = cols]

czdst.top <- lapply(senv, \(ev) {
  # Convert data
  cols.ev <- grep(ev, cols, value = TRUE)
  m <- czdst[, .SD, .SDcols = c("cazyme", cols.ev)] %>% 
    as.matrix(rownames = "cazyme")
  
  # Subset function
  ftop <- function(x, n) {
    sort_x <- sort(x[x > 0], decreasing = TRUE)
    names(head(sort_x, n))
  }
  top <- apply(m, 2, \(j) ftop(j, 10))
  
  # Unwrap
  top <- if (is.list(top))
    unwrap(top)
  else
    unique(as.vector(top))
  
  # Subset
  split_matrix(m[top, ], remove_zero_rowsum = TRUE)
})

czdst.fill_max <- max(sapply(czdst.top, max))

czdst.plt <- lapply(czdst.top, \(m) {
  # Transform
  m <- ifelse(m > 0, log2(m), m)
  fill_max <- log2(czdst.fill_max)
  
  # Ordering
  d <- dist(m)
  h <- dendro_data(hclust(d))
  
  # Convert data
  dt <- as.data.table(m, keep.rownames = 'cazyme') %>% 
    melt(measure.vars = colnames(m),
         variable.name = 'sample',
         value.name = 'stpm') %>% 
    .[, cazyme := factor(cazyme, levels = label(h)$label)]
  
  # Draw hclust
  p1 <- ggplot(segment(h)) +
    geom_segment(aes(x = -y,
                     y = x,
                     xend = -yend,
                     yend = xend)) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0.5)) +
    theme_dendro() +
    theme(plot.margin = unit(c(0, 0, 0, 5), "mm"))
  
  # Draw heatmap
  p2 <- ggplot(dt) +
    geom_tile(aes(x = sample,
                  y = cazyme,
                  fill = stpm)) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0),
                     position = 'right') +
    scale_fill_viridis_c(name = "log2 TPM",
                         limits = c(0, fill_max)) +
    theme_bw() +
    theme(plot.margin = unit(rep(0, 4), "mm"),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_text(size = 7),
          axis.text.x = element_text(angle = 90, 
                                     hjust = 1, 
                                     vjust = 0.5))
  
  p1 + p2
  
})

wrap_plots(czdst.plt, 
           widths = c(9, 21),
           guides = 'collect')

czdst.pltwp <- wrap_plots(czdst.plt, 
                          widths = c(9, 21),
                          guides = 'collect') +
  plot_annotation(caption = "Top 10 (per sample) transcriptionally active CAZymes")

walk(fig.fmt, \(fmt) {
  ggsave(filename = glue("figures/top10_transcribed_cazymes.{fmt}"),
         plot = czdst.pltwp, 
         width = 297, 
         height = 210, 
         units = 'mm', 
         scale = 1)
})


# Pick out terrestrial/riverine (pectin, cellulose, xylan), all (chitin), marine (laminarin, alginate, and ulvan) degrading CAZymes for functional profile and transcription spatial distribution


# Simpson's D of CAZyme transcripts 
# Aggregated at family-EC level per MAG
# Heuristic: Competition for substrate
fcomp <- function(x) {
  if (all(x == 0)) return(0)
  diversity(x, index = "invsimpson")
}
comp <- merge(cdt, 
              as.data.table(tpm, keep.rownames = "node"),
              by = "node") %>% 
  .[, lapply(.SD, fcomp), by = "cazyme", .SDcols = cols]

comp.top <- lapply(senv, \(ev) {
  # Convert data
  cols.ev <- grep(ev, cols, value = TRUE)
  m <- comp[, .SD, .SDcols = c("cazyme", cols.ev)] %>% 
    as.matrix(rownames = "cazyme")
  
  # Subset function
  ftop <- function(x, n) {
    sort_x <- sort(x[x > 0], decreasing = TRUE)
    names(head(sort_x, n))
  }
  top <- apply(m, 2, \(j) ftop(j, 10))
  
  # Unwrap
  top <- if (is.list(top))
    unwrap(top)
  else
    unique(as.vector(top))
  
  # Subset
  split_matrix(m[top, ], remove_zero_rowsum = TRUE)
})

comp.fill_max <- max(sapply(comp.top, max))

comp.plt <- lapply(comp.top, \(m) {
  # Transform
  m <- ifelse(m > 0, log2(m) + 1, 0)
  fill_max <- log2(comp.fill_max) + 1
  
  # Ordering
  d <- dist(m)
  h <- dendro_data(hclust(d))
  
  # Convert data
  dt <- as.data.table(m, keep.rownames = 'cazyme') %>% 
    melt(measure.vars = colnames(m),
         variable.name = 'sample',
         value.name = 'stpm') %>% 
    .[, cazyme := factor(cazyme, levels = label(h)$label)]
  
  # Draw hclust
  p1 <- ggplot(segment(h)) +
    geom_segment(aes(x = -y,
                     y = x,
                     xend = -yend,
                     yend = xend)) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0.5)) +
    theme_dendro() +
    theme(plot.margin = unit(c(0, 0, 0, 5), "mm"))
  
  # Draw heatmap
  p2 <- ggplot(dt) +
    geom_tile(aes(x = sample,
                  y = cazyme,
                  fill = stpm)) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0),
                     position = 'right') +
    scale_fill_viridis_c(name = "log2 TPM",
                         limits = c(0, fill_max)) +
    theme_minimal() +
    theme(plot.margin = unit(rep(0, 4), "mm"),
          panel.grid = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_text(size = 7),
          axis.text.x = element_text(angle = 90, 
                                     hjust = 1, 
                                     vjust = 0.5))
  
  p1 + p2
  
})

wrap_plots(comp.plt, 
           widths = c(9, 21),
           guides = 'collect')

comp.pltwp <- wrap_plots(comp.plt, 
                          widths = c(9, 21),
                          guides = 'collect') +
  plot_annotation(caption = "Top 10 (per sample) CAZymes by Simpson's diversity")

walk(fig.fmt, \(fmt) {
  ggsave(filename = glue("figures/top10_competitive_cazymes.{fmt}"),
         plot = comp.pltwp, 
         width = 297, 
         height = 210, 
         units = 'mm', 
         scale = 1)
})

# Which populations are most invested in these substrates (see above)?
# Proportion TPM normalised by MAG coverage
inv.cazyme <- lapply(senv, \(ev) {
  Reduce(union, 
         list(unique(fv.env[[ev]]$label),
              rownames(czdst.top[[ev]]),
              rownames(comp.top[[ev]])
              )
         ) %>% 
    unwrap()
})
# Remove peptidoglycan
peptidoglycan <- paste0("GH", c(19, 23, 73, 102, 103, 171)) %>% 
  paste(collapse = "|")
inv.cazyme <- lapply(inv.cazyme, \(x) grep(peptidoglycan, x, value = T, invert = T))

# Only endolytic
inv.cazyme <- c(
  "GH3 | 3.2.1.6",
  "GH5 | 3.2.1.4",
  "GH10 | 3.2.1.8",
  "GH13 | 3.2.1.1",
  "GH13 | 3.2.1.68",
  "GH16 | 3.2.1.6",
  "GH16 | 3.2.1.39",
  "GH16 | 3.2.1.73",
  "GH17 | 3.2.1.39",
  "GH43 | 3.2.1.99",
  "GH53 | 3.2.1.89",
  "GH57 | 3.2.1.1",
  "GH74 | 3.2.1.151",
  "GH81 | ",
  "GH148 | 3.2.1.4",
  "GH148 | 3.2.1.-",
  "GH158 | ",
  "PL1 | 4.2.2.2",
  "PL3 | 4.2.2.2",
  "PL6 | 4.2.2.-",
  "PL6 | 4.2.2.11"
)

cinv <- merge(cdt,
              as.data.table(inv, keep.rownames = "node"),
              by = 'node') %>% 
  .[, lapply(.SD, sum), by = c("bin", "cazyme"), .SDcols = cols]

cinv.top <- lapply(senv, \(ev) {
  # Convert data
  cols.ev <- grep(ev, cols, value = TRUE)
  m <- cinv[
    , rlab := paste0(bin, " | ", cazyme)
  ][
    , .SD, .SDcols = c("rlab", cols.ev)
  ] %>% 
    as.matrix(rownames = "rlab")
  
  # Subset function
  ftop <- function(x, n) {
    sort_x <- sort(x[x > 0], decreasing = TRUE)
    names(head(sort_x, n))
  }
  top <- apply(m, 2, \(j) ftop(j, 10))
  
  # Unwrap
  top <- if (is.list(top))
    unwrap(top)
  else
    unique(as.vector(top))
  
  # Subset
  split_matrix(m[top, ], remove_zero_rowsum = TRUE)
})

cinv.top <- lapply(senv, \(ev) {
  # Convert data
  cols.ev <- grep(ev, cols, value = TRUE)
  m <- cinv[
    cazyme %in% inv.cazyme
  ][
    , rlab := paste0(bin, " | ", cazyme)
  ][
    , .SD, .SDcols = c("rlab", cols.ev)
  ] %>% 
    as.matrix(rownames = "rlab")

  # Subset function
  ftop <- function(x, n) {
    sort_x <- sort(x[x > 0], decreasing = TRUE)
    names(head(sort_x, n))
  }
  top <- apply(m, 2, \(j) ftop(j, 10))

  # Unwrap
  top <- if (is.list(top))
    unwrap(top)
  else
    unique(as.vector(top))

  # Subset
  split_matrix(m[top, ], remove_zero_rowsum = TRUE)
})

cinv.fill_max <- max(sapply(cinv.top, max))

cinv.plt <- lapply(cinv.top, \(m) {
  # Transform
  m <- ifelse(m > 0, log2(m) + 1, 0)
  fill_max <- log2(cinv.fill_max) + 1
  
  # Ordering
  d <- dist(m)
  h <- dendro_data(hclust(d))
  
  # Convert data
  dt <- as.data.table(m, keep.rownames = 'cazyme') %>% 
    melt(measure.vars = colnames(m),
         variable.name = 'sample',
         value.name = 'stpm') %>% 
    .[, cazyme := factor(cazyme, levels = label(h)$label)]
  
  # Draw hclust
  p1 <- ggplot(segment(h)) +
    geom_segment(aes(x = -y,
                     y = x,
                     xend = -yend,
                     yend = xend)) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0.5)) +
    theme_dendro() +
    theme(plot.margin = unit(c(0, 0, 0, 5), "mm"))
  
  # Draw heatmap
  p2 <- ggplot(dt) +
    geom_tile(aes(x = sample,
                  y = cazyme,
                  fill = stpm)) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0),
                     position = 'right') +
    scale_fill_viridis_c(name = "log2 percent TPM",
                         limits = c(0, fill_max)) +
    theme_minimal() +
    theme(plot.margin = unit(rep(0, 4), "mm"),
          panel.grid = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_text(size = 7),
          axis.text.x = element_text(angle = 90, 
                                     hjust = 1, 
                                     vjust = 0.5))
  
  p1 + p2
  
})

wrap_plots(cinv.plt, 
           widths = c(9, 21),
           guides = 'collect')

cinv.pltwp <- wrap_plots(cinv.plt, 
                         widths = c(9, 21),
                         guides = 'collect') +
  plot_annotation(caption = "Top 10 MAG investment (coverage corrected) in highly active endolytic CAZyme-encoding genes for highly competitive substrates.")

walk(fig.fmt, \(fmt) {
  ggsave(filename = glue("figures/top10_invested_endo_cazymes.{fmt}"),
         plot = cinv.pltwp, 
         width = 297, 
         height = 210, 
         units = 'mm', 
         scale = 1)
})
