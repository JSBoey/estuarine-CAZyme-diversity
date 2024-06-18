# Explorations

# Where does coverage, genome cover, and transcript mapped match per bin?

# Per sample thresholds:
threshold <- c(
  "bcov" = 0.5,   # Coverage depth 
  "gcov" = 0.1, # Genome cover
  "tmap" = 100  # Transcript mapped
) 
log_threshold <- log10(threshold)

bcov <- as.data.table(proportions(BINCOV$WGS.numreads, 2), keep.rownames = "bin")
gcov <- get_breadth_coverage(CONTIGCOV$WGS.covbases) %>% 
  as.data.table(keep.rownames = "bin")
if (nrow(COUNTS$WTS) == nrow(ANNOTATION)) {
  tmap <- merge(as.data.table(COUNTS$WTS, keep.rownames = "node"),
                ANNOTATION[, c("node", "bin")],
                by = "node")
  tmap <- tmap[, lapply(.SD, sum), by = "bin", .SDcols = patterns("Filt|Sed")]
}

for (i in c("bcov", "gcov", "tmap")) {
  mdt <- melt(get(i),
              id.vars = "bin",
              variable.name = "sample",
              value.name = i)
  assign(i, mdt)
  rm(i, mdt)
}

fmerge <- function(a, b) merge(a, b, by = c("bin", "sample"), all.x = TRUE)
all_dt <- Reduce(fmerge, list(gcov, bcov, tmap))
all_dt <- all_dt[, type := gsub("(Filt|Sed).*", "\\1", sample)]
log_dt <- all_dt[
  , c("bcov", "gcov", "tmap") := Map(\(j) ifelse(j > 0, log10(j), 0), .SD),
  .SDcols = c("bcov", "gcov", "tmap")
][
  , tmap_cut := tmap > log_threshold["tmap"]
]

ggplot(log_dt) +
  geom_point(aes(x = gcov, 
                 y = bcov, 
                 colour = tmap_cut,
                 alpha = tmap), 
             size = 1) +
  # geom_hline(yintercept = log_threshold["bcov"], 
  #            linetype = 2, 
  #            colour = "orangered") +
  geom_vline(xintercept = log_threshold["gcov"], 
             linetype = 2, 
             colour = "orangered") +
  scale_colour_viridis_d() +
  facet_wrap(~ sample) +
  labs(x = "% Genome covered",
       y = "Coverage depth",
       colour = "Transcripts mapped to genome",
       caption = "All values are log10 transformed") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "bottom")

