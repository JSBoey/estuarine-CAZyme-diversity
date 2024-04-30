genome_size <- ANNOTATION[
  , c("node", "bin")
][
  , contig := str_remove(node, "_\\d+$")
][
  , c("bin", "contig")
] %>% 
  unique() %>% 
  merge(CLEN, all.x = T) %>% 
  group_by(bin) %>% 
  summarise(binsize = sum(contig_length)) %>% 
  left_join(CHECKM, by = join_by("bin" == "Rename")) %>% 
  mutate(
    across(c(Completeness, Contamination, Strain_heterogeneity), ~ .x / 100),
    genome_size = (binsize / Completeness) - (binsize * Contamination)
  )

wts_sum <- as.numeric(SUMMARY$WTS[1, -1])

wts_mag <- as_tibble(COUNTS$WTS, rownames = "node") %>% 
  left_join(ANNOTATION[, c("node", "bin")]) %>% 
  group_by(bin) %>% 
  summarise(across(where(is.numeric), sum)) %>% 
  as.data.frame()
rownames(wts_mag) <- wts_mag$bin
wts_mag <- as.matrix(wts_mag[, -1])

frac_genome_size <- proportions(genome_size$genome_size)
relative_activity <- sweep(sweep(wts_mag, 2, wts_sum, "/"), 1, frac_genome_size, "/")

apply(relative_activity, 2, \(x) sum(x > 1))

test <- relative_activity * ifelse(BINCOV[, 1:30] >= 0.5, 1, 0)

apply(test, 2, \(x) sum(x > 0.5))

# ----

x <- log10(MAG_relativeActivity)
ymax <- max(apply(x, 2, \(j) max(density(j)$y)))
ymin <- min(apply(x, 2, \(j) min(density(j)$y)))
xmax <- max(apply(x, 2, \(j) max(density(j)$x)))
xmin <- min(apply(x, 2, \(j) min(density(j)$x)))
plot(density(x), type = "n", 
     ylim = c(0, ymax) * 1.2, 
     xlim = c(xmin, xmax) * 1.2)
for (j in 1:ncol(x)) {
  lines(density(x[, j]), 
        col = ifelse(
          grepl("Filt", colnames(x)[j]), 
          "cornflowerblue", "firebrick3"
        ))
  rm(j, ymax, ymin, xmax, xmin)
}
lines(density(x), col = "orange", lwd = 6)

xlim <- c(min(x[is.finite(x)]), max(x)) * 1
ylim <- c(0, nrow(x)) * 1
plot(NA, xlim = xlim, ylim = ylim, 
     ylab = "Cumulative frequency", 
     xlab = "Log10 MAG relative activity")
for (j in 1:ncol(x)) {
  sortx <- sort(x[, j])
  breaks <- pretty(sortx, n = 100)
  cf <- cumsum(table(cut(sortx, breaks)))
  color <- ifelse(grepl("Filt", colnames(x)[j]), "cornflowerblue", "firebrick3")
  lines(c(0, cf) ~ breaks, col = color)
  rm(j, sortx, breaks, cf, color, xlim, ylim)
}
abline(v = -3, col = "grey50")

sortx <- sort(x)
breaks <- pretty(sortx, n = 100)
cfx <- cumsum(table(cut(sortx, breaks = breaks)))
plot(c(0, cfx) ~ breaks)
lines(c(0, cfx) ~ breaks)


# ----


bin_dt <- endozyme$beta_mannan$WTS$water[
  , .SD,
  .SDcols = patterns("bin|Sed|Filt")
][
  , lapply(.SD, sum), by = "bin"
]
col_req <- grep("Sed|Filt", colnames(bin_dt), value = T)
binary_dt <- copy(bin_dt)
binary_dt[, (col_req) := lapply(.SD, \(x) ifelse(x > 0, 1, 0)),
          .SDcols = col_req]
mags <- MAG_active[rownames(MAG_active) %in% bin_dt$bin, 
                   colnames(MAG_active) %in% colnames(binary_dt)]
if (length(dim(mags)) < 2)
  mags <- matrix(mags, nrow = 1, ncol = length(mags), 
                 dimnames = list(bin_dt$bin, names(mags)))

mags <- mags[order(match(rownames(mags), bin_dt$bin)), ]
binary_m <- as.matrix(binary_dt, rownames = "bin")
colSums(mags * binary_m)

# ----

sapply(
  c(0.01, seq(0.05, 0.5, by = 0.05)) %>% set_names(., .),
  \(i) {
    apply(MAG_activity, 2, \(j) {
      q <- quantile(j, i)
      sum(j > i)
    })
  }
)

A <- MAG_present[, colnames(MAG_present) %in% overlap_samples]
B <- MAG_relativeActivity
for (j in 1:ncol(B)) {
  V <- B[, j]
  B[, j] <- ifelse(V < quantile(V, 0.1), 0, 1)
}
A - B


library(Rfast)
library(Rfast2)
library(data.table)
library(magrittr)
library(stringr)

b2c <- fread("results/allbins_pred.annotation_table.tsv.gz")[
  , .(bin, node)
][
  , `:=`(
    contig = str_remove(node, "_\\d+$"),
    contig.length = as.numeric(str_replace(node, ".*length_(\\d+)_.*", "\\1"))
  )
][
  , node := NULL
] %>% 
  unique()

bcov_alt <- as.matrix(fread("data/bin_coverage.txt"), rownames = "bin")
bcov_alt <- bcov_alt[order(rownames(bcov_alt)), ]

covfiles <- list.files("data/", pattern = "contig_coverage", full.names = T)
names(covfiles) <- gsub(".*(WGS|WTS).*", "\\1", covfiles)
basecov <- lapply(covfiles, \(s) {
  # Subset by available contigs
  a <- fread(s)[
    rname %chin% b2c$contig, 
    .SD, 
    .SDcols = patterns("rname|bases")
  ]
})

bincov <- lapply(basecov, \(dt) {
  # Sum per bin
  b <- merge(b2c, dt, by.x = "contig", by.y = "rname", all.x = T)
  nm.j <- colnames(b)[sapply(b, is.numeric)]
  b <- b[, lapply(.SD, sum), .SDcols = nm.j, by = "bin"]
})

pbincov <- lapply(bincov, \(dt) {
  nm.j <- grep("covbases", colnames(dt), value = T)
  copy(dt)[, (nm.j) := lapply(.SD, \(x) x * 100 / contig.length), .SDcols = nm.j]
})

pmbincov <- lapply(pbincov, \(dt) {
  m <- as.matrix(copy(dt)[, contig.length := NULL], rownames = "bin")
  m[order(rownames(m)), ]
})

gcor.res <- numeric(nrow(bcov_alt))
names(gcor.res) <- rownames(bcov_alt)

tcor.res <- gcor.res

for (i in 1:nrow(bcov_alt)) {
  gcor.res[i] <- cor(
    bcov_alt[i, ], pmbincov$WGS[i, ], method = "s"
  )
}

for (i in 1:nrow(pmbincov$WTS)) {
  tcor.res[i] <- cor(
    bcov_alt[i, 1:ncol(pmbincov$WTS)], pmbincov$WTS[i, ], method = "s"
  )
}

par(mfrow = c(1, 2))
xlim <- c(
  min(density(gcor.res)$x, density(tcor.res)$x),
  max(density(gcor.res)$x, density(tcor.res)$x)
)
ylim <- c(
  min(density(gcor.res)$y, density(tcor.res)$y),
  max(density(gcor.res)$y, density(tcor.res)$y)
)
plot(density(gcor.res), xlim = xlim, ylim = ylim, col = "cornflowerblue",
     xlab = "Spearman's rho", ylab = "Density")
lines(density(tcor.res), col = "salmon")
plot(tcor.res ~ gcor.res, xlim = c(0, 1), ylim = c(0, 1),
     xlab = "Gene bases",
     ylab = "Transcript bases")
abline(h = 0.5, v = 0.5)



plot(NA, xlim = c(0, 100), ylim = c(0, 100), 
     xlab = "Percent covered (WGS)", ylab = "Percent covered (WTS)")
plot(log2(pmbincov$WTS) ~ log2(pmbincov$WGS[, 1:30]), 
     # xlim = c(0, 100), ylim = c(0, 100), 
     xlab = "Percent covered (WGS)", ylab = "Percent covered (WTS)")
for (i in 1:nrow(pmbincov$WTS)) {
  points(
    x = pmbincov$WGS[i, 1:30],
    y = pmbincov$WTS[i, ]
  )
}
abline(h = 30, v = 30)

sapply(seq(30, 95, 5) %>% set_names(., .), \(i) apply(pmbincov$WGS, 2, \(x) sum(x > i)))

sapply(
  c(.1, seq(0.5, 10, .5)) %>% 
    set_names(., .), 
  \(i) {
    apply(proportions(bcov_alt, 2) * 100, 2, \(x) {
      sum(x > i)
    })
  })
round(proportions(bcov_alt) * 100, 2)

sapply(seq(30, 95, 1) %>% set_names(., .), \(i) apply(pmbincov$WTS, 2, \(x) sum(x > i)))


# Has to have a proportion of relative activity >= 1% and has transcripts mapped to >= 10% of predicted genes
A <- ifelse(
  proportions(MAG_relativeActivity, 2) * 100 >= 1,
  1, 0
)
colSums(A)
B <- as.data.table(ifelse(COUNTS$WTS > 5, 1, 0), keep.rownames = "node")
B <- merge(B, ANNOTATION[, c("bin", "node")], by = "node", all.x = TRUE)
B <- B %>% 
  group_by(bin) %>% 
  summarise(across(where(is.numeric), ~ sum(.x) * 100 / length(.x)))
mB <- as.matrix(as.data.table(B), rownames = "bin")
colSums(mB > 10)
colSums(MAG_present)
colSums(proportions(MAG_relativeActivity, 2) > .025)
colSums(pmbincov$WGS > 20)
colSums(pmbincov$WTS > 10)
apply(COUNTS$WTS, 2, \(x) summary(x[x>0]))
apply(COUNTS$WTS, 2, \(x) sum(x>0))
