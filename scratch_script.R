par(mfrow = c(3, 3))

for (i in 1:ncol(TPM$WTS)) {
  cond <- CPM$WTS[, i] > 0 & COUNTS$WTS[, i] > 0
  y <- log10(CPM$WTS[cond, i])
  x <- log10(COUNTS$WTS[cond, i])
  mod <- lm(y ~ x)
  plot(y = y, x = x, 
       main = colnames(CPM$WTS)[i], 
       sub = str_glue("intercept = {round(coef(mod)[[1]], 2)}; slope = {round(coef(mod)[[2]], 2)}"),
       ylab = "Log CPM", xlab = "Log count", col = "grey50")
  abline(a = 0, b = 1, col = "red")
  abline(mod, col = "blue")
  rm(i, cond, y, x, mode)
}

plot(
  apply(CPM$WTS, 2, \(x) sum(x > 5)) ~ apply(COUNTS$WTS, 2, \(x) sum(x > 0)),
  ylab = "CPM > 1", xlab = "COUNTS > 0"
)

library(edgeR)
CPM1 <- filterByExpr(COUNTS$WTS, group = METADATA$type, 
                     min.total.count = 5, min.count = 1)

dim(COUNTS$WTS[rowSums(COUNTS$WTS) > 0, ])
sum(CPM1)

CPM1 <- ifelse(CPM$WTS > 0, 1, 0)

tapply(as.data.frame(t(COUNTS$WTS)), INDEX = as.factor(METADATA$type[1:30]), FUN = sum)
split(t(COUNTS$WTS), as.factor(METADATA$type[1:30])) %>% head()
head(TPM$WTS)






