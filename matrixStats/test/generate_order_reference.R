#!/usr/bin/env Rscript
# Regenerate the deterministic matrixStats 1.5.0 ordering/weighted oracle.
# Run from the package root: Rscript test/generate_order_reference.R
suppressPackageStartupMessages(library(matrixStats))
stopifnot(as.character(packageVersion("matrixStats")) == "1.5.0")
options(matrixStats.center.onUse = "ignore")
output <- commandArgs(trailingOnly = TRUE)
if (!length(output)) output <- "test/order_reference.txt"
con <- file(output[[1L]], "wt")
count <- 0L
emit <- function(label, value) {
  stopifnot(length(value) == 1L)
  rendered <- if (is.na(value)) "NaN" else sprintf("%.17g", value)
  writeLines(paste(label, rendered), con)
  count <<- count + 1L
}
key <- function(prefix, ...) paste(c(prefix, ...), collapse = "_")
samples <- cbind(c(8, -2, 4, 4, 9, 1, 7, 0, 3),
                 c(NaN, 8, 1, NaN, -1, 6, 6, 2, 3),
                 c(Inf, 3, -Inf, 2, 8, 8, 0, 9, 1), rep(NaN, 9))
probs <- c(0, 0.01, 0.1, 0.25, 1/3, 0.5, 0.7, 0.9, 0.99, 1)
for (s in seq_len(ncol(samples))) {
  x <- samples[, s]
  for (type in 1:9) {
    q <- rowQuantiles(matrix(x, nrow = 1), probs = probs, type = type, na.rm = TRUE)
    for (p in seq_along(probs)) emit(key("q", s, type, p), q[p])
    emit(key("iqr", s, type), iqr(x, na.rm = TRUE, type = type))
  }
  for (na in 0:1) {
    emit(key("median", s, na), rowMedians(matrix(x, nrow = 1), na.rm = as.logical(na)))
    emit(key("mad", s, na), rowMads(matrix(x, nrow = 1), na.rm = as.logical(na)))
    emit(key("mad_center", s, na), rowMads(matrix(x, nrow = 1), center = 2,
                                          constant = 2, na.rm = as.logical(na)))
  }
  for (method in c("max", "min", "average", "first", "last", "dense")) {
    ranks <- rowRanks(matrix(x, nrow = 1), ties.method = method)
    for (i in seq_along(x)) emit(key(paste0("rank_", method), s, i), ranks[i])
  }
  if (!anyNA(x)) {
    for (i in seq_along(x)) emit(key("order", s, i), rowOrderStats(matrix(x, nrow = 1), which = i))
  }
}
xs <- cbind(c(-4, 7, 2, 9, 2, 1), c(NaN, 7, 2, 9, 2, 1), c(NaN, 7, 2, NaN, 2, 1),
            c(Inf, 7, 2, 9, 2, 1), c(-Inf, 7, 2, Inf, 2, 1), rep(NaN, 6))
ws <- cbind(rep(1, 6), c(0, 2, 4, 1, 0, 3), c(-2, 2, -1, 1, 0, 3),
            c(NaN, 2, 4, 1, 0, 3), c(Inf, 2, 4, Inf, 0, 3), rep(0.1, 6),
            rep(0, 6), c(Inf, 2, -Inf, 1, 0, 3))
for (s in seq_len(ncol(xs))) for (t in seq_len(ncol(ws))) for (na in 0:1) {
  x <- xs[, s]
  w <- ws[, t]
  remove <- as.logical(na)
  emit(key("wmean", s, t, na), weightedMean(x, w, na.rm = remove))
  emit(key("wrefine", s, t, na), weightedMean(x, w, na.rm = remove, refine = TRUE))
  emit(key("wmedian", s, t, na), weightedMedian(x, w, na.rm = remove))
  for (method in c("weighted", "min", "max", "mean")) {
    emit(key(paste0("wmedian_", method), s, t, na), weightedMedian(x, w, na.rm = remove, ties = method))
  }
  emit(key("wvar", s, t, na), suppressWarnings(weightedVar(x, w, na.rm = remove)))
  emit(key("wsd", s, t, na), suppressWarnings(weightedSd(x, w, na.rm = remove)))
  emit(key("wmad", s, t, na), weightedMad(x, w, na.rm = remove))
  emit(key("wvar_center", s, t, na), suppressWarnings(weightedVar(x, w, na.rm = remove, center = 2)))
  emit(key("wmad_center", s, t, na), weightedMad(x, w, na.rm = remove, center = 2, constant = 2))
}
for (s in seq_len(ncol(xs))) for (na in 0:1) {
  x <- xs[, s]
  remove <- as.logical(na)
  emit(key("unit_mean", s, na), weightedMean(x, na.rm = remove))
  emit(key("unit_median", s, na), weightedMedian(x, na.rm = remove))
  emit(key("unit_var", s, na), weightedVar(x, na.rm = remove))
  emit(key("unit_mad", s, na), weightedMad(x, na.rm = remove))
}
# Several lengths exercise quicksort partitioning and its unstable order of ties.
for (n in 1:80) {
  i <- seq_len(n)
  x <- as.double((7 * i + 3 * n) %% 11 - 5)
  w <- as.double((3 * i + 2 * n) %% 7)
  emit(key("tie_median", n), weightedMedian(x, w))
  emit(key("tie_mad", n), weightedMad(x, w))
  for (method in c("weighted", "min", "max", "mean")) {
    emit(key(paste0("tie_", method), n), weightedMedian(x, w, ties = method))
  }
}
x <- c(4, 1, 7, 1, 9, 2)
w <- c(1, 2, 0, 3, 4, 2)
emit("unequal_tie_mad", weightedMad(x, w))
close(con)
cat(sprintf("Wrote %d values from matrixStats %s and R %s to %s\n", count,
            packageVersion("matrixStats"), getRversion(), output[[1L]]))
