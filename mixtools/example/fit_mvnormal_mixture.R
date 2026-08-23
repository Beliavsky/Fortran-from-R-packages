#!/usr/bin/env Rscript

# Fit a two-component multivariate-normal mixture with the R mixtools package.
# The first optional argument is the observation file and the second is the
# report file.

overall_start <- Sys.time()
suppressPackageStartupMessages(library(mixtools))

script_directory <- function() {
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) == 0L) {
    return(getwd())
  }
  dirname(normalizePath(sub("^--file=", "", file_arg[[1L]])))
}

args <- commandArgs(trailingOnly = TRUE)
data_file <- if (length(args) >= 1L) {
  args[[1L]]
} else {
  file.path(script_directory(), "mvnormal_mixture_data.txt")
}
report_file <- if (length(args) >= 2L) {
  args[[2L]]
} else {
  file.path(script_directory(), "mvnormal_mixture_fit_r.txt")
}

read_start <- Sys.time()
connection <- file(data_file, open = "rt")
dimensions <- scan(connection, what = integer(), nmax = 2L, quiet = TRUE)
values <- scan(connection, what = double(), quiet = TRUE)
close(connection)
read_seconds <- as.numeric(difftime(Sys.time(), read_start, units = "secs"))

if (length(dimensions) != 2L || any(dimensions <= 0L)) {
  stop("The first line must contain positive row and column counts.")
}
if (length(values) != prod(dimensions)) {
  stop("The number of observations does not match the file header.")
}

x <- matrix(values, nrow = dimensions[[1L]], ncol = dimensions[[2L]], byrow = TRUE)
k <- 2L
initial_rows <- c(1L, 1L + nrow(x) %/% k)
pooled_covariance <- stats::cov(x)

fit_start <- Sys.time()
fit <- mixtools::mvnormalmixEM(
  x,
  lambda = rep(1 / k, k),
  mu = lapply(initial_rows, function(i) x[i, ]),
  sigma = replicate(k, pooled_covariance, simplify = FALSE),
  k = k,
  epsilon = 1e-8,
  maxit = 1000,
  verb = FALSE
)
fit_seconds <- as.numeric(difftime(Sys.time(), fit_start, units = "secs"))
overall_seconds <- as.numeric(difftime(Sys.time(), overall_start, units = "secs"))

mean_matrix <- do.call(cbind, fit$mu)
component_order <- order(mean_matrix[1L, ])

report <- capture.output({
  cat("R mixtools mvnormalmixEM fit\n")
  cat("observations:", nrow(x), " dimensions:", ncol(x), "\n")
  cat("loglik:", sprintf("%.12g", fit$loglik), "\n")
  cat("iterations:", length(fit$all.loglik) - 1L, "\n")
  cat("read seconds:", sprintf("%.6f", read_seconds), "\n")
  cat("fit seconds:", sprintf("%.6f", fit_seconds), "\n")
  cat("overall seconds:", sprintf("%.6f", overall_seconds), "\n")
  for (position in seq_along(component_order)) {
    component <- component_order[[position]]
    cat("component", position, "\n")
    cat("  weight:", sprintf("%.12g", fit$lambda[[component]]), "\n")
    cat("  mean:", paste(sprintf("%.12g", fit$mu[[component]]), collapse = " "), "\n")
    cat("  covariance:\n")
    for (row in seq_len(nrow(fit$sigma[[component]]))) {
      cat(
        paste(sprintf("%.12g", fit$sigma[[component]][row, ]), collapse = " "),
        "\n"
      )
    }
  }
})

writeLines(report, report_file)
cat(paste(report, collapse = "\n"), "\n")
cat("Wrote fit report to", report_file, "\n")
