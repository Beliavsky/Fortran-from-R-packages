#!/usr/bin/env Rscript

# Generate a reproducible two-component mixture of bivariate normal samples.
# By default the observations are written beside this script. An alternative
# output path may be supplied as the first command-line argument.

suppressPackageStartupMessages(library(mixtools))

script_directory <- function() {
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) == 0L) {
    return(getwd())
  }
  dirname(normalizePath(sub("^--file=", "", file_arg[[1L]])))
}

args <- commandArgs(trailingOnly = TRUE)
output_file <- if (length(args) >= 1L) {
  args[[1L]]
} else {
  file.path(script_directory(), "mvnormal_mixture_data.txt")
}

set.seed(20260823)

component_sizes <- c(240L, 360L)
means <- list(c(-2.0, -1.0), c(2.0, 2.0))
covariances <- list(
  matrix(c(1.0, 0.45, 0.45, 0.70), nrow = 2L, byrow = TRUE),
  matrix(c(0.80, -0.30, -0.30, 1.20), nrow = 2L, byrow = TRUE)
)

# Keeping the component blocks together makes the deterministic initial means
# used by the Fortran example start in different components. The component
# labels themselves are not written to the data file.
observations <- rbind(
  mixtools::rmvnorm(component_sizes[[1L]], means[[1L]], covariances[[1L]]),
  mixtools::rmvnorm(component_sizes[[2L]], means[[2L]], covariances[[2L]])
)

connection <- file(output_file, open = "wt")
on.exit(close(connection), add = TRUE)
writeLines(sprintf("%d %d", nrow(observations), ncol(observations)), connection)
writeLines(
  apply(observations, 1L, function(row) {
    paste(sprintf("%.17g", row), collapse = " ")
  }),
  connection
)

cat("Wrote", nrow(observations), "observations to", output_file, "\n")
cat("True weights: 0.4 0.6\n")
cat("True means: (-2, -1) and (2, 2)\n")
