#!/usr/bin/env Rscript
# TimeIntegers.R: Compare base R vs Randompack for integer draws (ns/value)

library(stats)
library(randompack)

compute_reps <- function(chunk) {
  return(max(1, as.integer(1000000 / chunk)))
}

time_int_range <- function(rng, chunk, bench_time, m, n, use_randompack) {
  reps <- compute_reps(chunk)
  calls <- 0
  t0 <- as.numeric(Sys.time())
  t <- t0
  range_size <- n - m + 1
  use_large <- range_size > 2^31 - 1
  while ((t - t0) < bench_time) {
    for (i in 1:reps) {
      if (use_randompack) {
        rng$int(chunk, m, n)
      }
      else if (use_large) {
        floor(runif(chunk, min = m, max = n + 1))
      }
      else {
        sample.int(as.integer(range_size), chunk, replace = TRUE)
      }
    }
    calls <- calls + reps
    t <- as.numeric(Sys.time())
  }
  total <- calls * chunk
  if (total <= 0) return(NaN)
  return(1e9 * (t - t0) / total)
}

time_perm <- function(rng, bench_time, n, use_randompack) {
  reps <- max(1, as.integer(100000 / n))
  calls <- 0
  t0 <- as.numeric(Sys.time())
  t <- t0
  while ((t - t0) < bench_time) {
    for (i in 1:reps) {
      if (use_randompack) {
        rng$perm(n)
      }
      else {
        sample(n)
      }
    }
    calls <- calls + reps
    t <- as.numeric(Sys.time())
  }
  total <- calls * n
  if (total <= 0) return(NaN)
  return(1e9 * (t - t0) / total)
}

time_sample <- function(rng, bench_time, n, k, use_randompack) {
  reps <- max(1, as.integer(100000 / n))
  calls <- 0
  t0 <- as.numeric(Sys.time())
  t <- t0
  while ((t - t0) < bench_time) {
    for (i in 1:reps) {
      if (use_randompack) {
        rng$sample(n, k)
      }
      else {
        sample(n, k, replace = FALSE)
      }
    }
    calls <- calls + reps
    t <- as.numeric(Sys.time())
  }
  total <- calls * k
  if (total <= 0) return(NaN)
  return(1e9 * (t - t0) / total)
}

warmup <- function(seconds) {
  rp_rng <- randompack_rng()
  t0 <- as.numeric(Sys.time())
  while ((as.numeric(Sys.time()) - t0) < seconds) {
    sample.int(1000, 1024, replace = TRUE)
    sample(1000)
    sample(1000, 100, replace = FALSE)
    rp_rng$int(1024, 1, 1000)
    rp_rng$perm(1000)
    rp_rng$sample(1000, 100)
  }
}

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  bench_time <- 0.2
  chunk <- 4096
  seed <- NULL
  engine <- "x256++simd"

  i <- 1
  while (i <= length(args)) {
    arg <- args[i]
    if (arg == "-h" || arg == "--help") {
      cat("Usage: Rscript TimeIntegers.R [-t sec] [-c chunk] [-s seed] [-e engine]\n")
      cat("  -t sec     time per case in seconds (default: 0.2)\n")
      cat("  -c chunk   number of draws per call (default: 4096)\n")
      cat("  -s seed    fixed seed (default random seed per case)\n")
      cat("  -e engine  Randompack engine (default: x256++simd)\n")
      return(invisible())
    }
    else if (arg == "-t") {
      if (i == length(args)) stop("missing value for -t")
      i <- i + 1
      bench_time <- as.numeric(args[i])
      if (bench_time <= 0) stop("time per case must be positive")
    }
    else if (arg == "-c") {
      if (i == length(args)) stop("missing value for -c")
      i <- i + 1
      chunk <- as.integer(args[i])
      if (chunk <= 0) stop("chunk must be positive")
    }
    else if (arg == "-s") {
      if (i == length(args)) stop("missing value for -s")
      i <- i + 1
      seed <- as.integer(args[i])
    }
    else if (arg == "-e") {
      if (i == length(args)) stop("missing value for -e")
      i <- i + 1
      engine <- args[i]
    }
    else {
      stop(paste("Unknown option:", arg))
    }
    i <- i + 1
  }

  warmup(0.1)

  cat("time per value:   ns/value\n")
  cat(sprintf("bench_time:       %.3f s per case\n", bench_time))
  cat(sprintf("chunk:            %d\n", chunk))
  cat(sprintf("engine:           %s\n", engine))
  cat("\n")
  cat(sprintf("%-18s %10s %11s %8s\n", "Benchmark", "R", "Randompack", "Factor"))

  # int 1-10
  if (!is.null(seed)) set.seed(seed) else set.seed(as.integer(runif(1, 0, 2^31)))
  rp_rng <- randompack_rng(engine)
  rp_rng$seed(if (!is.null(seed)) seed else as.integer(runif(1, 0, 2^31)))
  r_ns <- time_int_range(NULL, chunk, bench_time, 1, 10, FALSE)
  rp_ns <- time_int_range(rp_rng, chunk, bench_time, 1, 10, TRUE)
  factor <- r_ns / rp_ns
  cat(sprintf("%-18s %10.2f %11.2f %8.2f\n", "int 1-10", r_ns, rp_ns, factor))

  # int 1-1e5
  if (!is.null(seed)) set.seed(seed) else set.seed(as.integer(runif(1, 0, 2^31)))
  rp_rng$seed(if (!is.null(seed)) seed else as.integer(runif(1, 0, 2^31)))
  r_ns <- time_int_range(NULL, chunk, bench_time, 1, 100000, FALSE)
  rp_ns <- time_int_range(rp_rng, chunk, bench_time, 1, 100000, TRUE)
  factor <- r_ns / rp_ns
  cat(sprintf("%-18s %10.2f %11.2f %8.2f\n", "int 1-1e5", r_ns, rp_ns, factor))

  # int 1-2e9
  if (!is.null(seed)) set.seed(seed) else set.seed(as.integer(runif(1, 0, 2^31)))
  rp_rng$seed(if (!is.null(seed)) seed else as.integer(runif(1, 0, 2^31)))
  r_ns <- time_int_range(NULL, chunk, bench_time, 1, 2000000000, FALSE)
  rp_ns <- time_int_range(rp_rng, chunk, bench_time, 1, 2000000000, TRUE)
  factor <- r_ns / rp_ns
  cat(sprintf("%-18s %10.2f %11.2f %8.2f\n", "int 1-2e9", r_ns, rp_ns, factor))

  # perm 100
  if (!is.null(seed)) set.seed(seed) else set.seed(as.integer(runif(1, 0, 2^31)))
  rp_rng$seed(if (!is.null(seed)) seed else as.integer(runif(1, 0, 2^31)))
  r_ns <- time_perm(NULL, bench_time, 100, FALSE)
  rp_ns <- time_perm(rp_rng, bench_time, 100, TRUE)
  factor <- r_ns / rp_ns
  cat(sprintf("%-18s %10.2f %11.2f %8.2f\n", "perm 100", r_ns, rp_ns, factor))

  # perm 100000
  if (!is.null(seed)) set.seed(seed) else set.seed(as.integer(runif(1, 0, 2^31)))
  rp_rng$seed(if (!is.null(seed)) seed else as.integer(runif(1, 0, 2^31)))
  r_ns <- time_perm(NULL, bench_time, 100000, FALSE)
  rp_ns <- time_perm(rp_rng, bench_time, 100000, TRUE)
  factor <- r_ns / rp_ns
  cat(sprintf("%-18s %10.2f %11.2f %8.2f\n", "perm 100000", r_ns, rp_ns, factor))

  # sample 20/1000
  if (!is.null(seed)) set.seed(seed) else set.seed(as.integer(runif(1, 0, 2^31)))
  rp_rng$seed(if (!is.null(seed)) seed else as.integer(runif(1, 0, 2^31)))
  r_ns <- time_sample(NULL, bench_time, 1000, 20, FALSE)
  rp_ns <- time_sample(rp_rng, bench_time, 1000, 20, TRUE)
  factor <- r_ns / rp_ns
  cat(sprintf("%-18s %10.2f %11.2f %8.2f\n", "sample 20/1000", r_ns, rp_ns, factor))

  # sample 500/1000
  if (!is.null(seed)) set.seed(seed) else set.seed(as.integer(runif(1, 0, 2^31)))
  rp_rng$seed(if (!is.null(seed)) seed else as.integer(runif(1, 0, 2^31)))
  r_ns <- time_sample(NULL, bench_time, 1000, 500, FALSE)
  rp_ns <- time_sample(rp_rng, bench_time, 1000, 500, TRUE)
  factor <- r_ns / rp_ns
  cat(sprintf("%-18s %10.2f %11.2f %8.2f\n", "sample 500/1000", r_ns, rp_ns, factor))

  # sample 980/1000
  if (!is.null(seed)) set.seed(seed) else set.seed(as.integer(runif(1, 0, 2^31)))
  rp_rng$seed(if (!is.null(seed)) seed else as.integer(runif(1, 0, 2^31)))
  r_ns <- time_sample(NULL, bench_time, 1000, 980, FALSE)
  rp_ns <- time_sample(rp_rng, bench_time, 1000, 980, TRUE)
  factor <- r_ns / rp_ns
  cat(sprintf("%-18s %10.2f %11.2f %8.2f\n", "sample 980/1000", r_ns, rp_ns, factor))
}

main()
