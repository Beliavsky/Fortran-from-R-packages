# TimeDist.R
# Compare base R, dqrng, and randompack distributions (ns/value)

if (!requireNamespace("randompack", quietly = TRUE)) {
  stop("Package 'randompack' is required. Install it before running.")
}
if (!requireNamespace("dqrng", quietly = TRUE)) {
  stop("Package 'dqrng' is required. Install it before running.")
}

args <- commandArgs(trailingOnly = TRUE)
chunk <- 4096L
engine <- ""
bitexact <- FALSE
bench_time <- 0.2
seed <- NA_integer_

print_help <- function() {
  cat("Usage: Rscript TimeDist.R [-h] [-b] [-c chunk] [-e engine] [-t sec] [-s seed] [engine]\n")
  cat("  -h         Show this help message\n")
  cat("  -b         Enable bitexact mode\n")
  cat("  -c chunk   Set chunk size (default 4096)\n")
  cat("  -e engine  Set RNG engine (default x256++simd)\n")
  cat("  -t sec     Set time per case in seconds (default 0.2)\n")
  cat("  -s seed    Set fixed seed (default random seed per case)\n")
}

i <- 1L
while (i <= length(args)) {
  arg <- args[[i]]
  if (arg == "-h") {
    print_help()
    quit(save="no", status=0)
  } else if (arg == "-b") {
    bitexact <- TRUE
  } else if (arg == "-c") {
    if (i == length(args)) {
      stop("Use -c <chunk> with a following positive integer.")
    }
    i <- i + 1L
    chunk <- as.integer(args[[i]])
    if (is.na(chunk) || chunk <= 0L) {
      stop("Chunk size must be a positive integer.")
    }
  } else if (arg == "-e") {
    if (i == length(args)) {
      stop("Use -e <engine> with a following engine name.")
    }
    i <- i + 1L
    engine <- args[[i]]
  } else if (arg == "-t") {
    if (i == length(args)) {
      stop("Use -t <sec> with a following positive number.")
    }
    i <- i + 1L
    bench_time <- as.numeric(args[[i]])
    if (is.na(bench_time) || bench_time <= 0) {
      stop("Time per case must be a positive number.")
    }
  } else if (arg == "-s") {
    if (i == length(args)) {
      stop("Use -s <seed> with a following integer.")
    }
    i <- i + 1L
    seed <- as.integer(args[[i]])
    if (is.na(seed)) {
      stop("Seed must be an integer.")
    }
  } else if (startsWith(arg, "-")) {
    stop(sprintf("Unknown option: %s", arg))
  } else if (!nzchar(engine)) {
    engine <- arg
  } else {
    stop(sprintf("Unexpected argument: %s", arg))
  }
  i <- i + 1L
}

if (!nzchar(engine)) {
  engine <- "x256++simd"
}

rng <- randompack::randompack_rng(engine=engine, bitexact=bitexact)

reps = max(1, floor(1e6 / chunk))

cat(sprintf("Platform:  %s\n", R.version$platform))
cat(sprintf("Engine:    %s\n", engine))
cat(sprintf("Chunk:     %d\n", chunk))
cat(sprintf("Time/case: %.3f s\n", bench_time))

time_dist = function(f, chunk, reps, bench_time) {
  calls = 0
  t0 = proc.time()[["elapsed"]]
  while ((proc.time()[["elapsed"]] - t0) < bench_time) {
    for (i in 1:reps) {
      x = f()
      sink = x[chunk]
    }
    calls = calls + reps
  }
  t1 = proc.time()[["elapsed"]]
  1e9 * (t1 - t0) / (calls * chunk)
}

warmup = function(seconds) {
  t0 <- proc.time()[["elapsed"]]
  while ((proc.time()[["elapsed"]] - t0) < seconds) {
    rng$unif(1000)
    runif(1000)
  }
  proc.time()[["elapsed"]] - t0
}

warm = warmup(0.1)
cat(sprintf("Warmup:    %.3f s\n\n", warm))

dist_names = c(
  "u01",
  "unif(2,5)",
  "norm",
  "normal(2,3)",
  "exp(1)",
  "exp(2)",
  "lognormal(0,1)",
  "skew_norm(0,1,2)",
  "gumbel(0,1)",
  "pareto(1,2)"
)
dist_width = max(nchar(c("DISTRIBUTION", dist_names)))
header_fmt = paste0("%-", dist_width, "s %10s %11s %10s %8s %8s\n")
row_fmt = paste0("%-", dist_width, "s %10.2f %11.2f %10s %8.2f %8s\n")
rp_only_fmt = paste0("%-", dist_width, "s %10s %11.2f %10s %8s %8s\n")

cat(sprintf(header_fmt, "DISTRIBUTION", "BASE-R", "RANDOMPACK", "DQRNG",
            "FACTOR-B", "FACTOR-D"))

old_kind <- RNGkind()
draw_seed <- function() {
  if (!is.na(seed)) {
    return(seed)
  }
  sample.int(.Machine$integer.max, 1L)
}

restore_rng <- function(case_seed) {
  RNGkind(kind=old_kind[1], normal.kind=old_kind[2], sample.kind=old_kind[3])
  set.seed(case_seed)
  rng$seed(case_seed)
}

run_case = function(name, f_base, f_rp, f_dqrng=NULL, use_dqset=FALSE) {
  case_seed = draw_seed()
  restore_rng(case_seed)
  rp_ns = time_dist(f_rp, chunk, reps, bench_time)
  restore_rng(case_seed)
  base_ns = time_dist(f_base, chunk, reps, bench_time)
  factor_b = base_ns / rp_ns
  if (is.null(f_dqrng)) {
    dqrng_str = sprintf("%10s", "")
    factor_d_str = sprintf("%8s", "")
  } else {
    if (use_dqset) dqrng::dqset.seed(case_seed)
    dqrng_ns = time_dist(f_dqrng, chunk, reps, bench_time)
    dqrng_str = sprintf("%10.2f", dqrng_ns)
    factor_d_str = sprintf("%8.2f", dqrng_ns / rp_ns)
  }
  cat(sprintf(row_fmt,
              name, base_ns, rp_ns, dqrng_str, factor_b, factor_d_str))
}

run_case("u01",
         function() runif(chunk),
         function() rng$unif(chunk),
         function() dqrng::dqrunif(chunk),
         use_dqset=TRUE)

run_case("unif(2,5)",
         function() runif(chunk, min=2, max=5),
         function() rng$unif(chunk, 2, 5),
         function() dqrng::dqrunif(chunk, min=2, max=5),
         use_dqset=TRUE)

run_case("norm",
         function() rnorm(chunk),
         function() rng$normal(chunk),
         function() dqrng::dqrnorm(chunk),
         use_dqset=TRUE)

run_case("normal(2,3)",
         function() rnorm(chunk, mean=2, sd=3),
         function() rng$normal(chunk, 2, 3),
         function() dqrng::dqrnorm(chunk, mean=2, sd=3),
         use_dqset=TRUE)

run_case("exp(1)",
         function() rexp(chunk, rate=1),
         function() rng$exp(chunk, 1))

run_case("exp(2)",
         function() rexp(chunk, rate=1/2),
         function() rng$exp(chunk, 2))

run_case_rp_only = function(name, f_rp) {
  case_seed = draw_seed()
  restore_rng(case_seed)
  rp_ns = time_dist(f_rp, chunk, reps, bench_time)
  dqrng_str = sprintf("%10s", "")
  factor_d_str = sprintf("%8s", "")
  cat(sprintf(rp_only_fmt,
              name, "-", rp_ns, dqrng_str, "-", factor_d_str))
}

run_case("lognormal(0,1)",
         function() rlnorm(chunk, meanlog=0, sdlog=1),
         function() rng$lognormal(chunk, 0, 1))

run_case_rp_only("skew_norm(0,1,2)",
                 function() rng$skew_normal(chunk, 0, 1, 2))

run_case_rp_only("gumbel(0,1)",
                 function() rng$gumbel(chunk, 0, 1))

run_case_rp_only("pareto(1,2)",
                 function() rng$pareto(chunk, 1, 2))

run_case("gamma(2,3)",
         function() rgamma(chunk, shape=2, scale=3),
         function() rng$gamma(chunk, 2, 3))

run_case("beta(2,5)",
         function() rbeta(chunk, shape1=2, shape2=5),
         function() rng$beta(chunk, 2, 5))

run_case("chi2(5)",
         function() rchisq(chunk, df=5),
         function() rng$chi2(chunk, 5))

run_case("t(10)",
         function() rt(chunk, df=10),
         function() rng$t(chunk, 10))

run_case("F(5,10)",
         function() rf(chunk, df1=5, df2=10),
         function() rng$f(chunk, 5, 10))

run_case("weibull(2,3)",
         function() rweibull(chunk, shape=2, scale=3),
         function() rng$weibull(chunk, 2, 3))

run_case("weibull(3,4)",
         function() rweibull(chunk, shape=3, scale=4),
         function() rng$weibull(chunk, 3, 4))
