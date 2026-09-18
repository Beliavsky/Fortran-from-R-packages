# inst/examples/DrawDemo.R
# Purpose: deterministic demo that prints a few unif draws for selected engines.
# Matching C program: examples/DrawDemo.c

suppressMessages(library(randompack))

fmt1 <- "%s seed=%d spawnkey=%d unif=%.17g %.17g\n"
fmt2 <- "%s counter=%d key=%d unif=%.17g %.17g\n"

u64_halves <- function(x) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x))
    stop("x must be a length-1 numeric")
  x <- as.numeric(x)
  if (x < 0 || x > 18446744073709551615)
    stop("x must be in [0, 2^64-1]")
  lo <- x %% 4294967296
  hi <- floor(x / 4294967296)
  c(lo, hi)
}

draw_u01 <- function(engine, seed, key) {
  rng <- randompack_rng(engine)
  if (is.null(rng)) stop("engine creation failed")
  rng$seed(seed, as.integer(key))
  x <- rng$unif(2)
  cat(sprintf(fmt1, engine, seed, key, x[1], x[2]))
  invisible(NULL)
}

draw_philox <- function(counter0, key0) {
  rng <- randompack_rng("philox")
  if (is.null(rng)) stop("engine creation failed")
  state <- c(u64_halves(counter0), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
  key <- c(u64_halves(key0), 0, 0)
  rng$set_state(state)
  rng$philox_set_key(key)
  x <- rng$unif(2)
  cat(sprintf(fmt2, "philox", counter0, key0, x[1], x[2]))
  invisible(NULL)
}

draw_squares <- function(counter, key) {
  rng <- randompack_rng("squares")
  state <- c(u64_halves(counter), 0, 0)
  key_vec <- u64_halves(key)
  rng$set_state(state)
  rng$squares_set_key(key_vec)
  x <- rng$unif(2)
  cat(sprintf(fmt2, "squares", counter, key, x[1], x[2]))
  invisible(NULL)
}

main <- function() {
  draw_u01("x256++", 42, 2)
  draw_u01("pcg64", 42, 2)
  draw_philox(1, 42)
  draw_u01("chacha20", 42, 2)
  draw_squares(1, 123456789)
  invisible(NULL)
}

main()
