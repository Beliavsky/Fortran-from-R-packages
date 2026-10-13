# Reference assertions for test_draws.f90, evaluated against posterior 1.7.0.
# Run: Rscript test/reference_draws.R
suppressPackageStartupMessages(library(posterior))
stopifnot(as.character(packageVersion("posterior")) == "1.7.0")
checks <- 0L
same <- function(actual, expected, tolerance = 2e-12) {
  stopifnot(isTRUE(all.equal(as.numeric(actual), as.numeric(expected),
                           tolerance = tolerance, check.attributes = FALSE)))
  checks <<- checks + 1L
}
raw <- function(x) as.vector(x)
x <- as_draws_array(array(1:20, c(5, 2, 2)))
small <- draws_array(a = 1:4, b = 5:8)
same(as_draws_matrix(x), 1:20)
same(as_draws_array(as_draws_matrix(x)), x)
same(merge_chains(x), 1:20)
same(suppressWarnings(split_chains(x)),
     c(1, 2, 6, 7, 3, 4, 8, 9, 11, 12, 16, 17, 13, 14, 18, 19))
stopifnot(identical(dim(suppressWarnings(split_chains(x[1, , ]))), c(0L, 4L, 2L)))
same(bind_draws(x, x, along = "iteration"),
     array(rep(c(1:5, 1:5, 6:10, 6:10, 11:15, 11:15, 16:20, 16:20)), c(10, 2, 2)))
stopifnot(identical(dim(bind_draws(x, x, along = "chain")), c(5L, 4L, 2L)))
y <- x
variables(y) <- c("other1", "other2")
stopifnot(identical(dim(bind_draws(x, y, along = "variable")), c(5L, 2L, 4L)))
same(subset_draws(x, iteration = c(4, 2, 2), chain = c(2, 1), unique = FALSE),
     c(2, 2, 4, 7, 7, 9, 12, 12, 14, 17, 17, 19))
same(subset_draws(x, iteration = c(4, 2, 2), variable = rev(variables(x)))[, 1, 1], c(12, 14))
same(subset_draws(x, iteration = c(2, 5), exclude = TRUE)[, 1, 1], c(1, 3, 4))
same(subset_draws(x, draw = c(9, 2, 2), unique = FALSE)[, 1, 1], c(2, 2, 9))
same(thin_draws(x, 1.5)[, 1, 1], c(1, 2, 4))
same(thin_draws(x, 2.5)[, 1, 1], c(1, 4))
y <- x
dimnames(y)[[1]] <- as.character(5:1)
dimnames(y)[[2]] <- as.character(2:1)
same(order_draws(y)[, 1, 1], 10:6)
same(resample_draws(small, weights = 1:4, ndraws = 7, method = "deterministic")[, 1, 1],
     c(1, 2, 3, 3, 4, 4, 4))
set.seed(867)
same(runif(4), c(0.703724499791861, 0.708184336544946, 0.749903351534158, 0.468353937380016))
set.seed(867)
same(resample_draws(small, weights = 1:4, ndraws = 7, method = "stratified")[, 1, 1],
     c(2, 2, 3, 3, 4, 4, 4))
same(extract_variable(x, variables(x)[2]), 11:20)
same(extract_variable_matrix(x, variables(x)[2]), matrix(11:20, 5))
w <- weight_draws(small, c(0, 1, 2, 3))
same(weights(w), c(0, 1, 2, 3) / 6)
same(weights(w, normalize = FALSE), c(0, 1, 2, 3))
w <- weight_draws(w, c(1000, 1000, -Inf, 1000), log = TRUE)
same(weights(w), c(1, 1, 0, 1) / 3)
stopifnot(all(is.nan(weights(weight_draws(small, rep(0, 4))))))
rows <- rbind(c(1, 3, 5, 7), c(2, 4, 6, 8), c(-1, 0, 1, 2))
r <- rvar(rows)
expected <- list(mean = c(4, 5, .5), median = c(4, 5, .5), sum = c(16, 20, 2),
                 prod = c(105, 384, 0), min = c(1, 2, -1), max = c(7, 8, 2),
                 var = c(20, 20, 5) / 3, sd = sqrt(c(20, 20, 5) / 3),
                 mad = c(2.9652, 2.9652, 1.4826))
for (method in names(expected)) {
  f <- get(paste0("rvar_", method))
  same(draws_of(f(r)), expected[[method]])
}
same(draws_of(rvar_mad(r, constant = 1)), c(2, 2, 1))
same(draws_of(rvar_quantile(r, probs = c(0, .25, .5, 1))),
     c(1, 2, -1, 2.5, 3.5, -.25, 4, 5, .5, 7, 8, 2))
same(draws_of(rvar_range(r)), c(1, 2, -1, 7, 8, 2))
same(E(r), c(2, 7, 12, 17) / 3)
same(var(rvar(t(rows))), c(20, 20, 5) / 3)
same(variance(rvar(t(rows))), c(20, 20, 5) / 3)
same(sd(rvar(t(rows))), sqrt(c(20, 20, 5) / 3))
same(mad(rvar(t(rows))), c(2.9652, 2.9652, 1.4826))
special <- rvar(rbind(c(1, NaN, 3), c(Inf, Inf, Inf), c(-Inf, Inf, 3)))
same(draws_of(rvar_mean(special)), c(NaN, Inf, NaN))
same(draws_of(rvar_median(special)), c(NaN, Inf, 3))
same(draws_of(rvar_mad(special)), c(NaN, NaN, Inf))
stopifnot(all(is.na(draws_of(rvar_var(special)))))
q <- draws_of(rvar_quantile(special, probs = c(0, .25, .5, .75, 1)))
stopifnot(all(is.na(q[1, ])))
same(q[2, ], rep(Inf, 5))
same(q[3, ], c(-Inf, -Inf, 3, Inf, Inf))
same(draws_of(rvar_quantile(special, probs = c(0, .5, 1), na.rm = TRUE))[1, ], c(1, 2, 3))
empty <- rvar(matrix(numeric(), 1, 0))
same(draws_of(rvar_sum(empty)), 0)
same(draws_of(rvar_prod(empty)), 1)
same(suppressWarnings(draws_of(rvar_min(empty))), Inf)
same(suppressWarnings(draws_of(rvar_max(empty))), -Inf)
stopifnot(is.na(draws_of(rvar_mean(empty))))
logical_rows <- rbind(c(FALSE, NA, TRUE), c(FALSE, NA, FALSE), c(TRUE, NA, TRUE))
l <- rvar(logical_rows)
same(draws_of(rvar_all(l)), c(FALSE, FALSE, NA))
same(draws_of(rvar_any(l)), c(TRUE, NA, TRUE))
same(draws_of(rvar_all(l, na.rm = TRUE)), c(FALSE, FALSE, TRUE))
same(draws_of(rvar_any(l, na.rm = TRUE)), c(TRUE, FALSE, TRUE))
same(draws_of(rvar_any(r)), c(TRUE, TRUE, TRUE))
same(Pr(l, na.rm = TRUE), c(1 / 3, NaN, 2 / 3))
same(cdf(rvar(c(1, 2, NaN, Inf)), c(-Inf, 1, 1.5, 2, Inf, NaN)), c(0, 1 / 3, 1 / 3, 2 / 3, 1, NaN))
stopifnot(sum(draws_of(rvar_is_finite(special))) == 3,
          sum(draws_of(rvar_is_infinite(special))) == 5,
          sum(draws_of(rvar_is_nan(special))) == 1,
          sum(draws_of(rvar_is_na(special))) == 1)
left <- array(0, c(2, 2, 2))
left[1, , ] <- matrix(1:4, 2)
left[2, , ] <- 2 * left[1, , ]
product <- draws_of(rvar(left) %**% diag(c(1, 2)))
same(product[1, , ], c(1, 2, 6, 8))
same(product[2, , ], c(2, 4, 12, 16))
same(draws_of(as_rvar(c(10, 20, 30))), c(10, 20, 30))
y <- draws_of(rvar_ifelse(l, c(10, 20, 30), -1))
same(y[1, ], c(-1, NA, 30))
same(y[2, ], c(-1, NA, -1))
same(y[3, ], c(10, NA, 30))
cat("posterior 1.7.0 draw/rvar reference: ", checks, " numerical comparisons passed\n", sep = "")
