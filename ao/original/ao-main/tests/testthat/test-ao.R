check_ao_output <- function(x, target = NULL, add_details = TRUE) {
  has_target <- !is.null(target)
  names <- c("estimate", "value", "seconds", "stopping_reason")
  if (add_details) {
    names <- oeli::insert_vector_entry(names, "details", 2)
  }
  if (has_target) {
    names <- oeli::insert_vector_entry(names, "estimate_split", 1)
  }
  checkmate::expect_list(x, len = length(names))
  checkmate::expect_names(names(x), identical.to = names)
  if (has_target) {
    checkmate::expect_list(x$estimate_split, len = 3)
    checkmate::expect_names(names(x$estimate_split), identical.to = target)
  }
}

test_that("ao works", {
  himmelblau <- function(x) (x[1]^2 + x[2] - 11)^2 + (x[1] + x[2]^2 - 7)^2
  ao_out <- ao(f = himmelblau, initial = c(0, 0))
  check_ao_output(ao_out)
  ao_out <- ao(f = himmelblau, initial = c(0, 0), add_details = FALSE)
  check_ao_output(ao_out, add_details = FALSE)
})

test_that("ao with custom gradient works", {
  himmelblau <- function(x) (x[1]^2 + x[2] - 11)^2 + (x[1] + x[2]^2 - 7)^2
  gradient <- function(x) {
    c(
      4 * x[1] * (x[1]^2 + x[2] - 11) + 2 * (x[1] + x[2]^2 - 7),
      2 * (x[1]^2 + x[2] - 11) + 4 * x[2] * (x[1] + x[2]^2 - 7)
    )
  }
  ao_out <- ao(f = himmelblau, initial = c(0, 0), gradient = gradient)
  check_ao_output(ao_out)
  ao_random_out <- ao(
    f = himmelblau, initial = c(0, 0), partition = "random", gradient = gradient
  )
  check_ao_output(ao_random_out)
})

test_that("ao with custom Hessian works", {
  f <- function(x) {
    # print("huhu, my name is f")
    x^2
  }
  g <- function(x) {
    # print("huhu, my name is g")
    2*x
  }
  h <- function(x) {
    # print("huhu, my name is h")
    matrix(2) # important to have matrix(2) instead of 2!
  }
  ao_out <- ao(
    f = f, initial = 10, gradient = g, hessian = h,
    # switch to nlm because stats::optim does not support Hessian
    base_optimizer = optimizeR::Optimizer$new("stats::nlm")
  )
  check_ao_output(ao_out)
})

test_that("ao with additional parameters works", {
  himmelblau <- function(x, a, b, c) {
    (x[1]^a + x[2] - b)^a + (x[1] + x[2]^a - c)^a
  }
  expect_error(
    ao(f = himmelblau, initial = c(0, 0)),
    "Function argument `a` is required but not specified yet."
  )
  ao_out <- ao(f = himmelblau, initial = c(0, 0), a = 2, b = 11, c = 7)
  check_ao_output(ao_out)
})

test_that("ao with targets works", {
  himmelblau <- function(x, a, b, c) {
    (x[1]^a + x[2] - b)^a + (x[1] + x[2]^a - c)^a
  }
  target <- c("x", "b", "c")
  ao_out <- ao(
    f = himmelblau, initial = 1:4, target = target, npar = c(2, 1, 1), a = 1
  )
  check_ao_output(ao_out, target = target)
})

test_that("ao with NULL values for fixed arguments works", {
  f <- function(x, a, b, ind) {
    if (is.null(ind)) {
      (x[1]^2 + x[2] + a)^2 + (x[1] + x[2]^2 + b)^2 + (x[3] - 1)^2
    }
  }
  ao_out <- ao(f = f, initial = c(0, 0, 0), a = -11, b = -7, ind = NULL)
  check_ao_output(ao_out)
})

test_that("ao with a different base optimizer works", {
  himmelblau <- function(x) (x[1]^2 + x[2] - 11)^2 + (x[1] + x[2]^2 - 7)^2
  ao_out <- ao(
    f = himmelblau, initial = c(0, 0),
    base_optimizer = optimizeR::Optimizer$new("stats::nlm")
  )
  check_ao_output(ao_out)
})

test_that("ao with custom partition works", {
  f <- function(x) (x[1]^2 + x[2])^2 + (x[1] + x[2]^2)^2 + (x[3] + x[4]^2)^2
  ao_custom <- ao(
    f = f, initial = c(1, 1, 1, 1), partition = list(1, 2, 3:4)
  )
  check_ao_output(ao_custom)
})

test_that("ao with parameter bounds works", {
  rosenbrock <- function(x) (1 - x[1])^2 + (x[2] - x[1]^2)^2
  lower <- 1.5
  upper <- 2.5
  initial <- runif(2, min = lower, max = upper)
  ao_bounds <- ao(
    f = rosenbrock, initial = initial, lower = lower, upper = upper
  )
  check_ao_output(ao_bounds)
  expect_true(all(ao_bounds$estimate <= upper))
  expect_true(all(ao_bounds$estimate >= lower))
})

test_that("ao rejects inconsistent target and npar inputs", {
  f <- function(x, y) sum(x^2) + sum(y^2)
  expect_error(ao(f = f, initial = 1:2, target = c("x", "y")))
})

test_that("ao rejects invalid custom partitions", {
  f <- function(x) sum(x^2)
  expect_error(ao(f = f, initial = 1:2, partition = list(1, 3)))
  expect_error(ao(f = f, initial = 1:2, partition = list(1, NA_integer_)))
})

test_that("ao rejects starting values outside bounds", {
  f <- function(x) sum(x^2)
  expect_error(ao(f = f, initial = 0, lower = 1))
  expect_error(ao(f = f, initial = 2, upper = 1))
})

test_that("multiple ao processes work", {
  himmelblau <- function(x) (x[1]^2 + x[2] - 11)^2 + (x[1] + x[2]^2 - 7)^2
  out_multi <- ao(
    f = himmelblau,
    initial = list(c(0, 0), c(1, 1)),
    partition = list("random", list(2, 1), "none", "sequential"),
    base_optimizer = list(
      optimizeR::Optimizer$new("stats::nlm"),
      optimizeR::Optimizer$new("stats::optim")
    ),
    add_details = FALSE
  )
  check_ao_output(out_multi, add_details = FALSE)
})

test_that("multiple ao processes can omit process details", {
  f <- function(x) sum((x - 1)^2)
  out_multi <- ao(
    f = f,
    initial = list(c(0, 0), c(2, 2)),
    partition = "none",
    iteration_limit = 1,
    add_details = FALSE
  )
  check_ao_output(out_multi, add_details = FALSE)
  expect_false(
    any(
      c("details", "estimates", "values", "seconds_each") %in%
        names(out_multi)
    )
  )
})

test_that("ao fails graciously", {
  f <- function(x) x[1]^2 + x[2]^2
  expect_silent(ao(f = f, initial = 1)) # initial is mis-specified
})

test_that("random partition works with ucminf optimizer", {
  skip_if_not_installed("ucminf")
  set.seed(1)

  # simple convex objective in 6D
  f <- function(x) sum((x - 1)^2)
  out <- ao(
    f = f,
    initial = rep(0, 6),
    partition = "random",
    new_block_probability = 0.35,
    minimum_block_number = 2,
    base_optimizer = optimizeR::Optimizer$new("ucminf::ucminf"),
    iteration_limit = 3, # keep it fast for CRAN
    add_details = TRUE
  )
  check_ao_output(out, add_details = TRUE)
  d <- out$details
  bcols <- grep("^b[0-9]+$", names(d), value = TRUE)
  expect_length(bcols, 6)

  # drop init row if present
  d2 <- d[d$iteration > 0, , drop = FALSE]

  # each step should activate at least one parameter and not exceed npar
  active_counts <- rowSums(as.matrix(d2[, bcols, drop = FALSE]) == 1)
  expect_all_true(active_counts >= 1)
  expect_all_true(active_counts <= 6)

  # block indicators must be binary (0/1)
  expect_all_true(unlist(d2[, bcols, drop = FALSE]) %in% c(0, 1))
})
