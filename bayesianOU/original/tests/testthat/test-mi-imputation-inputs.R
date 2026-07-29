# Unified external multiple imputation (D-IMPL-13.5): the MI driver accepts
# X / V_value / CAPITAL_TOTAL / COM either as a fixed matrix (legacy path,
# bit-identical) or as a per-imputation set (list of D matrices / [T, S, D]
# array) paired with phi_draws. These tests are Stan-free: they exercise the
# input-normalization helpers, the picker, the up-front length validation, and
# the per-imputation wiring of the fit loop (via a mocked engine).

test_that(".as_imp_input classifies NULL / matrix / data.frame as fixed-or-null", {
  expect_null(bayesianOU:::.as_imp_input(NULL, "V_value", 3L))

  m  <- matrix(as.numeric(1:6), 3, 2)
  im <- bayesianOU:::.as_imp_input(m, "X", 3L)
  expect_false(im$varies)
  expect_null(im$list)
  expect_identical(im$fixed, m)

  df  <- as.data.frame(m)                       # a data.frame is a list: must be fixed
  idf <- bayesianOU:::.as_imp_input(df, "X", 3L)
  expect_false(idf$varies)
  expect_equal(idf$fixed, as.matrix(df), ignore_attr = TRUE)
})

test_that(".as_imp_input accepts a list of D matrices and a [T, S, D] array", {
  T0 <- 4L; S0 <- 2L; D <- 3L
  lst <- lapply(seq_len(D), function(k) matrix(as.numeric(k), T0, S0))
  il  <- bayesianOU:::.as_imp_input(lst, "X", D)
  expect_true(il$varies)
  expect_length(il$list, D)
  expect_equal(il$list[[2L]], matrix(2, T0, S0))

  arr <- array(0, dim = c(T0, S0, D))
  for (k in seq_len(D)) arr[, , k] <- k
  ia  <- bayesianOU:::.as_imp_input(arr, "X", D)
  expect_true(ia$varies)
  expect_length(ia$list, D)
  expect_equal(ia$list[[3L]], matrix(3, T0, S0))
})

test_that(".as_imp_input rejects a per-imputation set whose length != D", {
  T0 <- 4L; S0 <- 2L
  lst_bad <- lapply(1:2, function(k) matrix(as.numeric(k), T0, S0))   # 2 != D = 3
  expect_error(bayesianOU:::.as_imp_input(lst_bad, "CAPITAL_TOTAL", 3L),
               "list of 2 matrices")
  arr_bad <- array(0, dim = c(T0, S0, 2L))                            # slice dim 2 != 3
  expect_error(bayesianOU:::.as_imp_input(arr_bad, "COM", 3L),
               "\\[T, S, 2\\] array")
})

test_that(".pick_imp returns NULL, the fixed matrix, or the m-th slice", {
  expect_null(bayesianOU:::.pick_imp(NULL, 1L))

  m   <- matrix(as.numeric(1:6), 3, 2)
  fix <- bayesianOU:::.as_imp_input(m, "X", 3L)
  expect_identical(bayesianOU:::.pick_imp(fix, 1L), m)
  expect_identical(bayesianOU:::.pick_imp(fix, 3L), m)        # same object for every draw

  lst <- lapply(1:3, function(k) matrix(as.numeric(k), 3, 2))
  var <- bayesianOU:::.as_imp_input(lst, "X", 3L)
  expect_equal(bayesianOU:::.pick_imp(var, 2L), matrix(2, 3, 2))
})

test_that("fit_ou_nested_mi rejects a per-imputation input of the wrong length before any fit", {
  D <- 3L; T0 <- 8L; S0 <- 2L
  ph <- replicate(D, matrix(stats::rnorm(T0 * S0), T0, S0), simplify = FALSE)
  X_bad <- replicate(2L, matrix(stats::rnorm(T0 * S0), T0, S0), simplify = FALSE)  # 2 != 3
  expect_error(
    fit_ou_nested_mi(phi_draws = ph, X = X_bad, TMG = stats::rnorm(T0),
                     COM = matrix(1, T0, S0), CAPITAL_TOTAL = matrix(1, T0, S0),
                     Gprime = stats::rnorm(T0), n_levels = 2L, M = D),
    "list of 2 matrices")
})

test_that("the fit loop pairs each per-imputation input with its draw (D-IMPL-13.5 wiring)", {
  skip_if_not(exists("local_mocked_bindings", where = asNamespace("testthat"),
                     inherits = FALSE),
              "testthat 3e local_mocked_bindings required")

  D <- 3L; T0 <- 8L; S0 <- 2L
  ph <- lapply(seq_len(D), function(k) matrix(as.numeric(10 * k), T0, S0))
  X  <- lapply(seq_len(D), function(k) matrix(as.numeric(100 * k), T0, S0))
  K  <- lapply(seq_len(D), function(k) matrix(as.numeric(200 * k), T0, S0))
  C  <- lapply(seq_len(D), function(k) matrix(as.numeric(300 * k), T0, S0))
  V  <- lapply(seq_len(D), function(k) matrix(as.numeric(400 * k), T0, S0))
  Gp <- stats::rnorm(T0); Tm <- stats::rnorm(T0)

  captured <- NULL
  fake_fit <- function(Y, X, TMG, COM, CAPITAL_TOTAL, Gprime = NULL,
                       V_value = NULL, n_levels = 1L, seed = 1234,
                       verbose = FALSE, ...) {
    captured <<- list(Y = Y, X = X, COM = COM, CAPITAL_TOTAL = CAPITAL_TOTAL,
                      V_value = V_value, TMG = TMG, Gprime = Gprime, seed = seed)
    stop("STOP_SENTINEL")
  }

  testthat::local_mocked_bindings(fit_ou_nested = fake_fit)
  # First imputation uses draw m = 1: every per-imputation input must be its slice 1,
  # the fixed aggregates pass through, and the seed is offset by the draw index.
  expect_error(
    fit_ou_nested_mi(phi_draws = ph, X = X, TMG = Tm, COM = C, CAPITAL_TOTAL = K,
                     Gprime = Gp, V_value = V, n_levels = 3L, M = D, seed = 1000),
    "STOP_SENTINEL")
  expect_equal(captured$Y,             ph[[1L]])
  expect_equal(captured$X,             X[[1L]])
  expect_equal(captured$COM,           C[[1L]])
  expect_equal(captured$CAPITAL_TOTAL, K[[1L]])
  expect_equal(captured$V_value,       V[[1L]])
  expect_identical(captured$TMG,    Tm)          # aggregates fixed across imputations
  expect_identical(captured$Gprime, Gp)
  expect_identical(captured$seed,   1001)        # seed + m, m = 1
})

test_that("the fixed-matrix path passes the same anchor verbatim (legacy, bit-identical)", {
  skip_if_not(exists("local_mocked_bindings", where = asNamespace("testthat"),
                     inherits = FALSE),
              "testthat 3e local_mocked_bindings required")

  D <- 3L; T0 <- 8L; S0 <- 2L
  ph <- lapply(seq_len(D), function(k) matrix(as.numeric(10 * k), T0, S0))
  Xf <- matrix(stats::rnorm(T0 * S0), T0, S0)
  Kf <- matrix(1, T0, S0); Cf <- matrix(1, T0, S0)

  captured <- NULL
  fake_fit <- function(Y, X, TMG, COM, CAPITAL_TOTAL, Gprime = NULL,
                       V_value = NULL, n_levels = 1L, seed = 1234,
                       verbose = FALSE, ...) {
    captured <<- list(X = X, V_value = V_value); stop("STOP_SENTINEL")
  }
  testthat::local_mocked_bindings(fit_ou_nested = fake_fit)
  expect_error(
    fit_ou_nested_mi(phi_draws = ph, X = Xf, TMG = stats::rnorm(T0), COM = Cf,
                     CAPITAL_TOTAL = Kf, Gprime = stats::rnorm(T0),
                     n_levels = 2L, M = D),
    "STOP_SENTINEL")
  expect_identical(captured$X, Xf)               # fixed matrix forwarded unchanged
  expect_null(captured$V_value)                  # NULL below Level 3
})
