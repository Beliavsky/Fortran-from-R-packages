# tests/testthat/test-ml-helpers.R

test_that("feature prep lags correctly (no look-ahead)", {
  data(sample_prices_weekly); data(sample_prices_daily)

  # Expected: manual lag of 12w momentum
  raw_lag <- panel_lag(calc_momentum(sample_prices_weekly, 12L), 1L)

  # Actual: helper output
  X <- ml_prepare_features(sample_prices_weekly, sample_prices_daily, include = c("mom12"))

  # Compare per-symbol series (Date-aligned)
  syms_exp <- setdiff(names(raw_lag), "Date")
  syms_act <- setdiff(names(X$mom12), "Date")
  syms     <- intersect(syms_exp, syms_act)
  expect_gt(length(syms), 0L)

  for (s in syms) {
    expect_equal(X$mom12[[s]], raw_lag[[s]], tolerance = 1e-12)
  }
})

test_that("panel ops preserve Date and symbol sets", {
  data(sample_prices_weekly)

  A <- panel_lag(calc_momentum(sample_prices_weekly, 12L), 1L)
  B <- panel_lag(calc_momentum(sample_prices_weekly, 26L), 1L)

  C_int <- ml_panel_op(A, B, op = `+`, how = "intersect")
  expect_true(inherits(C_int$Date, "Date"))
  expect_true(all(intersect(setdiff(names(A),"Date"), setdiff(names(B),"Date")) %in%
                    setdiff(names(C_int),"Date")))

  C_uni <- ml_panel_op(A, B, op = `+`, how = "union")
  expect_true(all(union(setdiff(names(A),"Date"), setdiff(names(B),"Date")) %in%
                    setdiff(names(C_uni),"Date")))
})

test_that("ml_add_interactions builds product panels", {
  data(sample_prices_weekly); data(sample_prices_daily)

  X0 <- ml_prepare_features(
    sample_prices_weekly, sample_prices_daily,
    include = c("mom12","rsi14")
  )
  X  <- ml_add_interactions(X0, list(mom_rsi = c("mom12","rsi14")))

  expect_true("mom_rsi" %in% names(X))
  syms <- intersect(setdiff(names(X0$mom12),"Date"), setdiff(names(X0$rsi14),"Date"))
  s <- syms[1]
  exp_vec <- X0$mom12[[s]] * X0$rsi14[[s]]
  act_vec <- X$mom_rsi[[s]]
  expect_equal(act_vec, exp_vec, tolerance = 1e-12)
})

test_that("ensemble is NA-safe and shape-safe", {
  set.seed(1)
  X <- as.data.frame(matrix(rnorm(200), 100, 2))
  y <- rnorm(100)
  ridge <- ml_make_model("ridge")
  rf    <- ml_make_model("rf", num.trees = 50)
  ens   <- ml_make_ensemble(ridge, rf)

  fit  <- ens$fit(X, y)
  pred <- ens$predict(fit, X)
  expect_length(pred, nrow(X))

  X_na <- X; X_na[1,1] <- NA_real_
  pred2 <- ens$predict(fit, X_na)
  expect_length(pred2, nrow(X_na))
})

test_that("sequence factory works in linear mode (no TF)", {
  m <- ml_make_seq_model("linear", steps = 4L)  # stays entirely in base R
  X <- matrix(rnorm(10 * 12), nrow = 10, ncol = 12)
  y <- rnorm(10)
  fit <- m$fit(X, y)
  p   <- m$predict(fit, X)
  expect_length(p, nrow(X))
  expect_true(is.numeric(p))
})


test_that("backtest invariants hold on active dates (smoke, linear)", {
  data(sample_prices_weekly); data(sample_prices_daily)

  X <- ml_prepare_features(sample_prices_weekly, sample_prices_daily)
  Y <- make_labels(sample_prices_weekly, horizon = 4L, type = "log")
  lin <- ml_make_model("linear")

  set.seed(1)
  res <- ml_backtest(
    features_list = X, labels = Y,
    fit_fn = lin$fit, predict_fn = lin$predict,
    schedule = list(is = 104L, oos = 4L, step = 4L),
    transform = "zscore",
    selection = list(top_k = 10L),
    weighting = list(method = "rank"),
    prices = sample_prices_weekly,
    initial_capital = 1e5
  )

  wt <- as.data.frame(res$weights)
  stopifnot("Date" %in% names(wt))

  # Detect long vs wide weights
  is_long <- "Symbol" %in% names(wt)

  # Build per-date weight sums on ACTIVE (non-empty) rows only
  if (is_long) {
    # Find numeric weight column robustly
    val_col <- if ("value" %in% names(wt)) {
      "value"
    } else if ("weight" %in% names(wt)) {
      "weight"
    } else {
      setdiff(names(wt)[vapply(wt, is.numeric, TRUE)], c("Date","Symbol"))[1]
    }
    # Keep rows with finite weights
    wt_active <- wt[is.finite(wt[[val_col]]), , drop = FALSE]
    sums <- tapply(wt_active[[val_col]], wt_active$Date, function(z) sum(z, na.rm = TRUE))
  } else {
    # Wide: sum across all numeric columns except Date
    num_cols <- setdiff(names(wt)[vapply(wt, is.numeric, TRUE)], "Date")
    row_sums <- rowSums(wt[, num_cols, drop = FALSE], na.rm = TRUE)
    # Active rows: at least one finite (non-zero) weight present
    active_idx <- is.finite(row_sums) & (row_sums > 0)
    sums <- row_sums[active_idx]
  }

  expect_true(length(sums) > 0L)
  expect_true(max(abs(sums - 1), na.rm = TRUE) < 1e-6)

  # No NA weights on active rows
  if (is_long) {
    val_col <- if ("value" %in% names(wt)) "value" else if ("weight" %in% names(wt)) "weight" else
      setdiff(names(wt)[vapply(wt, is.numeric, TRUE)], c("Date","Symbol"))[1]
    active_idx <- is.finite(wt[[val_col]])
    expect_false(any(is.na(wt[[val_col]][active_idx])))
  } else {
    num_cols <- setdiff(names(wt)[vapply(wt, is.numeric, TRUE)], "Date")
    active_idx <- rowSums(is.finite(as.matrix(wt[, num_cols, drop = FALSE]))) > 0
    expect_false(any(is.na(as.matrix(wt[active_idx, num_cols, drop = FALSE]))))
  }
})

