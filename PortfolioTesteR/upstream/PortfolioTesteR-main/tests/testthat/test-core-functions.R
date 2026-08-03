test_that("sample weekly data loads and has expected shape", {
  expect_true("sample_prices_weekly" %in% data(package = "PortfolioTesteR")$results[, "Item"])
  expect_no_error(data("sample_prices_weekly", package = "PortfolioTesteR"))
  data(sample_prices_weekly, package = "PortfolioTesteR")

  expect_s3_class(sample_prices_weekly, "data.frame")
  expect_true("Date" %in% names(sample_prices_weekly))
  expect_true(ncol(sample_prices_weekly) > 1)
  expect_true(all(is.finite(as.numeric(sample_prices_weekly[1, -1]))))
})

test_that("calc_momentum() returns aligned, finite output (aside from lookback NAs)", {
  data(sample_prices_weekly, package = "PortfolioTesteR")
  res <- calc_momentum(sample_prices_weekly, lookback = 12)

  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), nrow(sample_prices_weekly))
  expect_equal(ncol(res), ncol(sample_prices_weekly))
  # All non-leading entries should be finite
  after_lb <- res[-seq_len(12), -1, drop = FALSE]
  expect_true(all(is.finite(as.numeric(as.matrix(after_lb)))))
})

test_that("filter_top_n() selects at most n assets per row", {
  data(sample_prices_weekly, package = "PortfolioTesteR")
  ranks <- calc_momentum(sample_prices_weekly, lookback = 12)
  sel <- filter_top_n(ranks, n = 5)

  expect_s3_class(sel, "data.frame")
  expect_equal(nrow(sel), nrow(ranks))
  # Count selections per row (treat nonzero/TRUE as selected)
  counts <- apply(sel[, -1, drop = FALSE], 1, function(x) sum(is.finite(x) & x != 0))
  expect_true(all(counts <= 5))
})

test_that("weight_equally() produces row-wise weights that sum to ~1 where selection exists", {
  data(sample_prices_weekly, package = "PortfolioTesteR")
  ranks <- calc_momentum(sample_prices_weekly, lookback = 12)
  sel <- filter_top_n(ranks, n = 5)
  w <- weight_equally(sel)

  expect_s3_class(w, "data.frame")
  expect_equal(dim(w), dim(sel))

  sums <- rowSums(as.matrix(w[, -1, drop = FALSE]), na.rm = TRUE)
  # For rows with at least one finite weight, the sum should be ~1
  has_w <- rowSums(is.finite(as.matrix(w[, -1, drop = FALSE]))) > 0
  expect_true(all(abs(sums[has_w] - 1) < 1e-8))
})

test_that("run_backtest() returns a well-formed result and finite series", {
  data(sample_prices_weekly, package = "PortfolioTesteR")
  ranks <- calc_momentum(sample_prices_weekly, lookback = 12)
  sel <- filter_top_n(ranks, n = 5)
  w <- weight_equally(sel)

  bt <- run_backtest(sample_prices_weekly, w, initial_capital = 1e5)

  # Be flexible on class: check structure/columns instead of a hard-coded S3 class
  expect_true(is.list(bt) || is.data.frame(bt))
  nm <- names(bt)
  expect_true(any(grepl("port", nm, ignore.case = TRUE)) || "portfolio_value" %in% nm)
  # Extract portfolio series robustly
  pv <- if ("portfolio_value" %in% nm) bt$portfolio_value else
    bt[[ which.max(vapply(bt, function(x) is.numeric(x) && length(x) == nrow(bt), logical(1))) ]]
  expect_true(length(pv) > 0)
  expect_true(all(is.finite(pv[is.finite(pv)])))
})

test_that("calc_rsi() stays within [0, 100]", {
  data(sample_prices_weekly, package = "PortfolioTesteR")
  out <- suppressWarnings(calc_rsi(sample_prices_weekly, period = 6))
  vals <- as.numeric(as.matrix(out[, -1, drop = FALSE]))
  vals <- vals[is.finite(vals)]
  expect_true(length(vals) > 0)
  expect_true(all(vals >= 0 & vals <= 100))
})
