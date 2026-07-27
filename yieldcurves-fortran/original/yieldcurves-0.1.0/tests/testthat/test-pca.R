test_that("yc_pca returns valid object", {
  set.seed(42)
  n <- 50
  tenors <- c(1, 2, 5, 10, 30)
  curves <- matrix(rnorm(n * length(tenors), mean = 0.04, sd = 0.005),
                   n, length(tenors))
  colnames(curves) <- paste0(tenors, "Y")

  result <- yc_pca(curves)
  expect_s3_class(result, "yc_pca")
  expect_equal(result$n_components, 3)
  expect_equal(nrow(result$loadings), 5)
  expect_equal(ncol(result$loadings), 3)
  expect_equal(nrow(result$scores), 50)
})

test_that("yc_pca variance explained sums correctly", {
  set.seed(42)
  curves <- matrix(rnorm(100 * 5), 100, 5)
  colnames(curves) <- paste0("T", 1:5)

  result <- yc_pca(curves, n_components = 5)
  expect_equal(sum(result$variance_explained), 1, tolerance = 1e-10)
})

test_that("yc_pca first component is level-like", {
  set.seed(42)
  n <- 200
  tenors <- c(1, 2, 5, 10, 30)

  # Create curves dominated by parallel shifts (level factor)
  level <- cumsum(rnorm(n, 0, 0.002))
  curves <- matrix(0.04, n, length(tenors))
  for (i in seq_len(n)) curves[i, ] <- curves[i, ] + level[i]
  colnames(curves) <- paste0(tenors, "Y")

  result <- yc_pca(curves)
  # First component should explain > 90% of variance
  expect_true(result$variance_explained[1] > 0.9)
  # Level loading should be roughly equal across tenors
  loading_range <- diff(range(result$loadings[, 1]))
  expect_true(loading_range < 0.3)
})

test_that("yc_pca rejects too few rows", {
  expect_error(yc_pca(matrix(1:4, 2, 2)), "at least 3")
})

test_that("yc_pca rejects single column", {
  expect_error(yc_pca(matrix(1:10, 10, 1)), "at least 2")
})

test_that("yc_pca rejects NAs", {
  m <- matrix(1:20, 5, 4)
  m[1, 1] <- NA
  expect_error(yc_pca(m), "NA")
})

test_that("yc_pca reconstruction recovers original data", {
  set.seed(42)
  curves <- matrix(rnorm(50 * 5, mean = 0.04, sd = 0.005), 50, 5)
  colnames(curves) <- paste0("T", 1:5)

  result <- yc_pca(curves, n_components = 5)
  # Reconstruct: scores %*% t(loadings) + center
  reconstructed <- result$scores %*% t(result$loadings) +
    matrix(result$center, nrow = 50, ncol = 5, byrow = TRUE)
  expect_equal(reconstructed, curves, tolerance = 1e-8)
})

test_that("yc_pca print method works", {
  set.seed(42)
  curves <- matrix(rnorm(50 * 5, mean = 0.04), 50, 5)
  colnames(curves) <- paste0("T", 1:5)
  result <- yc_pca(curves)
  expect_no_error(print(result))
})

test_that("plot.yc_pca runs without error", {
  set.seed(42)
  curves <- matrix(rnorm(50 * 5, mean = 0.04), 50, 5)
  colnames(curves) <- paste0("T", 1:5)
  result <- yc_pca(curves)
  expect_no_error(plot(result))
})

test_that("summary.yc_pca runs without error", {
  set.seed(42)
  curves <- matrix(rnorm(50 * 5, mean = 0.04), 50, 5)
  colnames(curves) <- paste0("T", 1:5)
  result <- yc_pca(curves)
  expect_no_error(summary(result))
})
