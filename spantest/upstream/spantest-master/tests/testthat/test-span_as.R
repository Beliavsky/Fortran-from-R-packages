test_that("span_as returns aggregated p-values", {
  set.seed(1)
  bench <- matrix(rnorm(300), nrow = 100, ncol = 3)
  test  <- matrix(rnorm(200), nrow = 100, ncol = 2)

  result <- span_as(bench, test, control = list(ks = c(1/3), L = c(0)))

  expect_type(result, "list")
  expect_true(all(grepl("^CCT[ad]{1,2}_L\\d+_k\\d+$", names(result))))
  expect_true(all(vapply(result, function(x) is.numeric(x) && x >= 0 && x <= 1, logical(1))))
})

test_that("span_as is deterministic and names follow the L / k grid", {
  set.seed(7)
  bench <- matrix(rnorm(450), 150, 3)
  test  <- matrix(rnorm(600), 150, 4)

  r1 <- span_as(bench, test, control = list(ks = c(1/3), L = c(0, 2)))
  r2 <- span_as(bench, test, control = list(ks = c(1/3), L = c(0, 2)))

  expect_identical(r1, r2)  # fixed internal seed -> reproducible, no RNG side effect
  expect_setequal(
    names(r1),
    c("CCTd_L0_k1", "CCTd_L2_k1", "CCTad_L0_k1", "CCTad_L2_k1",
      "CCTa_L0_k1", "CCTa_L2_k1")
  )
})

test_that("span_as detects an alpha alternative", {
  set.seed(3)
  T <- 150; K <- 3; N <- 4
  bench <- matrix(rnorm(T * K), T, K)
  B <- matrix(runif(K * N, 0.2, 1), K, N)
  B <- sweep(B, 2, colSums(B), "/")           # columns sum to 1 -> delta = 0
  test <- 0.5 + bench %*% B + matrix(rnorm(T * N, sd = 0.3), T, N)  # alpha = 0.5

  result <- span_as(bench, test, control = list(ks = c(1/3), L = c(0)))
  expect_lt(result$CCTa_L0_k1, 0.05)   # alpha test should reject
})
