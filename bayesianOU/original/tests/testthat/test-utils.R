test_that("zscore_train standardizes correctly", {
  M <- matrix(1:20, nrow = 10, ncol = 2)
  T_train <- 7
  
  result <- zscore_train(M, T_train)
  
  expect_type(result, "list")
  expect_named(result, c("Mz", "mu", "sd"))
  expect_equal(dim(result$Mz), dim(M))
  expect_length(result$mu, 2)
  expect_length(result$sd, 2)
  
  train_means <- colMeans(result$Mz[1:T_train, ])
  expect_true(all(abs(train_means) < 1e-10))
  
  train_sds <- apply(result$Mz[1:T_train, ], 2, sd)
  expect_true(all(abs(train_sds - 1) < 1e-10))
})

test_that("zscore_train handles zero variance", {
  M <- matrix(c(rep(5, 10), 1:10), nrow = 10, ncol = 2)
  T_train <- 7
  
  result <- zscore_train(M, T_train)
  
  expect_true(all(is.finite(result$Mz)))
  expect_equal(result$sd[1], 1)
})

test_that("null coalescing operator works", {
  expect_equal(NULL %||% 5, 5)
  expect_equal(3 %||% 5, 3)
  expect_equal("a" %||% "b", "a")
})

test_that("align_columns reorders correctly", {
  A <- matrix(1:6, nrow = 2, ncol = 3)
  colnames(A) <- c("x", "y", "z")
  
  B <- matrix(7:12, nrow = 2, ncol = 3)
  colnames(B) <- c("z", "x", "y")
  
  result <- align_columns(A, B)
  
  expect_equal(colnames(result), colnames(A))
})

test_that("align_columns errors on missing columns", {
  A <- matrix(1:6, nrow = 2, ncol = 3)
  colnames(A) <- c("x", "y", "z")
  
  B <- matrix(1:4, nrow = 2, ncol = 2)
  colnames(B) <- c("x", "y")
  
  expect_error(align_columns(A, B), "Missing columns")
})

test_that("check_stan_backend returns valid values", {
  result <- check_stan_backend(verbose = FALSE)
  expect_true(result %in% c("cmdstanr", "rstan", "none"))
})

test_that("vmsg respects verbose flag", {
  expect_silent(vmsg("test", verbose = FALSE))
  expect_message(vmsg("test message", verbose = TRUE), "test message")
})
