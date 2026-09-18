test_that("int returns integer values in [min,max]", {
  rng <- randompack_rng("pcg64")
  x <- rng$int(10, -2, 3)
  expect_type(x, "integer")
  expect_length(x, 10)
  expect_true(all(x >= -2))
  expect_true(all(x <= 3))
})

test_that("perm returns a permutation of 1..n", {
  rng <- randompack_rng("pcg64")
  x <- rng$perm(10)
  expect_type(x, "integer")
  expect_length(x, 10)
  expect_identical(sort(x), 1:10)
})

test_that("sample returns unique values in 1..n", {
  rng <- randompack_rng("pcg64")
  x <- rng$sample(10, 5)
  expect_type(x, "integer")
  expect_length(x, 5)
  expect_identical(length(unique(x)), 5L)
  expect_true(all(x >= 1))
  expect_true(all(x <= 10))
})

test_that("raw returns a raw vector", {
  rng <- randompack_rng("pcg64")
  x <- rng$raw(16)
  expect_type(x, "raw")
  expect_length(x, 16)
})
