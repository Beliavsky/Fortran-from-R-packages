test_that("randompack_rng creates a working RNG", {
  rng <- randompack_rng("pcg64")
  expect_true(is.list(rng) || is.environment(rng))
})

test_that("engines returns a data frame", {
  engines <- randompack_engines()
  expect_true(is.data.frame(engines))
  expect_true(all(c("engine", "description") %in% names(engines)))
})

test_that("serialize and deserialize restore RNG state", {
  rng1 <- randompack_rng("pcg64")
  a <- rng1$unif(10)
  state <- rng1$serialize()
  b <- rng1$unif(10)
  rng2 <- randompack_rng("pcg64")
  rng2$deserialize(state)
  b2 <- rng2$unif(10)
  expect_identical(b, b2)
})

test_that("duplicate creates independent RNG with identical initial state", {
  rng1 <- randompack_rng("pcg64")
  rng2 <- rng1$duplicate()
  expect_identical(rng1$unif(5), rng2$unif(5))
  rng1$unif(1)
  expect_false(identical(rng1$unif(5), rng2$unif(5)))
})
