test_that("seed resets RNG deterministically", {
  rng <- randompack_rng("pcg64")
  rng$seed(123)
  a <- rng$unif(5)
  rng$seed(123)
  b <- rng$unif(5)
  expect_identical(a, b)
})

test_that("randomize runs and allows draws", {
  rng <- randompack_rng("pcg64")
  rng$randomize()
  x <- rng$unif(5)
  expect_type(x, "double")
  expect_length(x, 5)
  expect_true(all(is.finite(x)))
})

test_that("full_mantissa create runs and allows draws", {
  rng <- randompack_rng("pcg64", full_mantissa = TRUE)
  x <- rng$unif(5)
  expect_type(x, "double")
  expect_length(x, 5)
  expect_true(all(is.finite(x)))
})

test_that("set_state accepts numeric state words", {
  rng1 <- randompack_rng("pcg64")
  rng2 <- randompack_rng("pcg64")
  state <- c(1, 0, 2, 0, 3, 0, 4, 0)
  rng1$set_state(state)
  rng2$set_state(state)
  expect_identical(rng1$unif(5), rng2$unif(5))
})

test_that("advance matches jump for pcg64", {
  rng1 <- randompack_rng("pcg64")
  rng2 <- randompack_rng("pcg64")
  state <- c(
    2887678797, 1542324465,
    850656868, 2020473006,
    2185206843, 651939783,
    3901109050, 223234302
  )
  rng1$set_state(state)
  rng2$set_state(state)
  rng1$advance(c(0, 0, 65536, 0))
  rng2$jump(80)
  expect_identical(rng1$unif(5), rng2$unif(5))
  rng1$set_state(state)
  rng2$set_state(state)
  rng1$advance(0)
  expect_identical(rng1$unif(5), rng2$unif(5))
  expect_error(rng1$advance(c(1, 2, 3, 4, 5)))
  expect_error(randompack_rng("squares")$advance(c(1, 0)))
})

test_that("pcg64_set_inc makes draws repeatable", {
  rng1 <- randompack_rng("pcg64")
  rng2 <- randompack_rng("pcg64")
  state <- c(1, 0, 2, 0, 1, 0, 0, 0)
  inc <- c(3, 0, 5, 0)
  rng1$set_state(state)
  rng2$set_state(state)
  rng1$pcg64_set_inc(inc)
  rng2$pcg64_set_inc(inc)
  expect_identical(rng1$unif(5), rng2$unif(5))
  expect_error(rng1$pcg64_set_inc(c(2, 0, 5, 0)))
})

test_that("cwg128_set_weyl makes draws repeatable", {
  skip_if_not("cwg128" %in% randompack_engines()$engine)
  rng1 <- randompack_rng("cwg128")
  rng2 <- randompack_rng("cwg128")
  state <- c(1, 0, 0, 0, 7, 0, 0, 0, 11, 0, 0, 0, 13, 0, 0, 0)
  rng1$set_state(state)
  rng2$set_state(state)
  rng1$cwg128_set_weyl(c(3, 0, 5, 0))
  rng2$cwg128_set_weyl(c(3, 0, 5, 0))
  expect_equal(rng1$raw(100), rng2$raw(100))
  expect_error(rng1$cwg128_set_weyl(c(2, 0, 5, 0)))
})

test_that("sfc64_set_abc makes draws repeatable", {
  rng1 <- randompack_rng("sfc64")
  rng2 <- randompack_rng("sfc64")
  abc <- c(7, 0, 11, 0, 13, 0)
  rng1$set_state(c(1, 0, 2, 0, 3, 0, 17, 0))
  rng2$set_state(c(1, 0, 2, 0, 3, 0, 17, 0))
  rng1$sfc64_set_abc(abc)
  rng2$sfc64_set_abc(abc)
  expect_identical(rng1$unif(5), rng2$unif(5))
})

test_that("chacha_set_nonce makes draws repeatable", {
  rng1 <- randompack_rng("chacha20")
  rng2 <- randompack_rng("chacha20")
  nonce <- c(7, 11, 13)
  rng1$seed(123)
  rng2$seed(123)
  rng1$chacha_set_nonce(nonce)
  rng2$chacha_set_nonce(nonce)
  expect_identical(rng1$unif(5), rng2$unif(5))
})

test_that("philox_set_key makes draws repeatable", {
  rng1 <- randompack_rng("philox")
  rng2 <- randompack_rng("philox")
  key <- c(5, 0, 6, 0)
  rng1$set_state(c(1, 0, 2, 0, 3, 0, 4, 0, 0, 0, 0, 0))
  rng1$philox_set_key(key)
  rng2$set_state(c(1, 0, 2, 0, 3, 0, 4, 0, 0, 0, 0, 0))
  rng2$philox_set_key(key)
  expect_identical(rng1$unif(5), rng2$unif(5))
})

test_that("squares_set_key makes draws repeatable", {
  rng1 <- randompack_rng("squares")
  rng2 <- randompack_rng("squares")
  key <- c(2, 0)
  rng1$set_state(c(1, 0, 0, 0))
  rng1$squares_set_key(key)
  rng2$set_state(c(1, 0, 0, 0))
  rng2$squares_set_key(key)
  expect_identical(rng1$unif(5), rng2$unif(5))
})

test_that("jump advances the default engine", {
  rng1 <- randompack_rng()
  rng2 <- randompack_rng()
  rng1$seed(123)
  rng2$seed(123)
  a <- rng1$unif(5)
  rng2$jump(128)
  b <- rng2$unif(5)
  expect_false(identical(a, b))
  expect_error(rng1$jump(129))
  expect_error(rng1$jump(254))
})
