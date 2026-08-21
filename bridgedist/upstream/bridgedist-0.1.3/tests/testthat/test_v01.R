library(bridgedist)

test_that("Confirm unit variance +/- 0.1 for phi = 1/sqrt(1+3/pi^2)", {
  expect_equal(object = var(rbridge(1e5, phi = 1/sqrt(1+3/pi^2))),
               expected = 1,
               tolerance = 0.1,
               scale = 1)
})

test_that("Confirm p(q(x)) = x", {
  phi <- runif(1, 0.01, 0.99)
  x <- runif(1, 0.01, 0.99)
  expect_equal(object = pbridge(qbridge(x, phi=phi), phi=phi),
               expected = x,
               tolerance = 1e-6,
               scale = 1)
})

test_that("Confirm q(p(x)) = x", {
  phi <- runif(1, 0.01, 0.99)
  x <- rnorm(1, sd=2)
  expect_equal(object = qbridge(pbridge(x, phi=phi), phi=phi),
               expected = x,
               tolerance = 1e-6,
               scale = 1)
})


test_that("Confirm log = TRUE feature works in dbridge", {
  phi <- runif(1, 0.01, 0.99)
  x <- rnorm(1, sd=2)
  expect_equal(object =       dbridge(x, phi=phi, log=TRUE),
               expected = log(dbridge(x, phi=phi, log=FALSE)),
               tolerance = 1e-6,
               scale = 1)
})


test_that("Confirm log.p = TRUE feature works in pbridge", {
  phi <- runif(1, 0.01, 0.99)
  x <- rnorm(1, sd=2)
  expect_equal(object =       pbridge(x, phi=phi, log.p=TRUE),
               expected = log(pbridge(x, phi=phi, log.p=FALSE)),
               tolerance = 1e-6,
               scale = 1)
})

test_that("Confirm log.p = TRUE feature works in qbridge", {
  phi <- runif(1, 0.01, 0.99)
  x <- runif(1, 0.01, 0.99)
  expect_equal(object =       qbridge(log(x), phi=phi, log.p=TRUE),
               expected =     qbridge(    x , phi=phi, log.p=FALSE),
               tolerance = 1e-6,
               scale = 1)
})


test_that("Confirm lower.tail = FALSE feature works in pbridge", {
  phi <- runif(1, 0.01, 0.99)
  x <- rnorm(1, sd=2)
  expect_equal(object =       pbridge(x, phi=phi, lower.tail=TRUE),
               expected =   1-pbridge(x, phi=phi, lower.tail=FALSE),
               tolerance = 1e-6,
               scale = 1)
})

test_that("Confirm lower.tail = FALSE feature works in qbridge", {
  phi <- runif(1, 0.01, 0.99)
  x <- runif(1, 0.01, 0.99)
  expect_equal(object =       qbridge(  x, phi=phi, lower.tail=TRUE),
               expected =     qbridge(1-x, phi=phi, lower.tail=FALSE),
               tolerance = 1e-6,
               scale = 1)
})


test_that("Confirm lower.tail = FALSE feature works in pbridge with log.p=T", {
  phi <- runif(1, 0.01, 0.99)
  x <- rnorm(1, sd=2)
  expect_equal(object =             pbridge(x, phi=phi, lower.tail=TRUE , log.p=TRUE),
               expected = log(1-exp(pbridge(x, phi=phi, lower.tail=FALSE, log.p=TRUE))),
               tolerance = 1e-6,
               scale = 1)
})

test_that("Confirm lower.tail = FALSE feature works in qbridge with log.p=T", {
  phi <- runif(1, 0.01, 0.99)
  x <- runif(1, 0.01, 0.99)
  expect_equal(object =       qbridge(log(  x), phi=phi, lower.tail=TRUE , log.p=TRUE),
               expected =     qbridge(log(1-x), phi=phi, lower.tail=FALSE, log.p=TRUE),
               tolerance = 1e-6,
               scale = 1)
})



test_that("Confirm log.p = TRUE feature works in pbridge with lower.tail=F", {
  phi <- runif(1, 0.01, 0.99)
  x <- rnorm(1, sd=2)
  expect_equal(object =       pbridge(x, phi=phi, lower.tail=FALSE, log.p=TRUE),
               expected = log(pbridge(x, phi=phi, lower.tail=FALSE, log.p=FALSE)),
               tolerance = 1e-6,
               scale = 1)
})

test_that("Confirm log.p = TRUE feature works in qbridge with lower.tail=F", {
  phi <- runif(1, 0.01, 0.99)
  x <- runif(1, 0.01, 0.99)
  expect_equal(object =       qbridge(log(x), phi=phi, lower.tail=FALSE, log.p=TRUE),
               expected =     qbridge(    x , phi=phi, lower.tail=FALSE, log.p=FALSE),
               tolerance = 1e-6,
               scale = 1)
})

test_that("test first logical", {
  phi <- runif(1, 0.01, 0.99)
  x <- runif(1, 0.01, 0.99)
  expect_equal(object =       qbridge(log(x), phi=phi, lower.tail=FALSE, log.p=c(TRUE ,FALSE)),
               expected =     qbridge(    x , phi=phi, lower.tail=FALSE, log.p=c(FALSE,TRUE)),
               tolerance = 1e-6,
               scale = 1)
})


test_that("rbridgeRecycle n=1", {
  n <- 1
  phi <- runif(sample(2:10),0.01,0.99)
  expect_equal(object =   length(rbridge(n,phi))    ,
               expected =   1  ,
               tolerance = 1e-6,
               scale = 1)
})


test_that("rbridgeRecycle n>1 length(phi) > n", {
  n <- 5
  phi <- rep(c(0.01,0.99), length=6)
  expect_equal(object =   length(rbridge(n,phi))    ,
               expected =   5  ,
               tolerance = 1e-6,
               scale = 1)
})


test_that("rbridgeRecycle length(n)>1 length(phi) > length(n)", {
  n <- rep(NA,5)
  phi <- rep(c(0.01,0.99), length=6)
  expect_equal(object =   length(rbridge(n,phi))    ,
               expected =   5  ,
               tolerance = 1e-6,
               scale = 1)
})

