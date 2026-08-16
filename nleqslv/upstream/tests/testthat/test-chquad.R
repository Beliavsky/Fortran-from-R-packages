# Copyright (c) Rob Carnell 2026

# More, Garbow, Hillstrom: Testing Unconstrained Optimization Software
# ACM Trans. Math. Software, 7, March 1981, 17--41
# Function as shown in this paper is for optimizing not for non linear equations
# from the paper it is not clear what was used by the authors for their tests
# of non linear equation solvers
# function 35

# slightly modified from original
# this version taken from Julia's package NLsolve test
# Burkardt (https://people.sc.fsu.edu/~jburkardt/f_src/test_nls/test_nls.f90) problem 15
# http://ftp.mcs.anl.gov/pub/MINPACK-2/tprobs.92/dchqfj.f

chebyquad <- function(x) {
  n <- length(x)
  y <- numeric(n)

  for(j in 1:n) {
    t1 <- 1.0
    t2 <- 2.0*x[j] - 1.0
    temp <- 2.0 * t2

    for(i in 1:n) {
      y[i] <- y[i] + t2

      t3 <- temp * t2 - t1
      t1 <- t2
      t2 <- t3
    }
  }

  y <- y / n

  for(i in 1:n) {
    if ( i%%2 == 0 ) {
      y[i] <- y[i] + 1.0 / (i * i - 1)
    }
  }

  y
}

# from dchqfj.f

chebyjac <- function(x) {
  #         do 90 j = 1, n
  #            temp1 = one
  #            temp2 = two*x(j) - one
  #            temp = two*temp2
  #            temp3 = zero
  #            temp4 = two
  #            do 80 i = 1, m
  #               fjac(i,j) = dx*temp4
  #               ti = four*temp2 + temp*temp4 - temp3
  #               temp3 = temp4
  #               temp4 = ti
  #               ti = temp*temp2 - temp1
  #               temp1 = temp2
  #               temp2 = ti
  #   80       continue
  #   90    continue

  n <- length(x)
  y <- numeric(n)
  dx <- 1/n
  jac <- matrix(0,nrow=n,ncol=n)
  for( j in 1:n) {
    t1 <- 1.0
    t2 <- 2.0*x[j] - 1.0
    temp <- 2.0 * t2
    t3 <- 0.0
    t4 <- 2.0
    for( i in 1:n) {
      jac[i,j] <- dx*t4
      ti <- 4.0 * t2 + temp * t4 - t3
      t3 <- t4
      t4 <- ti
      ti <- temp*t2 - t1
      t1 <- t2
      t2 <- ti
    }
  }
  jac
}

chebyinit <- function(n) (1:n) / (n + 1)

test_that("Chebyquad solutions converge for n = 1:7,9 using default method", {
  for (n in c(1:7, 9)) {
    xstart <- chebyinit(n)

    zz <- nleqslv(
      xstart, chebyquad,
      global = "dbldog",
      control = list(
        ftol = 1e-8,
        xtol = 1e-8,
        trace = 0,
        btol = 0.01,
        delta = -2
      )
    )

    expect_true(all(abs(zz$fvec) <= 1e-8))
  }
})


test_that("Chebyquad solutions converge for n = 1:7,9 using Newton method", {
  for (n in c(1:7, 9)) {
    xstart <- chebyinit(n)

    zz <- nleqslv(
      xstart, chebyquad,
      global = "dbldog",
      method = "Newton",
      control = list(
        ftol = 1e-8,
        xtol = 1e-8,
        trace = 0,
        btol = 0.01,
        delta = -2
      )
    )

    expect_true(all(abs(zz$fvec) <= 1e-8))
  }

  for (n in c(1:7, 9)) {
    xstart <- chebyinit(n)

    zz <- nleqslv(
      xstart, chebyquad,
      chebyjac,
      method = "Newton",
      control = list(
        ftol = 1e-8,
        xtol = 1e-8,
        trace = 0,
        btol = 0.01,
        delta = -2,
        chkjac=TRUE
      )
    )

    expect_true(all(abs(zz$fvec) <= 1e-8))
  }
})

test_that("chebyquad solutions with testnslv", {
  xstart <- chebyinit(5)

  temp <- testnslv(xstart, chebyquad)
  expect_true(inherits(temp, "test.nleqslv"))
  expect_true(is.data.frame(temp$out))
  expect_true(all(temp$out$termcd %in% c(1,2)))
  expect_true(all(temp$out$Iter < 11))

  temp <- testnslv(xstart, chebyquad, chebyjac)
  expect_true(inherits(temp, "test.nleqslv"))
  expect_true(is.data.frame(temp$out))
  expect_true(all(temp$out$termcd %in% c(1,2)))
  expect_true(all(temp$out$Iter < 10))
})
