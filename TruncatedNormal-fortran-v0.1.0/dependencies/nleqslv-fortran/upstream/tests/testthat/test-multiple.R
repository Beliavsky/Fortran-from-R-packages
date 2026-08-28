# Copyright (c) Rob Carnell 2026

test_that("multiple functions with nleqslv", {
  func_list <- list(brown, brown, dcbval, dciequ, helval, pwlsng, rosbrk,
                    vardim, watson, watson, wood, pwlbsc, dslnex,
                    f_nocedal, sinexp, sinexp, troesch, kearfott,
                    hdp, twoip, broyden_f1, broyden_f1,
                    broyden_f2, broyden_f2, broyden_f2, broyden_f2,
                    broyden_f7, boggs, octxmpl)
  start_list <- list(brown_xstart(10), brown_xstart(20), dcbval_xstart(10),
                     dciequ_xstart(10), helval_xstart, pwlsng_xstart,
                     rosbrk_xstart, vardiminit(10), watson_xstart(6),
                     watson_xstart(9), wood_xstart, pwlbsc_xstart,
                     dslnex_xstart, nocedal_xstart,
                     sinexp_xstart(0.25), sinexp_xstart(0.75),
                     troesch_xstart(50), kearfott_xstart,
                     hdp_xstart, twoip_xstart, c(0, 0), c(1, 1),
                     c(1, 1, 0), c(-1, 1, 0), c(1, -1, 0), c(-1, -1, 0),
                     broyden_f7_xstart, boggs_xstart, octxmpl_xstart)

  expect_equal(length(func_list), length(start_list))

  for (i in seq_along(func_list)) {
    soln <- nleqslv(start_list[[i]], func_list[[i]], control = list(maxit=300))
    expect_true(soln$termcd %in% c(1, 4))
    expect_equal(soln$message, expectedMessage1)
    expect_true(all(abs(soln$fvec) <= 1e-6))
  }
})

test_that("multiple functions with testnslv", {
  func_list <- list(brown, dcbval, dciequ, helval, pwlsng, rosbrk,
                    vardim, watson, watson, wood, pwlbsc)
  start_list <- list(brown_xstart(20), dcbval_xstart(10), dciequ_xstart(10),
                     helval_xstart, pwlsng_xstart, rosbrk_xstart, vardiminit(10),
                     watson_xstart(6), watson_xstart(9), wood_xstart, pwlbsc_xstart)

  expect_equal(length(func_list), length(start_list))

  for (i in seq_along(func_list)) {
    soln <- testnslv(start_list[[i]], func_list[[i]])
    expect_true(inherits(soln, "test.nleqslv"))
    expect_true(is.data.frame(soln$out))
    expect_true(all(soln$out$termcd %in% c(1,4,5,6)))
    # expect_true(all(soln$out$Iter < 92)) # macos latest fails this test on 3/22/2026
  }
})

test_that("vardim", {
  soln <- nleqslv(vardiminit(10), vardim)
  expect_equal(vardimsol(10), soln$x)

  soln <- nleqslv(vardiminit(5), vardim)
  expect_equal(vardimsol(5), soln$x)
})

test_that("with jacobian", {
  func_list <- list(pwlbsc, chandraH, chandraH, broyden_f7, bkrust, li5,
                    neum1, pcex)
  start_list <- list(pwlbsc_xstart, chandra_xstart(50), chandra_xstart(100),
                     broyden_f7_xstart, bkrust_xstart, li5_xstart,
                     neum1_xstart, pcex_xstart)
  jac_list <- list(pwlbscjac, chandraH.jac, chandraH.jac, broyden_f7jac,
                   bkrjac, li5jac, neum1jac, pcexjac)

  expect_equal(length(func_list), length(start_list))
  expect_equal(length(func_list), length(jac_list))

  for (i in seq_along(func_list)) {
    if (identical(func_list[[i]], li5)) {
      soln <- nleqslv(start_list[[i]], func_list[[i]], jac_list[[i]], global = "gline")
    } else {
      soln <- nleqslv(start_list[[i]], func_list[[i]], jac_list[[i]])
    }
    expect_true(soln$termcd %in% c(1,3))
    expect_true(all(abs(soln$fvec) <= 1e-4))
  }
})

test_that("dslnex solutions", {
  z <- nleqslv(dslnex_xstart, dslnex, jacobian=TRUE)
  expect_equal(jacdsln(z$x), z$jac, tolerance = 0.1, ignore_attr = TRUE)

  z <- nleqslv(dslnex_xstart, dslnex, global="none", jacobian=TRUE)
  expect_equal(jacdsln(z$x), z$jac, tolerance = 0.1, ignore_attr = TRUE)

  z <- nleqslv(dslnex_xstart, dslnex, method="Newton", jacobian=TRUE)
  expect_equal(jacdsln(z$x), z$jac, tolerance = 0.1, ignore_attr = TRUE)

  z <- nleqslv(dslnex_xstart, dslnex, method="Newton", jacobian=TRUE, global="none")
  expect_equal(jacdsln(z$x), z$jac, tolerance = 0.1, ignore_attr = TRUE)

  dslnex_xstart <- c(2, .5)

  z <- nleqslv(dslnex_xstart, dslnex, jacobian=TRUE)
  expect_equal(jacdsln(z$x), z$jac, tolerance = 0.1, ignore_attr = TRUE)

  z <- nleqslv(dslnex_xstart, dslnex, method="Newton", jacobian=TRUE)
  expect_equal(jacdsln(z$x), z$jac, tolerance = 0.1, ignore_attr = TRUE)

  set.seed(11)

  dslnex_xstart <- matrix(runif(50, -3, 3), ncol=2)
  temp <- searchZeros(dslnex_xstart, dslnex)
  # two solutions of two values
  expect_equal(nrow(temp$x), 2)
  expect_equal(ncol(temp$x), 2)
  expect_true(all(apply(temp$x, 1, dslnex) < 1E-9))

})

test_that("nocedal", {
  z <- nleqslv(nocedal_xstart, f_nocedal, global="none")
  expect_equal(z$fvec, c(0, 0), tolerance = 1E-9)
  z <- nleqslv(nocedal_xstart, f_nocedal, method="Newton")
  expect_equal(z$fvec, c(0, 0), tolerance = 1E-9)
  z <- nleqslv(nocedal_xstart, f_nocedal, jac_nocedal, method="Newton")
  expect_equal(z$fvec, c(0, 0), tolerance = 1E-9)
  expect_equal(jac_nocedal(c(0,1)), jac_nocedal(z$x), tolerance = 1E-9)
})

test_that("sinexp searchZeros", {
  N <- 5
  set.seed(11)
  xstart <- cbind(runif(N, min=0.25, max=1), runif(N, min=1.5, max=2*pi))

  temp <- searchZeros(xstart, sinexp, method="Broyden", global="dbldog")

  expect_equal(nrow(temp$x), 2)
  expect_equal(ncol(temp$x), 2)
  expect_true(all(abs(apply(temp$x, 1, sinexp)) < 1E-8))
})

test_that("circle", {
  act <- circle.intersect(xc1,xc2,r0,r1)
  temp1 <- nleqslv(c(0,0), fcircle, xc1=xc1, xc2=xc2, r0=r0, r1=r1)
  temp2 <- nleqslv(c(2,2), fcircle, xc1=xc1, xc2=xc2, r0=r0, r1=r1)

  expect_equal(act$P2, temp1$x)
  expect_equal(act$P1, temp2$x)

  xstart <- matrix(c(0,0,2,2),nrow=2,byrow=TRUE)
  sz <- searchZeros(xstart, fcircle, xc1=xc1, xc2=xc2, r0=r0, r1=r1)
  expect_true(all(abs(act$P1 - sz$x[1,]) < 1E-6) |
                all(abs(act$P1 - sz$x[2,]) < 1E-6))
  expect_true(all(abs(act$P2 - sz$x[1,]) < 1E-6) |
                all(abs(act$P2 - sz$x[2,]) < 1E-6))

})

test_that("cutlip", {
  # paper has wrong order of parameters
  # use the Fortran program to get the correct values

  # parameter set 1
  k1 <-  31.24
  k2 <-  2.062
  k3 <-  303.03
  kr1 <-  0.272
  kr2 <-  0.02

  # initial estimate and solutions page 1464
  xinit1 <- c(0.99, 0.05, 0.05, 0.99, 0.05, 0)
  xinit2 <- c(0.05, 0.99, 0.05, 0.05, 0.99, 0)
  xinit3 <- c(0.05, 0.5, 2, 0.05, 0.99, 0.5)

  # solution Table IV (parameter set; used here)
  xsolA <- c(0.97007,0.98005,0.05985,0.99003,0.0000997,0.009873)
  xsolB <- c(0.035614, 0.35708, 1.9288, 0.035997,0.088409, 0.87559)
  xsolC <- c(1.03329,1.02220,-0.066597,-0.000109,1.00114, -0.00103)

  temp1 <- nleqslv(xinit1, cutlip, k1=k1, k2=k2, k3=k3, kr1=kr1, kr2=kr2)
  expect_equal(temp1$x, xsolA, tolerance = 1E-5)
  temp2 <- nleqslv(xinit2, cutlip, k1=k1, k2=k2, k3=k3, kr1=kr1, kr2=kr2)
  expect_equal(temp2$x, xsolC, tolerance = 1E-5)
  temp3 <- nleqslv(xinit3, cutlip, k1=k1, k2=k2, k3=k3, kr1=kr1, kr2=kr2,
                   method="Broyden", global="dbldog")
  expect_equal(temp3$x, xsolB, tolerance = 1E-4)

  # parameter set 2
  k1 <-  17.721
  k2 <-  3.483
  k3 <-  505.051
  kr1 <- 0.118
  kr2 <- 0.033

  # initial estimate and solutions page 1464
  xinit1 <- c(0.99, 0.05, 0.05, 0.99, 0.05, 0)
  xinit2 <- c(0.05, 0.99, 0.05, 0.05, 0.99, 0)

  # solution Table V (parameter set 2; used here)

  xsolA <- c(0.94994,0.96662,0.10013,0.98999,0.0001001,0.009913)
  xsolB <- c(0.11766,0.41177,1.7647 ,0.003045,0.57362, 0.42334)
  xsolC <- c(1.0698,1.04655,-0.13967,-0.0001377,1.003822, -0.003684)

  temp1 <- nleqslv(xinit1, cutlip, k1=k1, k2=k2, k3=k3, kr1=kr1, kr2=kr2)
  expect_equal(temp1$x, xsolA, tolerance = 1E-5)
  temp2 <- nleqslv(xinit2, cutlip, k1=k1, k2=k2, k3=k3, kr1=kr1, kr2=kr2)
  expect_equal(temp2$x, xsolC, tolerance = 1E-4)

  # parameter set 3
  k1 <-  17.721
  k2 <-  6.996
  k3 <-  505.051
  kr1 <- 0.118
  kr2 <- 333.333

  # initial estimate and solutions page 1464
  xinit1 <- c(0.99, 0.05, 0.05, 0.99, 0.05, 0)

  # solution Table VI (parameter set; used here)
  xsolA <- c(0.94994,0.96662,0.10013,0.98999,0.0001001,0.009913)

  # xinit2 does not converge; iteration limit exceeded
  temp1 <- nleqslv(xinit1, cutlip, k1=k1, k2=k2, k3=k3, kr1=kr1, kr2=kr2)
  expect_equal(temp1$x, xsolA, tolerance = 1E-3)
})

test_that("kuno seader", {
  N <- 100
  xstart <- matrix(runif(2*N, min=-5, max=5), N, 2)

  z <- searchZeros(xstart, kunoseader_f)
  expect_equal(nrow(z$x), 7)

  Nrep <- 90
  Ncol <- 3
  xstart <- matrix(runif(Ncol*Nrep, min=0, max=1), Nrep, Ncol)

  z <- searchZeros(xstart, kunoseader_f34)
  expect_equal(nrow(z$x), 4)
})

test_that("econ with const", {
  const <- rep(0, 20)
  xstart <- rep(1, 20) / 20

  temp <- nleqslv(xstart, econ, const = const)
  expect_true(all(temp$fvec < 1E-8))
  test_temp <- testnslv(xstart, econ, const = const,
                        global = c("cline", "qline", "pwldob", "dbldog"))
  expect_true(all(test_temp$out$Fnorm < 1E-12))
})

test_that("hiebert", {
  xstart <- c(1, 1, 10, 1, 1, 1, 0, 0, 0, 0)

  temp <- nleqslv(xstart, hiebert, R=10)
  expect_true(all(temp$fvec < 1E-8))
  test_temp <- testnslv(xstart, hiebert, R=10, method = "Newton")
  expect_true(all(test_temp$out$Fnorm < 1E-12))
})

test_that("pce", {
  soln <- nleqslv(pcex_xstart, pcex, pcexjac, method="Broyden", global="hook",
                  control=list(trace=0, delta="cauchy", ftol=1e-6, xtol=1e-5))
  expect_equal(soln$x, pcex_xsol, tolerance = 1E-3)
})

test_that("frac", {
  for( k in 1:7 ) {
    a <- k/10
    xstart <- c(a,1-a, max(0.4-0.5*a,0.1))
    z <- nleqslv(xstart, frac2, control=list(stepmax=.1))
    expect_true(all(abs(z$fvec) < 1E-8))
  }
})


