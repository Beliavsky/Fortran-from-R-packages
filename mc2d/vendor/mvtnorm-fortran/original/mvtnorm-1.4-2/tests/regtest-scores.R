
library("mvtnorm")

set.seed(29081975)

chk <- function(...) 
    stopifnot(all.equal(..., tol = 1e-5, check.attributes = FALSE))

EVAL <- function(...) {}

if (require("numDeriv", quietly = TRUE) && 
    require("qrng", quietly = TRUE))
  EVAL <- eval

N <- 10
M <- 10000
MM <- M / N

prb <- 1:3 / 4

### chol 
thischeck <- expression({
  J <- cJ + dJ
  W <- NULL
  if (dJ > 1)
    W <- t(ghalton(M, d = dJ - 1))

  prm <- matrix(runif(J * (J - 1) / 2), ncol = 1)
  C <- ltMatrices(prm, byrow = BYROW)
  Z <- matrix(rnorm(J * N), ncol = N)
  Y <- Mult(C, Z)
  obs <- NULL
  if (cJ)
      obs <- Y[1:cJ,,drop = FALSE]
  lwr <- upr <- NULL
  if (dJ) {
    lwr <- t(apply(Y[cJ + (1:dJ),,drop = FALSE], 1, function(y) {
      qy <- quantile(y, prob = prb)
      c(-Inf, qy)[cut(y, breaks = c(-Inf, qy, Inf))]
    }))
    upr <- t(apply(Y[cJ + (1:dJ),,drop = FALSE], 1, function(y) {
      qy <- quantile(y, prob = prb)
      c(qy, Inf)[cut(y, breaks = c(-Inf, qy, Inf))]
    }))
  }

  ll <- function(prm) {
    C <- ltMatrices(prm, byrow = BYROW)
    -ldpmvnorm(obs = obs, lower = lwr, upper = upr, chol = C, M = MM, w = W)
  }

  sc <- function(prm) {
    C <- ltMatrices(prm, byrow = BYROW)
    ret <- sldpmvnorm(obs = obs, lower = lwr, upper = upr, chol = C,
                      M = MM, w = W)$chol
    -rowSums(Lower_tri(ret))
  }

  theta <- runif(J * (J - 1) / 2)
  print(ll(theta))
  chk(grad(ll, theta), sc(theta))
})

BYROW <- FALSE
cJ <- 4
dJ <- 4
EVAL(thischeck)

cJ <- 1
dJ <- 4
EVAL(thischeck)

cJ <- 4
dJ <- 1
EVAL(thischeck)

cJ <- 0
dJ <- 4
EVAL(thischeck)

cJ <- 4
dJ <- 0
EVAL(thischeck)


BYROW <- TRUE
cJ <- 4
dJ <- 4
EVAL(thischeck)

cJ <- 1
dJ <- 4
EVAL(thischeck)

cJ <- 4
dJ <- 1
EVAL(thischeck)

cJ <- 0
dJ <- 4
EVAL(thischeck)

cJ <- 4
dJ <- 0
EVAL(thischeck)

### invchol
thischeck <- expression({
  J <- cJ + dJ
  W <- NULL
  if (dJ > 1)
    W <- t(ghalton(M, d = dJ - 1))

  prm <- matrix(runif(J * (J - 1) / 2), ncol = 1)
  C <- ltMatrices(prm, byrow = BYROW)
  Z <- matrix(rnorm(J * N), ncol = N)
  Y <- Mult(C, Z)
  obs <- NULL
  if (cJ)
      obs <- Y[1:cJ,,drop = FALSE]
  lwr <- upr <- NULL
  if (dJ) {
    lwr <- t(apply(Y[cJ + (1:dJ),,drop = FALSE], 1, function(y) {
      qy <- quantile(y, prob = prb)
      c(-Inf, qy)[cut(y, breaks = c(-Inf, qy, Inf))]
    }))
    upr <- t(apply(Y[cJ + (1:dJ),,drop = FALSE], 1, function(y) {
      qy <- quantile(y, prob = prb)
      c(qy, Inf)[cut(y, breaks = c(-Inf, qy, Inf))]
    }))
  }

  ll <- function(prm) {
    L <- ltMatrices(prm, byrow = BYROW)
    -ldpmvnorm(obs = obs, 
              lower = lwr, upper = upr, invchol = L, M = MM, w = W)
  }

  sc <- function(prm) {
    L <- ltMatrices(prm, byrow = BYROW)
    ret <- sldpmvnorm(obs = obs, 
                      lower = lwr, upper = upr, invchol = L,
                      M = MM, w = W)$invchol
    -rowSums(Lower_tri(ret))
  }

  theta <- runif(J * (J - 1) / 2)
  C <- ltMatrices(matrix(theta, ncol = 1), byrow = BYROW)
  theta <- Lower_tri(solve(C))
  print(ll(theta))
  chk(grad(ll, theta), sc(theta))
})


BYROW <- FALSE
cJ <- 4
dJ <- 4
EVAL(thischeck)

cJ <- 1
dJ <- 4
EVAL(thischeck)

cJ <- 4
dJ <- 1
EVAL(thischeck)

cJ <- 0
dJ <- 4
EVAL(thischeck)

cJ <- 4
dJ <- 0
EVAL(thischeck)


BYROW <- TRUE
cJ <- 4
dJ <- 4
EVAL(thischeck)

cJ <- 1
dJ <- 4
EVAL(thischeck)

cJ <- 4
dJ <- 1
EVAL(thischeck)

cJ <- 0
dJ <- 4
EVAL(thischeck)

cJ <- 4
dJ <- 0
EVAL(thischeck)

### chol standardized
thischeck <- expression({
  J <- cJ + dJ
  W <- NULL
  if (dJ > 1)
    W <- t(ghalton(M, d = dJ - 1))

  prm <- matrix(runif(J * (J - 1) / 2), ncol = 1)
  C <- ltMatrices(prm)
  C <- ltMatrices(C, byrow = BYROW)
  Z <- matrix(rnorm(J * N), ncol = N)
  Y <- Mult(C, Z)
  obs <- NULL
  if (cJ)
      obs <- Y[1:cJ,,drop = FALSE]
  lwr <- upr <- NULL
  if (dJ) {
    lwr <- t(apply(Y[cJ + (1:dJ),,drop = FALSE], 1, function(y) {
      qy <- quantile(y, prob = prb)
      c(-Inf, qy)[cut(y, breaks = c(-Inf, qy, Inf))]
    }))
    upr <- t(apply(Y[cJ + (1:dJ),,drop = FALSE], 1, function(y) {
      qy <- quantile(y, prob = prb)
      c(qy, Inf)[cut(y, breaks = c(-Inf, qy, Inf))]
    }))
  }

  ll <- function(prm) {
    C <- ltMatrices(prm, byrow = BYROW)
    Cs <- standardize(chol = C)
    -ldpmvnorm(obs = obs, lower = lwr, upper = upr, chol = Cs, M = MM, w = W)
  }

  sc <- function(prm) {
    C <- ltMatrices(prm, byrow = BYROW)
    Cs <- standardize(chol = C)
    ret <- sldpmvnorm(obs = obs, lower = lwr, upper = upr, chol = Cs,
                      M = MM, w = W)$chol
    ret <- destandardize(chol = C, score_schol = ret)
    -rowSums(Lower_tri(ret))
  }

  theta <- runif(J * (J - 1) / 2)
  print(ll(theta))
  chk(grad(ll, theta), sc(theta))
})

BYROW <- FALSE
cJ <- 4
dJ <- 4
EVAL(thischeck)

cJ <- 1
dJ <- 4
EVAL(thischeck)

cJ <- 4
dJ <- 1
EVAL(thischeck)

cJ <- 0
dJ <- 4
EVAL(thischeck)

cJ <- 4
dJ <- 0
EVAL(thischeck)


BYROW <- TRUE
cJ <- 4
dJ <- 4
EVAL(thischeck)

cJ <- 1
dJ <- 4
EVAL(thischeck)

cJ <- 4
dJ <- 1
EVAL(thischeck)

cJ <- 0
dJ <- 4
EVAL(thischeck)

cJ <- 4
dJ <- 0
EVAL(thischeck)

### invchol standardized
thischeck <- expression({
  J <- cJ + dJ
  W <- NULL
  if (dJ > 1)
    W <- t(ghalton(M, d = dJ - 1))

  prm <- matrix(runif(J * (J - 1) / 2), ncol = 1)
  C <- ltMatrices(prm, byrow = BYROW)
  Z <- matrix(rnorm(J * N), ncol = N)
  Y <- Mult(C, Z)
  obs <- NULL
  if (cJ)
      obs <- Y[1:cJ,,drop = FALSE]
  lwr <- upr <- NULL
  if (dJ) {
    lwr <- t(apply(Y[cJ + (1:dJ),,drop = FALSE], 1, function(y) {
      qy <- quantile(y, prob = prb)
      c(-Inf, qy)[cut(y, breaks = c(-Inf, qy, Inf))]
    }))
    upr <- t(apply(Y[cJ + (1:dJ),,drop = FALSE], 1, function(y) {
      qy <- quantile(y, prob = prb)
      c(qy, Inf)[cut(y, breaks = c(-Inf, qy, Inf))]
    }))
  }

  ll <- function(prm) {
    L <- ltMatrices(prm, byrow = BYROW)
    Ls <- standardize(invchol = L)
    -ldpmvnorm(obs = obs, 
              lower = lwr, upper = upr, invchol = Ls, M = MM, w = W)
  }

  sc <- function(prm) {
    L <- ltMatrices(prm, byrow = BYROW)
    Cs <- standardize(chol = solve(L))
    ret <- sldpmvnorm(obs = obs, 
                      lower = lwr, upper = upr, chol = Cs,
                      M = MM, w = W)$chol
    ret <- destandardize(invchol = L, score_schol = ret)
    -rowSums(Lower_tri(ret))
  }

  theta <- runif(J * (J - 1) / 2)
  C <- ltMatrices(matrix(theta, ncol = 1), byrow = BYROW)
  theta <- Lower_tri(solve(C))
  print(ll(theta))
  chk(grad(ll, theta), sc(theta))
})


BYROW <- FALSE
cJ <- 4
dJ <- 4
EVAL(thischeck)

cJ <- 1
dJ <- 4
EVAL(thischeck)

cJ <- 4
dJ <- 1
EVAL(thischeck)

cJ <- 0
dJ <- 4
EVAL(thischeck)

cJ <- 4
dJ <- 0
EVAL(thischeck)


BYROW <- TRUE
cJ <- 4
dJ <- 4
EVAL(thischeck)

cJ <- 1
dJ <- 4
EVAL(thischeck)

cJ <- 4
dJ <- 1
EVAL(thischeck)

cJ <- 0
dJ <- 4
EVAL(thischeck)

cJ <- 4
dJ <- 0
EVAL(thischeck)
